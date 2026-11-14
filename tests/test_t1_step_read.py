"""experiments/t1_m1/step_read.py (try T44 of experiments/harness/TRIES.md) on synthetic data: the scaled arms and the
reading. No released file, statistics file or score is read."""

from __future__ import annotations

import importlib.util
import json
import sys

import numpy as np
import pytest
from test_g_mom_verify import _small_inputs
from test_t1_state_transport import _gmom_ctx

from vec import harness, t1_scope
from vec.boards import REPO_ROOT

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "step_read.py"
GROUPS = {"V-CM": "CM", "Endocardium": "Endocardium", "Blood": "Blood", "Foregut": "other:Foregut"}


@pytest.fixture
def sr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("step_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def test_the_arms_fractions_are_one_and_a_half_and_two_times_the_base(sr):
    assert sr.MULTIPLIERS == {"m1_r_gmom_x15": 1.5, "m1_r_gmom_x2": 2.0}
    fr = sr.arm_fractions(0.2)
    assert fr["m1_r_gmom"] == 0.2
    assert fr["m1_r_gmom_x15"] == 0.2 * 1.5 and fr["m1_r_gmom_x2"] == 0.2 * 2
    assert fr["m1_r_gmom_x15"] == pytest.approx(0.3) and fr["m1_r_gmom_x2"] == pytest.approx(0.4)


def test_an_arm_differs_from_the_base_build_in_the_fraction_only(sr, monkeypatch):
    calls = []

    def spy(ctx, directions, groups, fraction, held_m1, held_factors):
        calls.append({"ctx": ctx, "directions": directions, "groups": groups, "fraction": fraction,
                      "held_m1": held_m1, "held_factors": held_factors})
        return {"m1_r": ("m1_r", fraction), "m1_r_gmom": ("m1_r_gmom", fraction)}

    monkeypatch.setattr(sr.gc, "m1_r_arms", spy)
    ctx, dirs, groups, held_m1, held_f = object(), {"CM": np.zeros(3)}, {"V-CM": "CM"}, ["a", "b"], ["a"]
    arms = sr.step_arms(ctx, dirs, groups, 0.2, held_m1, held_f)
    assert set(arms) == {"m1_r_gmom", "m1_r_gmom_x15", "m1_r_gmom_x2"}
    assert len(calls) == 3
    base = calls[0]
    assert base["fraction"] == 0.2
    assert [c["fraction"] for c in calls] == [0.2, 0.2 * 1.5, 0.2 * 2]
    assert [c["fraction"] / base["fraction"] for c in calls] == pytest.approx([1.0, 1.5, 2.0], abs=1e-15)
    for c in calls:
        for key in ("ctx", "directions", "groups", "held_m1", "held_factors"):
            assert c[key] is base[key], f"an arm differs from the base build in {key}"
    assert arms["m1_r_gmom"] == ("m1_r_gmom", 0.2) and arms["m1_r_gmom_x2"] == ("m1_r_gmom", 0.4)


def _arms_world():
    harness._G_MOM_PREP.clear()
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start)
    groups = {lab: GROUPS.get(lab, "x") for lab in ref.obs[t1_scope.LABEL_KEY].astype(str).unique()}
    rng = np.random.default_rng(5)
    dirs = {g: rng.normal(0, 0.4, ref.n_vars) for g in ("CM", "Endocardium")}
    return ctx, groups, dirs, list(ctx.derived["t1_held_flat"])


def test_the_base_arm_is_t42s_m1_r_gmom_and_a_multiplier_of_one_is_the_base_bit_for_bit(sr):
    """The idle case: at a multiplier of 1 the arm is the base build bit for bit; the base is gmom_combo's
    m1_r_gmom at the base fraction; the real multipliers move the build."""
    ctx, groups, dirs, held = _arms_world()
    arms = sr.step_arms(ctx, dirs, groups, 0.2, held, held, multipliers={"idle": 1.0, **sr.MULTIPLIERS})
    direct = sr.gc.m1_r_arms(ctx, dirs, groups, 0.2, held, held)["m1_r_gmom"]
    assert arms["m1_r_gmom"].X.tobytes() == direct.X.tobytes()
    assert arms["idle"].X.tobytes() == arms["m1_r_gmom"].X.tobytes()
    for a in sr.MULTIPLIERS:
        assert arms[a].X.shape == arms["m1_r_gmom"].X.shape
        assert not np.array_equal(arms[a].X, arms["m1_r_gmom"].X)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _arm(score, mmd_u=0.1):
    return {"score": score, "groups": {"distribution": score / 100},
            "metrics": {"mmd_u": mmd_u, "variogram": 0.01},
            "skills": {"de_score": score / 200, "de_direction": 0.4, "mmd_u": 0.3, "variogram": 0.2}}


B0, BASE = 50.0, 53.0


def _task(gains, seeds=range(10), commit="c" * 40, clean=True, task=None, fraction=0.2, fractions=None,
          base_shift=0.0, identical=False, per_seed=None):
    """A task file. gains: arm -> mean gain over m1_r_gmom; per_seed: arm -> list of ten per-seed gains (overrides
    gains). A scaled arm's raw metric differs from the base's by 1e-9 unless `identical`."""
    rows = []
    for i, s in enumerate(seeds):
        arms = {"b0": _arm(B0), "m1_r_gmom": _arm(BASE + base_shift)}
        for a in ("m1_r_gmom_x15", "m1_r_gmom_x2"):
            g = per_seed[a][i] if per_seed and a in per_seed else gains.get(a, 0.0)
            arms[a] = _arm(BASE + base_shift + g, mmd_u=0.1 if identical else 0.1 + 1e-9)
        rows.append({"seed": s, "floor": 50.0, "arms": arms})
    fr = fractions or {"m1_r_gmom": fraction, "m1_r_gmom_x15": fraction * 1.5, "m1_r_gmom_x2": fraction * 2}
    return {"label": "synthetic", "task": task, "fraction": fraction, "fractions": fr, "seeds": list(seeds),
            "rows": rows, "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}


def _write(sr, tmp_path, cuts, js=-0.13):
    """cuts: cut key -> task record. T42's committed files, synthetic, beside them."""
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k, rec in cuts.items():
        rec = {**rec, "task": k if rec.get("task") is None else rec["task"]}
        (out / "tasks" / sr.CUT_FILES[k]).write_text(json.dumps(rec))
    t42 = tmp_path / "t42"
    t42.mkdir()
    for f in sr.CUT_FILES.values():
        rows = [{"seed": s, "arms": {"b0": _arm(B0), "m1_r_gmom": _arm(BASE), "m1_r_gmom_js": _arm(BASE + js)}}
                for s in range(10)]
        (t42 / f).write_text(json.dumps({"rows": rows}))
    return out, t42


def _read(sr, tmp_path, gains=None, per_seed=None, **over):
    gains = gains or {}
    cuts = {k: _task(gains, per_seed=per_seed) for k in sr.CUT_FILES}
    for k, kw in over.items():
        if kw is None:
            del cuts[getattr(sr.gc, k)]
        else:
            cuts[getattr(sr.gc, k)] = _task(gains, **kw)
    out, t42 = _write(sr, tmp_path, cuts)
    return sr.read(out, t42)["reading"]


@pytest.mark.parametrize("gains, verdict", [
    ({"m1_r_gmom_x15": 1.0}, "a larger step scores better on these cuts: m1_r_gmom_x15"),
    ({"m1_r_gmom_x2": 1.0}, "a larger step scores better on these cuts: m1_r_gmom_x2"),
    ({"m1_r_gmom_x15": 0.5, "m1_r_gmom_x2": 1.0}, "a larger step scores better on these cuts: m1_r_gmom_x15, "
                                                    "m1_r_gmom_x2"),
    ({"m1_r_gmom_x15": -0.1, "m1_r_gmom_x2": -0.2}, "a larger step loses"),
    ({"m1_r_gmom_x15": 0.1, "m1_r_gmom_x2": -0.2}, "the step's size is not separated on these cuts"),
    ({}, "the step's size is not separated on these cuts"),
])
def test_each_verdict_is_reached_by_a_table_built_for_it(sr, tmp_path, gains, verdict):
    r = _read(sr, tmp_path, gains)
    assert r["valid"] and r["outcome"] == verdict
    for a in sr.MULTIPLIERS:
        comp = r["comparisons"][a]
        assert comp["base"] == "m1_r_gmom" and comp["change"] == a
        cut = comp["cuts"][sr.gc.TASK_CHAMBER]
        want = gains.get(a, 0.0)
        assert cut["mean"] == pytest.approx(want) and cut["ahead"] == (10 if want > 0 else 0)
        assert cut["range"] == pytest.approx([want, want])
        assert set(cut["skills"]) == {"de_score", "de_direction", "mmd_u", "variogram"}
        assert cut["skills"]["de_score"] == pytest.approx(want / 200)
        assert r["mean_scores"][sr.gc.TASK_WHOLE][a] == pytest.approx(BASE + want)
    sc = r["scalar_control_t42"]
    assert sc["base"] == "m1_r_gmom" and sc["change"] == "m1_r_gmom_js"
    assert sc["cuts"][sr.gc.TASK_WHOLE]["mean"] == pytest.approx(-0.13)


def test_losing_needs_at_most_two_seeds_ahead_on_each_cut(sr, tmp_path):
    three_ahead = [0.05, 0.05, 0.05] + [-0.3] * 7
    two_ahead = [0.05, 0.05] + [-0.3] * 8
    r = _read(sr, tmp_path, per_seed={"m1_r_gmom_x15": two_ahead, "m1_r_gmom_x2": two_ahead})
    assert r["valid"] and r["outcome"] == sr.LOSES
    assert r["comparisons"]["m1_r_gmom_x2"]["cuts"][sr.gc.TASK_CHAMBER]["range"] == pytest.approx([-0.3, 0.05])
    r = _read(sr, tmp_path / "more", {"m1_r_gmom_x15": -0.1, "m1_r_gmom_x2": -0.2},
              TASK_WHOLE={"per_seed": {"m1_r_gmom_x15": three_ahead, "m1_r_gmom_x2": two_ahead}})
    assert r["valid"] and r["outcome"] == sr.NOT_SEPARATED


def test_a_gain_of_exactly_point_three_does_not_pass(sr, tmp_path):
    r = _read(sr, tmp_path, {"m1_r_gmom_x15": 0.3, "m1_r_gmom_x2": 0.3})
    assert r["valid"] and not any(c["passes"] for c in r["comparisons"].values())
    assert r["outcome"] == sr.NOT_SEPARATED
    r = _read(sr, tmp_path / "more", {"m1_r_gmom_x2": 0.31})
    assert r["comparisons"]["m1_r_gmom_x2"]["passes"] and r["outcome"] == f"{sr.BETTER}: m1_r_gmom_x2"


def test_a_loss_beyond_point_three_on_one_cut_vetoes(sr, tmp_path):
    cuts = {sr.gc.TASK_CHAMBER: _task({"m1_r_gmom_x2": 1.0}), sr.gc.TASK_WHOLE: _task({"m1_r_gmom_x2": -0.4})}
    out, t42 = _write(sr, tmp_path, cuts)
    r = sr.read(out, t42)["reading"]
    assert r["valid"] and not r["comparisons"]["m1_r_gmom_x2"]["passes"]
    assert r["outcome"] == sr.NOT_SEPARATED


@pytest.mark.parametrize("over, reason", [
    ({"TASK_WHOLE": None}, "missing cuts"),
    ({"TASK_CHAMBER": {"seeds": range(9)}}, "seeds"),
    ({"TASK_WHOLE": {"base_shift": 1e-6}}, "do not reproduce"),
    ({"TASK_WHOLE": {"fraction": 0.4}}, "fraction 0.4"),
    ({"TASK_WHOLE": {"fractions": {"m1_r_gmom": 0.2, "m1_r_gmom_x15": 0.3, "m1_r_gmom_x2": 0.6}}}, "wrong fraction"),
    ({"TASK_CHAMBER": {"commit": "d" * 40}}, "2 commits"),
    ({"TASK_WHOLE": {"clean": False}}, "not clean"),
    ({"TASK_CHAMBER": {"task": "T1/E9.5->E10.5"}}, "the file holds task"),
    ({"TASK_CHAMBER": {"identical": True}, "TASK_WHOLE": {"identical": True}}, "exactly as m1_r_gmom"),
])
def test_the_void_cases(sr, tmp_path, over, reason):
    r = _read(sr, tmp_path, {"m1_r_gmom_x2": 1.0}, **over)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]


def test_a_base_arm_off_its_committed_score_by_one_e_minus_six_voids_and_b0_too(sr, tmp_path):
    cuts = {k: _task({"m1_r_gmom_x2": 1.0}) for k in sr.CUT_FILES}
    for row in cuts[sr.gc.TASK_CHAMBER]["rows"][:1]:
        row["arms"]["b0"]["score"] += 1e-6
    out, t42 = _write(sr, tmp_path, cuts)
    r = sr.read(out, t42)["reading"]
    assert not r["valid"] and r["validity"]["reproduction"][sr.gc.TASK_CHAMBER]["largest_difference"]["b0"] > 1e-9
    r = _read(sr, tmp_path / "ok", {"m1_r_gmom_x2": 1.0}, TASK_WHOLE={"base_shift": 1e-12})
    assert r["valid"]


def test_a_scaled_arm_identical_on_one_cut_only_stays_valid(sr, tmp_path):
    r = _read(sr, tmp_path, {"m1_r_gmom_x2": 1.0}, TASK_WHOLE={"identical": True})
    assert r["valid"]


def test_a_missing_t42_file_voids(sr, tmp_path):
    out, t42 = _write(sr, tmp_path, {k: _task({}) for k in sr.CUT_FILES})
    (t42 / sr.CUT_FILES[sr.gc.TASK_WHOLE]).unlink()
    r = sr.read(out, t42)["reading"]
    assert not r["valid"] and any("do not reproduce" in x for x in r["validity"]["reasons"])


def test_the_reading_refuses_to_overwrite(sr, tmp_path):
    out, t42 = _write(sr, tmp_path, {k: _task({}) for k in sr.CUT_FILES})
    sr.read(out, t42)
    assert (out / "t44.json").exists() and (out / "t44.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        sr.read(out, t42)


def test_the_run_needs_its_row(sr, monkeypatch):
    monkeypatch.setattr(sr.n12, "tries_rows", lambda: {"T26", "T42"})
    monkeypatch.setattr(sr.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T44 has no row"):
        sr.main(["run", "--task", sr.gc.TASK_CHAMBER, "--out", "t44"])


def test_the_run_needs_n22(sr, monkeypatch):
    monkeypatch.setattr(sr.n12, "tries_rows", lambda: {"T44"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        sr.main(["run", "--task", sr.gc.TASK_CHAMBER, "--out", "t44", "--n22", "n22x"])


def test_the_run_refuses_an_existing_task_file(sr, monkeypatch):
    monkeypatch.setattr(sr.n12, "tries_rows", lambda: {"T44"})
    monkeypatch.setattr(sr.gc, "n22_passed", lambda path: True)
    tasks = sr.RESULTS / "t44" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / sr.CUT_FILES[sr.gc.TASK_CHAMBER]).write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        sr.main(["run", "--task", sr.gc.TASK_CHAMBER, "--out", "t44"])
