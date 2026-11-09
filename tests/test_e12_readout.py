"""experiments/harness/e12_readout.py (try E12 of experiments/harness/TRIES.md), on synthetic harness results only: the
row's kill list at its boundaries, the better-of-two, the mechanism read and the outputs (stamped, one run only)."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import adoption, harness
from vec.evaluate import SUBSAMPLE_FRAC

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def e12():
    spec = importlib.util.spec_from_file_location("harness_e12_readout",
                                                  REPO / "experiments" / "harness" / "e12_readout.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


METRICS = ["de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd"]
SEEDS = list(adoption.SEEDS)


def _skills(base: float, **over) -> dict:
    return {m: over.get(m, base) for m in METRICS}


def make_run(key: str, fgw, shape0, cand, *, commit="abc123", clean=True, skill_over=None, seeds=SEEDS) -> dict:
    """One harness task JSON as experiments/harness/run.py writes it (only what the readout and check_results read)."""
    skill_over = skill_over or {}
    rows = []
    for i, s in enumerate(seeds):
        cand_skills = _skills(0.6, **{m: v[i] for m, v in skill_over.items()})
        rows.append({"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s,
                     "scores": {"floor": 50.0, "fgw": float(fgw[i]), "fgw_shape0": float(shape0[i]),
                                "fgw_shape0_aligned": float(cand[i])},
                     "skills": {"floor": _skills(0.5), "fgw": _skills(0.6), "fgw_shape0": _skills(0.6),
                                "fgw_shape0_aligned": cand_skills}})
    return {"task": {"key": key}, "board": {"key": f"board:{key}"}, "candidates": ["floor", "fgw", "fgw_shape0",
                                                                                     "fgw_shape0_aligned"],
            "protocol": "bake", "build": "main", "anchors": "computed", "frac": SUBSAMPLE_FRAC, "seeds": list(seeds),
            "rows": rows,
            "provenance": {"commit": commit, "tree_clean": clean, "config_sha256": harness._sha256(harness.CONFIG),
                           "files_sha256": {f"T2/{key.split('/')[1][:5]}.h5ad": "x" * 64}, "derived_sha256": {},
                           "utc": "2026-11-09 12:00:00"}, "run_id": "e12"}


FGW_75, FGW_725 = np.full(10, 61.0), np.full(10, 58.8)
S0_75, S0_725 = np.full(10, 63.0), np.full(10, 53.7)
OK_CAND_75 = S0_75 + 0.6  # above +0.3 and ahead in every seed
OK_CAND_725 = S0_725 + 0.2


def runs_pair(e12, cand_75=OK_CAND_75, cand_725=OK_CAND_725, **kw):
    r75 = make_run(e12.TASK_E75, FGW_75, S0_75, cand_75, **kw)
    r725 = make_run(e12.TASK_E725, FGW_725, S0_725, cand_725)
    return {e12.TASK_E75: r75, e12.TASK_E725: r725}


GOOD_RECORD = {"entry": "fgw_shape0", "board": "T2:embryo:test_interp", "dry_run": False, "n_cells": "max_cells",
               "record": {"t_pos": 0.0, "move_coupling": "aligned"}}


def test_a_candidate_that_clears_every_kill(e12):
    kills = e12.read_kills(runs_pair(e12), GOOD_RECORD)
    assert not any(v["killed"] for v in kills.values())
    # the better-of-two: E7.5's max is the candidate (63.6), E7.25's is fgw (58.8), so 61.2, above 61.02
    assert kills["a"]["candidate"] == pytest.approx((63.6 + 58.8) / 2)
    assert kills["a"]["fgw_shape0"] == pytest.approx((63.0 + 58.8) / 2)
    assert kills["b"]["mean_gain"] == pytest.approx(0.6) and kills["b"]["seeds_ahead"] == 10
    assert e12.verdict_line(kills, {"ok": True, "max_abs_diff": {}}).startswith("CLEARS THE KILL LIST")


def test_kill_a_is_the_better_of_two_bar_61_02(e12):
    # E7.5 gain +0.1: the better-of-two is (63.1 + 58.8) / 2 = 60.95, below the bar
    kills = e12.read_kills(runs_pair(e12, cand_75=S0_75 + 0.1), GOOD_RECORD)
    assert kills["a"]["killed"] and kills["a"]["candidate"] == pytest.approx(60.95)
    assert e12.BAR_BETTER_OF_TWO == 61.02
    # the bar itself: a shade above it passes, a shade below it does not
    above = e12.read_kills(runs_pair(e12, cand_75=np.full(10, 2 * 61.03 - 58.8)), GOOD_RECORD)
    below = e12.read_kills(runs_pair(e12, cand_75=np.full(10, 2 * 61.01 - 58.8)), GOOD_RECORD)
    assert not above["a"]["killed"] and below["a"]["killed"]


def test_kill_a_counts_the_candidate_where_it_beats_fgw_on_e725(e12):
    # the better-of-two takes the larger of fgw's and the candidate's score seed by seed
    cand_725 = FGW_725.copy()
    cand_725[:5] = FGW_725[:5] + 2.0
    runs = runs_pair(e12, cand_725=cand_725)
    kills = e12.read_kills(runs, GOOD_RECORD)
    assert kills["a"]["candidate_per_task"][e12.TASK_E725] == pytest.approx(58.8 + 1.0)
    assert kills["a"]["candidate"] == pytest.approx((63.6 + 59.8) / 2)


@pytest.mark.parametrize("gains, killed", [
    ([0.6] * 10, False),
    ([0.3] * 10, True),  # the adoption rule's form needs a mean above 0.3 strictly
    ([1.2] * 7 + [-0.2] * 3, True),  # mean above 0.3, ahead in 7 of 10 only
    ([0.4] * 9 + [0.0], False),  # mean 0.36, ahead in 9
    ([0.4] * 8 + [-0.05] * 2, False),  # mean 0.31, ahead in 8
])
def test_kill_b_is_the_adoption_rules_form_on_e75(e12, gains, killed):
    cand = S0_75 + np.array(gains)
    kills = e12.read_kills(runs_pair(e12, cand_75=cand), GOOD_RECORD)
    assert kills["b"]["killed"] == killed
    assert kills["b"]["mean_gain"] == pytest.approx(np.mean(gains))


def _skill_diff(metric, values):
    """Candidate skills 0.6 plus a per-seed difference on one metric (the others identical to the entry's)."""
    return {metric: [0.6 + v for v in values]}


@pytest.mark.parametrize("metric", ["de_score", "de_direction", "mmd_u", "variogram"])
def test_kill_c_reads_every_expression_metric_on_both_tasks(e12, metric):
    tiny = [1e-7] + [0.0] * 9
    runs = runs_pair(e12)
    r75 = make_run(e12.TASK_E75, FGW_75, S0_75, OK_CAND_75, skill_over=_skill_diff(metric, tiny))
    assert not e12.read_kills({e12.TASK_E75: r75, e12.TASK_E725: runs[e12.TASK_E725]}, GOOD_RECORD)["c"]["killed"]
    big = [0.0] * 9 + [2e-6]
    r725 = make_run(e12.TASK_E725, FGW_725, S0_725, OK_CAND_725, skill_over=_skill_diff(metric, big))
    kills = e12.read_kills({e12.TASK_E75: runs[e12.TASK_E75], e12.TASK_E725: r725}, GOOD_RECORD)
    assert kills["c"]["killed"] and kills["c"]["max_abs_skill_diff"][metric] == pytest.approx(2e-6)


@pytest.mark.parametrize("change, killed", [(0.02, False), (-0.009, False), (-0.02, True)])
def test_kill_d_is_the_neighbourhood_skill_on_e75_against_the_entry(e12, change, killed):
    over = _skill_diff("neighborhood_mmd", [change] * 10)
    runs = runs_pair(e12)
    r75 = make_run(e12.TASK_E75, FGW_75, S0_75, OK_CAND_75, skill_over=over)
    kills = e12.read_kills({e12.TASK_E75: r75, e12.TASK_E725: runs[e12.TASK_E725]}, GOOD_RECORD)
    assert kills["d"]["killed"] == killed and kills["d"]["mean_change"] == pytest.approx(change)
    # E7.25 does not enter (d)
    over_725 = _skill_diff("neighborhood_mmd", [-0.5] * 10)
    r725 = make_run(e12.TASK_E725, FGW_725, S0_725, OK_CAND_725, skill_over=over_725)
    assert not e12.read_kills({e12.TASK_E75: runs[e12.TASK_E75], e12.TASK_E725: r725}, GOOD_RECORD)["d"]["killed"]


def test_a_non_finite_skill_is_a_kill_not_a_pass(e12):
    # kill (c) on any expression metric of either task, kill (d) on E7.5's neighbourhood skill
    runs = runs_pair(e12)
    nan_75 = make_run(e12.TASK_E75, FGW_75, S0_75, OK_CAND_75, skill_over={"mmd_u": [float("nan")] + [0.6] * 9})
    assert e12.read_kills({e12.TASK_E75: nan_75, e12.TASK_E725: runs[e12.TASK_E725]}, GOOD_RECORD)["c"]["killed"]
    nan_nb = make_run(e12.TASK_E75, FGW_75, S0_75, OK_CAND_75,
                      skill_over={"neighborhood_mmd": [float("nan")] + [0.6] * 9})
    assert e12.read_kills({e12.TASK_E75: nan_nb, e12.TASK_E725: runs[e12.TASK_E725]}, GOOD_RECORD)["d"]["killed"]


@pytest.mark.parametrize("record, killed", [
    (GOOD_RECORD, False),
    (None, True),
    ({**GOOD_RECORD, "entry": "fgw_aligned"}, True),
    ({**GOOD_RECORD, "board": "T2:embryo:val_interp"}, True),
    ({**GOOD_RECORD, "dry_run": True}, True),
    ({**GOOD_RECORD, "record": {"t_pos": 0.5, "move_coupling": "aligned"}}, True),
    ({**GOOD_RECORD, "record": {"t_pos": 0.0, "move_coupling": "fgw"}}, True),
])
def test_kill_e_is_the_test_builds_guard(e12, record, killed):
    assert e12.read_kills(runs_pair(e12), record)["e"]["killed"] == killed


def test_mechanism_reads_e725s_mirrored_seeds_separately(e12):
    gains = np.array([0.1, 0.2, 0.3, -1.0, -2.0, -3.0, 0.4, 0.5, -4.0, -5.0])  # the mirrored seeds 3, 4, 5, 8, 9 lose
    r = make_run(e12.TASK_E725, FGW_725, S0_725, S0_725 + gains)
    m = e12.mechanism(r)
    assert m["mirrored_seeds"]["seeds"] == [3, 4, 5, 8, 9] and m["other_seeds"]["seeds"] == [0, 1, 2, 6, 7]
    assert m["mirrored_seeds"]["mean_gain_over_fgw_shape0"] == pytest.approx(np.mean([-1, -2, -3, -4, -5]))
    assert m["other_seeds"]["seeds_ahead_of_fgw_shape0"] == 5 and m["mirrored_seeds"]["seeds_ahead_of_fgw_shape0"] == 0
    assert set(m["mirrored_seeds"]["mean_skill_change"]) == {"d2_shape", "occupancy_dice", "neighborhood_mmd"}


def test_reproduction_is_per_seed_against_the_committed_run(e12):
    runs = runs_pair(e12)
    ref = {k: json.loads(json.dumps(r)) for k, r in runs.items()}
    assert e12.reproduction(runs, ref)["ok"]
    ref[e12.TASK_E75]["rows"][4]["scores"]["fgw_shape0"] += 0.5
    out = e12.reproduction(runs, ref)
    assert not out["ok"] and max(out["max_abs_diff"].values()) == pytest.approx(0.5)
    assert not e12.reproduction(runs, None)["ok"]
    line = e12.verdict_line(e12.read_kills(runs, GOOD_RECORD), out)
    assert line.startswith("HOLD")
    # a kill does not outrank a failed reproduction: the run is held, not read
    assert e12.verdict_line(e12.read_kills(runs_pair(e12, cand_75=S0_75), GOOD_RECORD), out).startswith("HOLD")


def _write_run(folder: Path, e12, **kw) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    runs = runs_pair(e12, **kw)
    for key, r in runs.items():
        (folder / e12.FILES[key]).write_text(json.dumps(r))


def _no_counts(text: str, values: set[int]) -> None:
    text = re.sub(r"\b[0-9a-f]{7,64}\b", "<hex>", text)  # commits and sha256s hold digit runs by chance
    text = re.sub(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z?", "<utc>", text)
    found = {int(d) for d in re.findall(r"\d+", text)}
    assert not found & values, sorted(found & values)


def test_main_writes_stamped_outputs_refuses_a_rerun_and_states_the_try_count(e12, tmp_path, own_repo, capsys):
    runs_dir = tmp_path / "runs" / "emb_shape_aligned_v1"
    _write_run(runs_dir, e12)
    ref_dir = tmp_path / "ref"
    ref_dir.mkdir()
    for name in e12.FILES.values():
        (ref_dir / name).write_bytes((runs_dir / name).read_bytes())
    build_json = tmp_path / "builds" / "T2embryo_s2_aligned_e12.build.json"
    build_json.parent.mkdir()
    build_json.write_text(json.dumps(GOOD_RECORD))
    out = tmp_path / "out" / "e12_readout"
    argv = ["--runs", str(runs_dir), "--test-build-json", str(build_json), "--embryo-tries-before", "9",
            "--reference", str(ref_dir), "--out", str(out)]
    e12.main(argv)
    printed = capsys.readouterr().out
    rep = json.loads((out / "e12_readout.json").read_text())
    txt = (out / "e12_readout.txt").read_text()
    assert rep["verdict"].startswith("CLEARS THE KILL LIST") and printed == txt
    assert rep["reproduction"]["ok"] and rep["tries_before"] == 9 and rep["tries_after"] == 10
    assert "this run is one counted embryo try (9 before it, 10 after)" in txt
    assert set(rep["kills"]) == {"a", "b", "c", "d", "e"} and not any(v["killed"] for v in rep["kills"].values())
    assert rep["provenance"]["commit"] and rep["provenance"]["tree_clean"] is True and rep["provenance"]["utc"]
    hashes = rep["inputs"]
    assert hashes[str(build_json)] == harness._sha256(build_json) and len(hashes) == 5
    for name in ("e12_readout.json", "e12_readout.txt"):
        stamp = json.loads((out / f"{name}.stamp.json").read_text())
        assert stamp["sha256"] == harness._sha256(out / name) and stamp["clean"] is True
    # no cell count of any kind is echoed: the board's fields and the build record's limits never reach the outputs
    _no_counts(txt, {4321, 2468, 1357, 777, 1660, 830, 5000})
    # the mechanism read and the kill list are in the report; it refuses a second run
    assert "mirrored seeds 3, 4, 5, 8, 9" in txt and all(f"({k})" in txt for k in "abcde")
    with pytest.raises(SystemExit, match="a rerun needs a new TRIES row"):
        e12.main(argv)


def test_main_stops_on_a_kill_and_holds_when_the_reference_does_not_reproduce(e12, tmp_path, own_repo, capsys):
    runs_dir = tmp_path / "runs" / "e12"
    _write_run(runs_dir, e12, cand_75=S0_75 + 0.1)  # kills (a) and (b)
    ref_dir = tmp_path / "ref"
    ref_dir.mkdir()
    for name in e12.FILES.values():
        (ref_dir / name).write_bytes((runs_dir / name).read_bytes())
    build_json = tmp_path / "b.build.json"
    build_json.write_text(json.dumps(GOOD_RECORD))
    out1 = tmp_path / "o1"
    e12.main(["--runs", str(runs_dir), "--test-build-json", str(build_json), "--embryo-tries-before", "9",
              "--reference", str(ref_dir), "--out", str(out1)])
    rep = json.loads((out1 / "e12_readout.json").read_text())
    assert rep["verdict"].startswith("STOP, no change proposed") and "a" in rep["verdict"] and "b" in rep["verdict"]
    # a reference that the run does not reproduce holds the reading (fgw_shape0 per seed differs)
    r = json.loads((ref_dir / e12.FILES[e12.TASK_E75]).read_text())
    r["rows"][0]["scores"]["fgw_shape0"] -= 1.0
    (ref_dir / e12.FILES[e12.TASK_E75]).write_text(json.dumps(r))
    out2 = tmp_path / "o2"
    e12.main(["--runs", str(runs_dir), "--test-build-json", str(build_json), "--embryo-tries-before", "9",
              "--reference", str(ref_dir), "--out", str(out2)])
    assert json.loads((out2 / "e12_readout.json").read_text())["verdict"].startswith("HOLD")
    capsys.readouterr()


def test_main_refuses_a_run_that_mixes_commits_or_lacks_a_candidate(e12, tmp_path, own_repo):
    runs_dir = tmp_path / "runs" / "e12"
    _write_run(runs_dir, e12)
    build_json = tmp_path / "b.build.json"
    build_json.write_text(json.dumps(GOOD_RECORD))
    base = ["--runs", str(runs_dir), "--test-build-json", str(build_json), "--embryo-tries-before", "9"]
    r = json.loads((runs_dir / e12.FILES[e12.TASK_E725]).read_text())
    r["provenance"]["commit"] = "other-commit"
    (runs_dir / e12.FILES[e12.TASK_E725]).write_text(json.dumps(r))
    with pytest.raises(SystemExit, match="one clean commit"):
        e12.main([*base, "--out", str(tmp_path / "o1")])
    _write_run(runs_dir, e12)
    r = json.loads((runs_dir / e12.FILES[e12.TASK_E725]).read_text())
    r["candidates"] = ["floor", "fgw", "fgw_shape0"]
    (runs_dir / e12.FILES[e12.TASK_E725]).write_text(json.dumps(r))
    with pytest.raises(SystemExit, match="lacks"):
        e12.main([*base, "--out", str(tmp_path / "o2")])
    assert not (tmp_path / "o1").exists() and not (tmp_path / "o2").exists()


def test_main_refuses_a_run_missing_a_pseudo_task_file(e12, tmp_path, own_repo):
    runs_dir = tmp_path / "runs" / "e12"
    _write_run(runs_dir, e12)
    (runs_dir / e12.FILES[e12.TASK_E725]).unlink()
    build_json = tmp_path / "b.build.json"
    build_json.write_text(json.dumps(GOOD_RECORD))
    with pytest.raises(SystemExit, match="failed.txt"):
        e12.main(["--runs", str(runs_dir), "--test-build-json", str(build_json), "--embryo-tries-before", "9",
                  "--out", str(tmp_path / "o")])
    assert not (tmp_path / "o").exists()
