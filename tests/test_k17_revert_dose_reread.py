"""Try K17 (experiments/harness/TRIES.md; experiments/t3_decomp/k17_revert_dose_reread.py): K16's two runs read under
a guard that tolerates what the machines do not reproduce. Hand-made score tables only (no competition file, no score
of a fold)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from vec import harness
from vec.adoption import SEEDS

REPO = Path(__file__).resolve().parents[1]
T3_DIR = REPO / "experiments" / "t3_decomp"
if str(T3_DIR) not in sys.path:
    sys.path.insert(0, str(T3_DIR))
import k16_revert_dose as k16  # noqa: E402
import k17_revert_dose_reread as k17  # noqa: E402

FACTOR = {"groups": ["CM", "NCC"], "sd_log_phi": 0.25, "genes_at_clip": 3}


def _result(key: str, scores: dict[str, float], gated: frozenset = frozenset(), n_gated: int = 0,
            base_shift: dict | None = None, commit: str = "abc", factor: dict | None = None) -> dict:
    """Ten seeds; `scores` per arm (a constant offset per seed keeps every difference the same in every seed);
    `base_shift` {seed: by} moves every arm of that seed, so paired differences stay and the base leaves the
    committed one."""
    arms = {harness.FLOOR: 50.0, **scores}
    shift = base_shift or {}
    rows = []
    for s in SEEDS:
        rows.append({"seed": s, "scores": {a: v + 0.1 * s + shift.get(s, 0.0) for a, v in arms.items()},
                     "skills": {a: {"de_score": 0.5, "severity": 0.6} for a in arms},
                     "t3": {a: {"severity": {"state": "r2_gate" if a in gated and s < n_gated else "scored",
                                             "beta": 0.1 * (1 + s), "r2": 0.05}} for a in arms}})
    return {"task": {"key": key}, "rows": rows, "provenance": {"commit": commit, "tree_clean": True},
            "factor_record": FACTOR if factor is None else factor}


def _table(e075, e100, e150, p050=50.0, p100=50.0, base=50.0) -> dict[str, float]:
    return {"grp_e050": base, "grp_e075": e075, "grp_e100": e100, "grp_e150": e150, "perm_e050": p050,
            "perm_e100": p100}


def _refs(root):
    """The committed files the reading compares with: the base's scores (K8b on g1 and g2, K7b on m) and try K7's
    factor record for every task."""
    for n in k16.TASKS:
        files = {}
        folder, a = k16.REPRO[n]
        files.setdefault(folder, {})["rows"] = [{"seed": s, "scores": {a: 50.0 + 0.1 * s}} for s in SEEDS]
        files.setdefault(k17.FACTOR_REF, {})["factor_record"] = FACTOR
        for folder, content in files.items():
            p = root / folder / "tasks"
            p.mkdir(parents=True, exist_ok=True)
            (p / f"{n}.json").write_text(json.dumps({"rows": [], **content}))
    return root


def _runs(tmp_path, name, g, m=None, **kw):
    d = tmp_path / name
    d.mkdir()
    for n in ("g1", "g2"):
        (d / f"{n}.json").write_text(json.dumps(_result(k16.TASKS[n].key, g, **kw)))
    (d / "m.json").write_text(json.dumps(_result(k16.TASKS["m"].key, m or _table(50, 50, 50), **kw)))
    return d


RISE = _table(51.0, 52.0, 53.0)


def test_the_guard_constants_and_what_it_passes(tmp_path):
    assert (k17.GUARD_MEAN, k17.GUARD_SEED, k17.FACTOR_REF) == (0.1, 0.5, "k7_wt_revert")
    refs = _refs(tmp_path / "refs")
    exact = k17.read_run(_runs(tmp_path, "exact", RISE), refs)
    assert not exact["void"] and exact["base_distance"]["m"] == {"mean_diff": 0.0, "max_abs_diff": 0.0, "seeds_off": []}
    assert set(exact["task_sha256"]) == {"g1", "g2", "m"} and len(exact["task_sha256"]["m"]) == 64
    # what K16 voided on: two seeds off by a few thousandths
    small = k17.read_run(_runs(tmp_path, "small", RISE, base_shift={1: 0.004, 9: -0.003}), refs)
    assert not small["void"] and small["base_distance"]["g1"]["seeds_off"] == [1, 9]
    assert small["named_exponent"] == 1.5 and small["request"]
    assert k16.reading(tmp_path / "small", refs)["void"]  # K16's own guard still voids the same files
    # mixed signs: the mean of the per-seed differences is signed, and each seed is inside its limit
    mix = k17.read_run(_runs(tmp_path, "mix", RISE, base_shift={1: 0.45, 2: -0.45, 3: 0.45, 4: -0.45}), refs)
    assert not mix["void"] and mix["base_distance"]["g1"]["mean_diff"] == pytest.approx(0.0)
    assert mix["base_distance"]["g1"]["max_abs_diff"] == pytest.approx(0.45)


def test_the_guard_voids_a_base_that_drifted(tmp_path):
    refs = _refs(tmp_path / "refs")
    one_seed = k17.read_run(_runs(tmp_path, "one", RISE, base_shift={3: 0.6}), refs)
    assert one_seed["void"] and any("off the committed base" in w for w in one_seed["why"])
    up = k17.read_run(_runs(tmp_path, "up", RISE, base_shift=dict.fromkeys(SEEDS, 0.2)), refs)
    assert up["void"] and up["base_distance"]["g1"]["mean_diff"] == pytest.approx(0.2)
    down = k17.read_run(_runs(tmp_path, "down", RISE, base_shift=dict.fromkeys(SEEDS, -0.2)), refs)
    assert down["void"]
    # a committed base that is missing voids the run and reports no distance
    refs2 = _refs(tmp_path / "refs2")
    committed = json.loads((refs2 / "k7_wt_revert" / "tasks" / "m.json").read_text())
    committed["rows"] = committed["rows"][:-1]  # one seed short
    (refs2 / "k7_wt_revert" / "tasks" / "m.json").write_text(json.dumps(committed))
    short = k17.read_run(_runs(tmp_path, "short", RISE), refs2)
    assert short["void"] and any("committed base is missing or lacks a seed" in w for w in short["why"])
    assert short["base_distance"]["m"] is None
    # K16's reproduction check is restored after each read
    assert k16.reproduction.__name__ == "reproduction"


def test_the_factor_record_must_equal_k7s(tmp_path):
    refs = _refs(tmp_path / "refs")
    other = k17.read_run(_runs(tmp_path, "other", RISE, factor={**FACTOR, "sd_log_phi": 0.26}), refs)
    assert other["void"] and sum("factor record differs" in w for w in other["why"]) == 3
    refs2 = _refs(tmp_path / "refs2")
    (refs2 / "k7_wt_revert" / "tasks" / "g1.json").unlink()
    gone = k17.read_run(_runs(tmp_path, "gone", RISE), refs2)
    assert gone["void"] and any("committed factor record is missing" in w for w in gone["why"])


def test_the_outcome_needs_the_same_exponent_on_both_runs(tmp_path):
    refs = _refs(tmp_path / "refs")
    flat = _table(50.1, 50.1, 50.1)
    both = k17.reading(_runs(tmp_path, "l1", RISE), _runs(tmp_path, "d1", RISE, base_shift={2: 0.3}), refs)
    assert both["request"] and both["exponent"] == 1.5 and "priced question" in both["outcome"]
    assert "cold check" in both["outcome"]
    # named on the local run, nothing named on the DGX run
    split = k17.reading(_runs(tmp_path, "l2", RISE), _runs(tmp_path, "d2", flat), refs)
    assert not split["request"] and split["why"] == "the DGX run names no exponent"
    assert "stays at exponent 0.5" in split["outcome"]
    # the two runs name different exponents (1.5 here, 1.0 there): no request
    lower = k17.reading(_runs(tmp_path, "l3", RISE), _runs(tmp_path, "d3", _table(51.0, 52.0, 52.1)), refs)
    assert not lower["request"] and "different exponents (local 1.5, DGX 1)" in lower["why"]
    # nothing named locally
    none = k17.reading(_runs(tmp_path, "l4", flat), _runs(tmp_path, "d4", RISE), refs)
    assert not none["request"] and none["why"] == "the local run names no exponent"
    # not factor-specific on one run
    perm = _table(51.0, 52.0, 53.0, p100=51.8)
    local_perm = k17.reading(_runs(tmp_path, "l5", perm), _runs(tmp_path, "d5", RISE), refs)
    assert not local_perm["request"] and local_perm["why"] == "the local run fails the factor-specific check"
    dgx_perm = k17.reading(_runs(tmp_path, "l6", RISE), _runs(tmp_path, "d6", perm), refs)
    assert not dgx_perm["request"] and dgx_perm["why"] == "the DGX run fails the factor-specific check"
    # two commits
    two = k17.reading(_runs(tmp_path, "l7", RISE), _runs(tmp_path, "d7", RISE, commit="def"), refs)
    assert not two["request"] and "different commits" in two["why"]


def test_a_void_run_means_no_request_and_the_reading_still_writes(tmp_path):
    refs = _refs(tmp_path / "refs")
    local_void = k17.reading(_runs(tmp_path, "l1", RISE, base_shift={3: 0.6}), _runs(tmp_path, "d1", RISE), refs)
    assert not local_void["request"] and local_void["why"].startswith("the local run is void")
    off = _runs(tmp_path, "d2", RISE, base_shift=dict.fromkeys(SEEDS, 0.3))
    drifted = k17.reading(_runs(tmp_path, "l2", RISE), off, refs)
    assert not drifted["request"] and drifted["why"].startswith("the DGX run is void")
    # task files K16 itself voids: two missing, one without the base arm, one with no rows
    d = tmp_path / "d3"
    d.mkdir()
    broken = _result(k16.TASKS["g1"].key, RISE)
    for q in broken["rows"]:
        del q["scores"]["grp_e050"]
    (d / "g1.json").write_text(json.dumps(broken))
    empty = _result(k16.TASKS["m"].key, RISE)
    empty["rows"] = []
    (d / "m.json").write_text(json.dumps(empty))
    res = k17.reading(_runs(tmp_path, "l3", RISE), d, refs)
    assert not res["request"] and res["why"].startswith("the DGX run is void")
    assert res["dgx"]["base_distance"]["g1"] is None and res["dgx"]["base_distance"]["m"] is None
    t = k17.text(res)
    assert "no request" in t and "LOCAL RUN" in t and "DGX RUN" in t and "VOID" in t and "g1 n/a" in t
