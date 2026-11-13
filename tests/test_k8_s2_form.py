"""Tries K8a and K8b (experiments/harness/TRIES.md; experiments/t3_decomp/k8_s2_form.py): the second upload's own form
with the wild-type-only factor of tries K7. Synthetic data and committed records only (no competition file, no score
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
import k7_wt_revert as k7  # noqa: E402
import k8_s2_form as k8  # noqa: E402


def test_the_base_arm_is_the_frozen_second_upload():
    dec = json.loads((T3_DIR / "results" / "gata4" / "gata4_decisions.json").read_text())
    s2 = dec["submissions"]["refused"]["s1"]  # T3's second upload (docs/S2_RULES.md): the refused branch's s1
    assert s2["factor"] == {"source": "none", "exponent": 0.0}
    assert tuple(s2["composition"]["depleted_states"]) == k8.DEPLETED_STATES
    assert s2["composition"]["log2_magnitude"] == k8.LOG2_MAGNITUDE and s2["composition"]["kind"] == "central"
    spec = k8.spec()
    assert spec["factor"] == s2["factor"]
    assert {k: spec["composition"][k] for k in ("kind", "depleted_states", "log2_magnitude")} == {
        k: s2["composition"][k] for k in ("kind", "depleted_states", "log2_magnitude")}
    bs = harness._t3_head()
    bs.validate_spec(spec, bs.load_rules())  # the frozen builder accepts it as written


def test_tasks_and_arms():
    assert k8.TASKS == {n: k7.TASKS[n] for n in ("g1", "g2")}
    assert k8.ARMS == ("s2_form", "s2_form_revert", "s2_form_revert_grp") and k8.STREAM == "t3_s2_form"
    cands = k8.candidates()
    assert [c.name for c in cands] == [harness.FLOOR, *k8.ARMS]
    assert all(c.applies(t) for c in cands[1:] for t in k8.TASKS.values())
    assert not any(c.applies(k7.TASKS["m"]) for c in cands[1:])  # the scenario is defined on WT E8.75's vocabulary


def _result(key: str, base: float, pooled: float, grouped: float, clean: bool = True) -> dict:
    sc = {harness.FLOOR: 50.0, k8.BASE: base, k8.POOLED: pooled, k8.GROUPED: grouped}
    rows = [{"seed": s, "scores": dict(sc), "skills": {a: {"de_score": 0.5} for a in sc},
             "t3": {a: {"severity": {"state": "scored" if a != k8.BASE else "r2_gate"}} for a in sc}} for s in SEEDS]
    return {"task": {"key": key}, "rows": rows, "provenance": {"commit": "abc", "tree_clean": clean}}


def test_the_reading_asks_both_replicates(tmp_path):
    def write(d, **per):
        d.mkdir()
        for n, v in per.items():
            (d / f"{n}.json").write_text(json.dumps(_result(k8.TASKS[n].key, *v)))
        return k8.reading(d)

    both = write(tmp_path / "both", g1=(51.0, 55.0, 56.0), g2=(49.0, 52.0, 55.0))
    assert not both["void"] and both["carries"] == {k8.POOLED: True, k8.GROUPED: True}
    assert both["tasks"]["g2"]["paired"][f"{k8.GROUPED} - {k8.BASE}"]["mean"] == pytest.approx(6.0)
    assert "No request follows" in k8.text(both)
    one = write(tmp_path / "one", g1=(51.0, 55.0, 56.0), g2=(49.0, 49.2, 55.0))
    assert one["carries"] == {k8.POOLED: False, k8.GROUPED: True}
    assert write(tmp_path / "missing", g1=(51.0, 55.0, 56.0))["void"]
    d = tmp_path / "dirty"
    d.mkdir()
    for n in k8.TASKS:
        (d / f"{n}.json").write_text(json.dumps(_result(k8.TASKS[n].key, 51.0, 55.0, 56.0, clean=n != "g2")))
    assert k8.reading(d)["void"]
