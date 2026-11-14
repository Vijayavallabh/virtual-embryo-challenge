"""Reading N52 (experiments/harness/n52_defect_conditional_readings.py): committed scores split by row N51's flag, on
synthetic score tables and a synthetic flag table built here. Nothing is read from data/ or from a committed run."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HARNESS = Path(__file__).resolve().parents[1] / "experiments" / "harness"
RAMP = 0.1 * np.arange(10)


@pytest.fixture(scope="module")
def n52():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n52_defect_conditional_readings",
                                                  HARNESS / "n52_defect_conditional_readings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def names(n52):
    return [*n52.nt.DEEP, *n52.nt.OTHER]


def world(n52, flagged: dict | None = None, gains: dict | None = None):
    """Scores on the nine tasks: fgw at 60 + a ramp, the second upload at 61 everywhere, each arm at fgw + its
    per-seed gain (default: +2 on even seeds, -1 on odd ones). Flags: seeds 0-3 on every task unless given."""
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names(n52)}
    default = np.where(np.arange(10) % 2 == 0, 2.0, -1.0)
    scores = {}
    for t in tasks.values():
        s = {"floor": np.full(10, 50.0), n52.ENTRY: 60.0 + RAMP, n52.SECOND: np.full(10, 61.0)}
        for a in n52.ARMS:
            g = (gains or {}).get(a, default)
            if g is not None:
                s[a] = s[n52.ENTRY] + np.asarray(g, float)
        scores[t.key] = s
    flagged = flagged or {}
    n51 = {"commit": "c", "clean": True, "tasks": {
        n: {"key": t.key, "draws": [{"seed": s, "defect": s in flagged.get(n, (0, 1, 2, 3))} for s in range(10)]}
        for n, t in tasks.items()}}
    return scores, tasks, n51


def test_gains_and_counts_by_flag_by_hand(n52):
    scores, tasks, n51 = world(n52)
    res = n52.compute(scores, tasks, n51)
    assert not res["void"]
    t = res["arms"]["fgw_partial"]["tasks"]["half"]
    assert t["n_flagged"] == 4
    assert t["fgw_mean_flagged"] == pytest.approx(60.15)  # seeds 0-3
    assert t["fgw_mean_unflagged"] == pytest.approx(60 + np.mean(RAMP[4:]))
    assert t["gain_flagged"] == pytest.approx(0.5) and t["ahead_flagged"] == 2  # +2, -1, +2, -1
    assert t["gain_unflagged"] == pytest.approx(0.5) and t["ahead_unflagged"] == 3  # seeds 4-9
    # v = max(arm, 61) - max(fgw, 61); fgw < 61 on every seed
    fgw = 60 + RAMP
    arm = fgw + np.where(np.arange(10) % 2 == 0, 2.0, -1.0)
    v = np.maximum(arm, 61) - 61
    assert t["V_flag"] == pytest.approx(v[:4].mean()) and t["V_rest"] == pytest.approx(v[4:].mean())
    assert t["V_all"] == pytest.approx(v.mean())


def test_a_task_with_two_flagged_draws_is_left_out_and_named(n52):
    scores, tasks, n51 = world(n52, flagged={"x1": (0, 1), "e75": tuple(range(9))})
    r = n52.compute(scores, tasks, n51)["arms"]["fgw_damp"]
    assert "x1" not in r["flagged"]["tasks"] and r["flagged"]["left_out"] == ["x1"]
    assert r["unflagged"]["left_out"] == ["e75"]  # one unflagged draw
    assert "x1" in r["unflagged"]["tasks"]


def test_s_flag_by_hand(n52):
    nt = n52.nt
    gains = {a: np.zeros(10) for a in n52.ARMS}
    scores, tasks, n51 = world(n52, flagged={"x2": (0, 1)}, gains=gains)
    for i, n in enumerate(names(n52)):  # the arm at 70 on flagged seeds 0-3: v = 70 - 61 = 9 there, plus i / 10
        scores[tasks[n].key]["fgw_still"][:4] = 70.0 + i / 10
    r = n52.compute(scores, tasks, n51)["arms"]["fgw_still"]["flagged"]
    vflag = {n: 9.0 + i / 10 for i, n in enumerate(names(n52))}
    vflag["x2"] = 9.0 + names(n52).index("x2") / 10  # two flagged draws: left out
    deep = [vflag[n] for n in nt.DEEP if n != "x2"]
    other = [vflag[n] for n in nt.OTHER]
    d, o = np.mean(deep), np.mean(other)
    assert r["D"] == pytest.approx(d) and r["O"] == pytest.approx(o)
    for p in nt.PS:
        assert r[f"S{p}"] == pytest.approx(p * d + (1 - p) * o)
    assert r["left_out"] == ["x2"]


def test_all_draw_values_equal_n37_price(n52):
    scores, tasks, n51 = world(n52)
    res = n52.compute(scores, tasks, n51)
    ref = n52.n37.price(scores, tasks, candidates=("fgw_partial",))["first_uploads"]["fgw_partial"]
    a = res["arms"]["fgw_partial"]["all_draws"]
    for k in ("D", "O", *(f"S{p}" for p in n52.nt.PS)):
        assert a[k] == pytest.approx(ref[k], abs=1e-12)
    assert res["check_max_abs_difference"] <= 1e-9


def test_a_score_altered_under_the_check_voids_the_reading(n52, monkeypatch):
    scores, tasks, n51 = world(n52)
    real = n52.n37.price

    def altered(s, t, candidates):
        s = {k: dict(v) for k, v in s.items()}
        s[tasks["e725"].key]["fgw_damp_m"] = s[tasks["e725"].key]["fgw_damp_m"] + 0.5
        return real(s, t, candidates=candidates)

    monkeypatch.setattr(n52.n37, "price", altered)
    res = n52.compute(scores, tasks, n51)
    assert res["void"] and any("fgw_damp_m" in w for w in res["why"])
    assert "VOID" in n52.text(res)


def test_an_arm_missing_on_a_task_is_not_read(n52):
    scores, tasks, n51 = world(n52, gains={"fgw_damp_m_grp": None})
    del scores[tasks["new2"].key]["fgw_part_rand"]
    res = n52.compute(scores, tasks, n51)
    assert not res["void"]
    assert set(res["not_read"]) == {"fgw_damp_m_grp", "fgw_part_rand"}
    assert "fgw_part_rand" not in res["arms"] and "fgw_partial" in res["arms"]


def test_upload_seed_listing(n52):
    scores, tasks, n51 = world(n52, flagged={"e75": (1,), "half": (0,)})
    r = n52.compute(scores, tasks, n51)["arms"]["fgw_partial"]["upload_seed"]
    assert r["e75"]["defect"] is True and r["half"]["defect"] is False
    # seed 1: fgw 60.1, arm 59.1 (an odd seed loses 1): v = 61 - 61 = 0
    assert r["e75"]["v"] == pytest.approx(0.0)
    scores[tasks["e75"].key]["fgw_partial"][1] = 64.0
    r = n52.compute(scores, tasks, n51)["arms"]["fgw_partial"]["upload_seed"]
    assert r["e75"]["v"] == pytest.approx(3.0)


def test_a_bad_flag_table_voids(n52):
    scores, tasks, n51 = world(n52)
    n51["clean"] = False
    del n51["tasks"]["x3"]
    res = n52.compute(scores, tasks, n51)
    assert res["void"] and any("dirty" in w for w in res["why"]) and any("x3" in w for w in res["why"])


def test_text_says_it_decides_nothing(n52):
    scores, tasks, n51 = world(n52, flagged={"x1": (0, 1)})
    txt = n52.text(n52.compute(scores, tasks, n51))
    assert "decides nothing" in txt and "left out" in txt


def test_refuses_an_existing_output(n52, tmp_path, monkeypatch):
    monkeypatch.setattr(n52, "reading", lambda p: {"void": True, "why": ["synthetic"]})
    n51 = tmp_path / "n51.json"
    n51.write_text(json.dumps({}))
    out = tmp_path / "out"
    assert n52.main(["--n51", str(n51), "--out", str(out)]) == 0
    assert set(p.name for p in out.iterdir()) == {"n52.json", "n52.txt"}
    assert "provenance" in json.loads((out / "n52.json").read_text())
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        n52.main(["--n51", str(n51), "--out", str(out)])
