"""experiments/s2_hedge/payoff.py: the heart task lookup, PLAN 4.3 step 5's readout and summarise's refusals."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def payoff():
    spec = importlib.util.spec_from_file_location("s2_hedge_payoff", ROOT / "experiments" / "s2_hedge" / "payoff.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_heart_tasks(payoff):
    assert payoff.heart_task(payoff.HEART_TASK).counted
    t = payoff.heart_task("T2h/E8.25_late->E8.5")  # step 5's second transition: checked, never counted
    assert (t.kind, t.reference, t.target, t.panel_board) == (
        "extrap", "T2/E8.25_late.h5ad", "T2/E8.5.h5ad", "T2:heart:val_interp")
    assert not t.counted
    with pytest.raises(SystemExit):
        payoff.heart_task("T2h/E7.5->E8.5")


def _write(d: Path, task: str, rows: list[dict], name: str = "heart-scope.json", **extra) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    (d / name).write_text(json.dumps({"task": task, "axis": "heart-scope",
                                      "pred_values": ["hold", "a", "b", "c"], "rows": rows, **extra}))
    return d


def _row(seed: int, gains: dict[str, tuple[float, float]]) -> dict:
    own = {"hold": {"d2_shape": 0.3, "occupancy_dice": 0.5}}
    own.update({p: {"d2_shape": 0.3 + g2, "occupancy_dice": 0.5 + gd} for p, (g2, gd) in gains.items()})
    return {"seed": seed, "skills": {"own": own}}


def test_step5_bar(payoff, tmp_path):
    gains = {"a": (0.2, 0.01), "b": (0.2, -0.01), "c": (0.1, 0.02)}
    runs = [_write(tmp_path / "0-1", "T", [_row(0, gains), _row(1, gains)]),
            _write(tmp_path / "2-3", "T", [_row(2, gains), _row(3, gains)])]
    res = payoff.step5(runs, tmp_path / "out")
    assert res["seeds"] == [0, 1, 2, 3]
    assert {p: r["passes"] for p, r in res["preds"].items()} == {"a": True, "b": False, "c": False}
    assert res["preds"]["a"]["d2_gain"] == pytest.approx(0.2) and res["preds"]["a"]["d2_ahead"] == 4
    assert (tmp_path / "out" / "step5.txt").exists()


def test_mixed_tasks_refused(payoff, tmp_path):
    rows = [_row(0, {"a": (0.2, 0.0), "b": (0.0, 0.0), "c": (0.0, 0.0)})]
    runs = [_write(tmp_path / "x", "T", rows), _write(tmp_path / "y", "U", rows)]
    with pytest.raises(SystemExit):
        payoff.step5(runs)
    with pytest.raises(SystemExit):
        payoff.summarise(runs, tmp_path / "out")


def test_hedge_extra_priors(payoff, tmp_path):
    rows = [{"seed": 0, "scores": {"own": {"floor": 50.0, "s1": 60.0, "x": 55.0},
                                   "far": {"floor": 50.0, "s1": 40.0, "x": 70.0}}},
            {"seed": 1, "scores": {"own": {"floor": 50.0, "s1": 62.0, "x": 64.0},
                                   "far": {"floor": 50.0, "s1": 42.0, "x": 50.0}}}]
    d = tmp_path / "r"
    d.mkdir()
    (d / "part.json").write_text(json.dumps({"axis": "test", "pred_values": ["s1", "x"],
                                             "truth_values": ["own", "far"], "rows": rows}))
    res = payoff.hedge([d], "s1", tmp_path / "out", {"released": ["own"], "far_only": ["far"]})
    x = res["candidates"]["x"]
    assert x["released"]["mean"] == pytest.approx((60 + 64) / 2)  # the better of s1 and x, per seed, on own
    assert x["far_only"]["mean"] == pytest.approx((70 + 50) / 2)
    assert x["uniform"]["mean"] == pytest.approx(((60 + 70) / 2 + (64 + 50) / 2) / 2)
    assert res["candidates"]["s1"]["released"]["mean"] == pytest.approx(61.0)  # s2 = s1 is s1 alone
    for bad in ({"bad": ["nope"]}, {"empty": []}, {"uniform": ["own"]}):
        with pytest.raises(SystemExit):
            payoff.hedge([d], "s1", None, bad)
