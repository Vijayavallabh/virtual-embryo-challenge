"""experiments/intake/d1_check.py (decision point D1, assumption A3): the fixed rule, on synthetic boards."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
from conftest import GENES

from vec.boards import Board

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def d1():
    spec = importlib.util.spec_from_file_location("intake_d1_check", REPO / "experiments" / "intake" / "d1_check.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _board(key, task, setting, ref_cells, truth_cells=1000, anchors=None, genes=GENES):
    return Board(key=key, task=task, setting=setting, genes=list(genes), needs_coords=task != "T1", min_cells=500,
                 max_cells=20000, truth_cells=truth_cells, ref_cells=ref_cells, anchors=anchors or {})


def _all(d1, **over):
    b = {
        "T1:test": _board("T1:test", "T1", None, d1.ten_percent(24979)),
        "T2:heart:test_extrap": _board("T2:heart:test_extrap", "T2", "heart", d1.ten_percent(167864)),
        "T2:embryo:test_interp": _board("T2:embryo:test_interp", "T2", "embryo", d1.ten_percent(11665)),
        "T3:bcatenin": _board("T3:bcatenin", "T3", None, d1.ten_percent(24826)),
        "T1:val": _board("T1:val", "T1", None, 1706, anchors={"de_score": {"floor": 0.0, "ceiling": 0.8}}),
    }
    b.update(over)
    return b


def test_expected_reference_stages_pass(d1):
    res = d1.check(_all(d1))
    assert {k: v["A3"] for k, v in res.items()} == dict.fromkeys(d1.TEST_BOARDS, "pass")
    assert res["T2:embryo:test_interp"]["ref_stage_hits"] == ["E7.5"]
    assert res["T1:test"]["genes_same_as_validation_twin"] is True


def test_another_stage_fails_and_unknown_counts_stay_unresolved(d1):
    res = d1.check(_all(d1, **{
        "T2:embryo:test_interp": _board("T2:embryo:test_interp", "T2", "embryo", d1.ten_percent(31671)),  # E8.0
        "T1:test": _board("T1:test", "T1", None, 1234),
        "T3:bcatenin": _board("T3:bcatenin", "T3", None, d1.ten_percent(24826), truth_cells=None),
    }))
    assert res["T2:embryo:test_interp"]["A3"] == "fail" and res["T2:embryo:test_interp"]["ref_stage_hits"] == ["E8.0"]
    assert res["T1:test"]["A3"] == "unresolved"
    assert res["T3:bcatenin"]["A3"] == "fail" and "truth_cells" in res["T3:bcatenin"]["why"]


def test_missing_board_and_anchor_refusal(d1):
    b = _all(d1)
    del b["T2:heart:test_extrap"]
    assert d1.check(b)["T2:heart:test_extrap"] == {"present": False, "A3": "unresolved",
                                                   "why": "board not in boards/index.json"}
    leaky = _all(d1, **{"T3:bcatenin": _board("T3:bcatenin", "T3", None, 2483,
                                              anchors={"mmd_u": {"floor": 0.1, "ceiling": 0.01}})})
    with pytest.raises(SystemExit, match="still carry anchors"):
        d1.check(leaky)
