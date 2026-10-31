"""vec.adoption.keep_change: the adoption rule approved on 2026-11-01."""

import math

from vec import adoption, t2_build
from vec.adoption import SEEDS, keep_change


def _task(gains, base=50.0):
    return {s: (base, base + g) for s, g in zip(SEEDS, gains, strict=True)}


CLEAR = [1.0, 0.9, 1.2, 0.8, 1.1, 0.7, 1.3, 0.9, 1.0, 0.6]  # mean 0.95, ahead in 10
NULL = [0.5, -0.5] * 5  # mean 0, ahead in 5
MARGINAL = [0.4] * 8 + [-0.1, -0.1]  # mean 0.30 exactly, not above it


def test_same_numbers_as_approved_rule_7():
    assert adoption.MIN_GAIN == t2_build.COUNT_RULE_MIN_GAIN
    assert adoption.MIN_AHEAD == t2_build.COUNT_RULE_MIN_AHEAD and len(SEEDS) == 10


def test_a_change_with_no_effect_is_not_kept():
    r = keep_change({"a": _task(NULL), "b": _task(NULL), "c": _task(NULL)})
    assert not r["keep"] and r["passed"] == 0 and r["need"] == 2


def test_two_clear_passes_of_three_keep_the_change():
    r = keep_change({"a": _task(CLEAR), "b": _task(CLEAR), "c": _task(NULL)})
    assert r["keep"] and r["passed"] == 2


def test_a_loss_on_any_pseudo_task_vetoes():
    loss = [-0.5] * 10
    r = keep_change({"a": _task(CLEAR), "b": _task(CLEAR), "c": _task(loss)})
    assert not r["keep"] and "loses more than 0.3" in r["reason"]
    small_loss = [-0.2] * 10  # a mean loss of 0.2 is not a veto
    assert keep_change({"a": _task(CLEAR), "b": _task(CLEAR), "c": _task(small_loss)})["keep"]


def test_margin_and_seed_count_both_bind():
    assert not keep_change({"a": _task(MARGINAL)})["keep"]  # mean not above 0.3
    seven_ahead = [2.0] * 7 + [-0.5] * 3  # mean 1.25 but only 7 of 10 ahead
    assert not keep_change({"a": _task(seven_ahead)})["keep"]


def test_fewer_than_three_pseudo_tasks_need_every_one():
    assert keep_change({"a": _task(CLEAR), "b": _task(CLEAR)})["keep"]
    r = keep_change({"a": _task(CLEAR), "b": _task(NULL)})
    assert not r["keep"] and r["need"] == 2
    assert keep_change({"a": _task(CLEAR)})["keep"]


def test_incomplete_or_non_finite_pseudo_tasks_keep_nothing():
    short = {s: v for s, v in _task(CLEAR).items() if s != 9}
    r = keep_change({"a": _task(CLEAR), "b": _task(CLEAR), "c": short})
    assert not r["keep"] and "incomplete" in r["reason"] and r["tasks"]["c"]["missing_seeds"] == [9]
    nan = _task(CLEAR)
    nan[3] = (50.0, math.nan)
    assert not keep_change({"a": _task(CLEAR), "b": nan})["keep"]
    assert not keep_change({})["keep"]
