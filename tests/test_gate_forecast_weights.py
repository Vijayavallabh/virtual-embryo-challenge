"""experiments/gate/forecast_weights.py: its committed output agrees with the committed forecast it restates and is
ordered as a sensitivity table must be. No draw is repeated here; the script is rerun by hand (its docstring gives the
command)."""

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "experiments" / "gate" / "results"
WEIGHTS = RES / "forecast_weights_2026-11-14.json"
BUILT = RES / "forecast_built_2026-11-13.json"
BUILT_SET = "the files as built: R7 and the heart trend (H11)"
LEAST = "every weight at the end that suits the files least"
MOST = "every weight at the end that suits the files most"


@pytest.fixture(scope="module")
def res():
    return json.loads(WEIGHTS.read_text())


def test_the_central_case_is_the_committed_forecast_for_the_files_as_built(res):
    built = json.loads(BUILT.read_text())["files"]["sets"][BUILT_SET]
    c = res["cases"]["as committed"]
    assert c["weights_moved"] == {} and c["t1_depth_weights"] is None
    assert c["mean"] == pytest.approx(built["mean"], abs=0.02)
    for b, v in built["boards_better_of_two"].items():
        assert c["boards_better_of_two"][b] == pytest.approx(v, abs=0.02), b
    for line, p in c["p_at_least"].items():
        assert built["p_at_least"][line] == pytest.approx(p, abs=0.002), line
    assert res["draws"] == 1_000_000


def test_it_calls_itself_a_sensitivity_and_not_what_a_reply_implies(res):
    assert "a sensitivity, not what any organiser reply implies" in res["label"]


def test_each_weight_moves_only_its_own_board(res):
    base = res["cases"]["as committed"]["boards_better_of_two"]
    for name, c in res["cases"].items():
        if name in ("as committed", LEAST, MOST):
            continue
        board = {"embryo": "embryo", "T1": "T1", "T3": "T3"}[name.split(",")[0]]
        for b, v in c["boards_better_of_two"].items():
            if b != board:
                assert v == pytest.approx(base[b], abs=1e-9), (name, b)


def test_the_total_is_monotone_in_each_weight(res):
    m = {n: c["mean"] for n, c in res["cases"].items()}
    emb = [m[n] for n in ("embryo, emb_p_deep 0 (never its deeper bracket's profile)",
                          "embryo, emb_p_deep at its low 0.30", "as committed", "embryo, emb_p_deep at its high 0.60",
                          "embryo, emb_p_deep 1 (always its deeper bracket's profile)")]
    t3 = [m[n] for n in ("T3, stand-in A 0", "T3, stand-in A 0.15", "T3, stand-in A at its low 0.25", "as committed",
                         "T3, stand-in A at its high 0.75", "T3, stand-in A 1")]
    assert emb == sorted(emb) and t3 == sorted(t3)
    assert m["T1, shallower than library 3 for certain"] < m["as committed"] < m["T1, depth 1.37 times library 3's"]
    assert m["T3, stand-in A 0 and the factor's carry at its low 0.35"] < m["T3, stand-in A 0"]


def test_the_two_ends_bracket_every_case_and_neither_reaches_the_leading_total(res):
    m = {n: c["mean"] for n, c in res["cases"].items()}
    assert m[LEAST] == min(m.values()) and m[MOST] == max(m.values())
    most = res["cases"][MOST]
    assert most["weights_moved"] == {"emb_p_deep": 1.0, "t3_standin_A": 1.0, "t3_q_carry": 1.0}
    assert most["p95"] < 215.41 and most["p_at_least"]["215.41"] < 0.01
    assert m[MOST] < 209.0 < 215.41
    for c in res["cases"].values():
        assert set(c["p_at_least"]) == {"203.00", "205.32", "209.00", "215.41"}
        assert c["p05"] <= c["median"] <= c["p95"]
        p = [c["p_at_least"][k] for k in ("203.00", "205.32", "209.00", "215.41")]
        assert p == sorted(p, reverse=True)
