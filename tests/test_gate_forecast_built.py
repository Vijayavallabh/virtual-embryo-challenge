"""experiments/gate/forecast_built.py: its committed output agrees with the committed forecasts it restates and with
the files it read. No draw is repeated here; the script itself is rerun by hand (its docstring gives the command)."""

import hashlib
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "experiments" / "gate" / "results"
BUILT = RES / "forecast_built_2026-11-13.json"
PAIR = RES / "forecast_t1_pair_2026-11-13.json"
SETS = ("files of 2026-11-13 (forecast.py, committed)", "with T1's second upload on the mid base (R7)",
        "the files as built: R7 and the heart trend (H11)")


@pytest.fixture(scope="module")
def built():
    return json.loads(BUILT.read_text())


def test_inputs_are_the_committed_files(built):
    assert built["inputs_sha256"], "the output names the files it read"
    for rel, sha in built["inputs_sha256"].items():
        assert not rel.startswith(("boards/", "data/", "runs/", "docs/research/")) and "index.json" not in rel
        assert hashlib.sha256((REPO / rel).read_bytes()).hexdigest() == sha, f"{rel} changed since the run"


def test_first_two_sets_are_the_committed_pair_forecast(built):
    pair = json.loads(PAIR.read_text())["table"]["central"]
    sets = built["files"]["sets"]
    assert tuple(sets) == SETS
    for name, arm in ((SETS[0], "built"), (SETS[1], "mid_s2")):
        assert sets[name]["mean"] == pytest.approx(pair[arm]["total_mean"], abs=1e-9)
        assert sets[name]["boards_better_of_two"]["T1"] == pytest.approx(pair[arm]["t1_mean_better_of_two"], abs=1e-9)
        for line in ("203.00", "205.32", "209.00"):
            assert sets[name]["p_at_least"][line] == pytest.approx(pair[arm][f"p_{line}"], abs=1e-9)


def test_heart_moves_by_the_trends_committed_gain(built):
    f = built["files"]
    gains = [v["mean"] for v in f["heart_trend_gain"].values()]
    assert len(gains) == 2 and all(v["ahead"] == v["seeds"] == 10 for v in f["heart_trend_gain"].values())
    assert f["heart_shift"] == pytest.approx(sum(gains) / 2, abs=1e-12)
    a, b, c = (f["sets"][n]["boards_better_of_two"] for n in SETS)
    assert c["heart"] == pytest.approx(b["heart"] + f["heart_shift"], abs=1e-9)
    assert all(c[k] == pytest.approx(b[k], abs=1e-12) for k in ("T1", "embryo", "T3"))
    assert all(b[k] == pytest.approx(a[k], abs=1e-12) for k in ("heart", "embryo", "T3"))
    assert f["sets"][SETS[2]]["mean"] == pytest.approx(f["sets"][SETS[1]]["mean"] + f["heart_shift"] / 2, abs=0.01)


def test_what_it_takes_is_consistent(built):
    w = built["what_it_takes"]
    weights = {"T1": 1.0, "heart": 0.5, "embryo": 0.5, "T3": 1.0}
    assert w["sum_of_maxima"] == pytest.approx(sum(weights[b] * v["max"] for b, v in w["boards"].items()), abs=1e-9)
    for b, v in w["boards"].items():
        q = v["q05_q25_q50_q75_q95"]
        assert q == sorted(q) and q[-1] <= v["max"] and 0 <= q[0], b
    base = built["files"]["sets"][SETS[2]]
    for name, s in w["board_shift"].items():
        board, delta = name.split(" +")
        assert s["mean"] == pytest.approx(base["mean"] + weights[board] * float(delta), abs=0.02), name
        assert all(s["p_at_least"][c] >= base["p_at_least"][c] for c in base["p_at_least"]), name


def test_no_test_score_or_leaderboard_figure_is_an_input(built):
    assert "not a try" in built["label"] and "no test score" in built["label"]
    assert set(built["lines"]) == {"197.00", "200.00", "203.00", "205.32", "209.00", "215.41"}
