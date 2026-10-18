import numpy as np
import pytest
from conftest import simulate_stage

from vec.boards import GROUPS, RANKED_METRICS
from vec.evaluate import SUBSAMPLE_FRAC, aggregate, score_board, skills_from_anchors, subsample


def test_skill_map_anchors():
    anchors = {"mmd_u": {"floor": 0.08, "ceiling": 0.004}, "de_score": {"floor": 0.0, "ceiling": 0.8},
               "de_direction": {"floor": 0.0, "ceiling": 0.8}, "variogram": {"floor": 0.005, "ceiling": 0.0002}}
    at_floor = skills_from_anchors("T1", {"mmd_u": 0.08, "de_score": 0.0, "de_direction": 0.0,
                                          "variogram": 0.005}, anchors)
    at_ceiling = skills_from_anchors("T1", {"mmd_u": 0.004, "de_score": 0.8, "de_direction": 0.8,
                                            "variogram": 0.0002}, anchors)
    assert all(abs(v - 0.5) < 1e-9 for v in at_floor.values())
    assert all(abs(v - 1.0) < 1e-9 for v in at_ceiling.values())
    worse = skills_from_anchors("T1", {"mmd_u": 1.0, "de_score": -0.5, "de_direction": -0.5,
                                       "variogram": 1.0}, anchors)
    assert all(0.0 < v < 0.5 for v in worse.values())


def test_abs_metrics_use_magnitude():
    anchors = {"severity_slope": {"floor": -6.9, "ceiling": -0.01}}
    over = skills_from_anchors("T3", {"severity_slope": 0.5}, anchors)["severity_slope"]
    under = skills_from_anchors("T3", {"severity_slope": -0.5}, anchors)["severity_slope"]
    assert over == pytest.approx(under)
    assert 0.5 < over < 1.0


@pytest.mark.parametrize("task", ["T1", "T2", "T3"])
def test_aggregate_scale(task):
    half = dict.fromkeys(RANKED_METRICS[task], 0.5)
    full = dict.fromkeys(RANKED_METRICS[task], 1.0)
    assert aggregate(task, half)[0] == pytest.approx(50.0)
    assert aggregate(task, full)[0] == pytest.approx(100.0)
    assert sum(w for _, w, _ in GROUPS[task]) == pytest.approx(1.0)
    nan = dict.fromkeys(RANKED_METRICS[task], float("nan"))
    assert aggregate(task, nan)[0] == 0.0


def _stages(task, profiles, rng):
    base, later = profiles
    props_ref, props_tgt = [0.5, 0.3, 0.2], [0.3, 0.3, 0.4]
    spatial = task != "T1"
    centres = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0]], float)
    ref = simulate_stage(4000, base, props_ref, rng, coords_scale=2.0 if spatial else None,
                         coords_shift=centres if spatial else None)
    def target_like(seed):
        r = np.random.default_rng(seed)
        return simulate_stage(4000, later, props_tgt, r, coords_scale=3.0 if spatial else None,
                              coords_shift=1.5 * centres if spatial else None)
    return ref, target_like(1), target_like(2)


@pytest.mark.parametrize("task", ["T1", "T2", "T3"])
def test_floor_scores_fifty_and_good_prediction_beats_it(task, boards, profiles, rng):
    ref, target, good = _stages(task, profiles, rng)
    # The server's floor row is the reference's 10% working copy, which scores 50 by construction.
    working_copy = subsample(ref, SUBSAMPLE_FRAC, seed=1)
    floor = score_board(boards[task], working_copy, target, ref, seed=0)
    assert floor.score == pytest.approx(50.0, abs=1e-6)
    resubmitted = score_board(boards[task], ref, target, ref, seed=0)
    if task == "T2":
        # neighborhood_mmd builds 15-NN neighbourhoods over all submitted cells, so a 10x denser cloud than the
        # truth half shifts it: cell count is a real lever on T2 even though it is not itself ranked.
        assert np.isfinite(resubmitted.score)
    else:
        # Elsewhere resubmitting the full stage differs from the floor only by sampling noise.
        assert resubmitted.score == pytest.approx(50.0, abs=3.0)
    better = score_board(boards[task], good, target, ref, seed=0)
    assert better.score > 75.0, better.group_skill
