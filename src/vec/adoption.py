"""Adoption rule for model changes (approved by the captain on 2026-11-01; experiments/adoption_rule/README.md).

A change (a new component, setting or model) is kept only if, on its board's pseudo-tasks scored with the calibrated
emulator over the same seeds (and, where the code allows, the same cell draws for both candidates):
  1. it passes approved rule 7's test, a paired mean gain above 0.3 board points and ahead in at least 8 of the 10
     seeds, on at least 2 pseudo-tasks, or on every pseudo-task when the board has fewer than 3; and
  2. on no pseudo-task does it lose more than 0.3 points on average.
Otherwise the simpler candidate stays. It replaces "improves >= 2 of 3 pseudo-tasks" (ASSUMPTIONS.md A17), under which
a change with no effect passes half the time. Decisions already committed under their own pre-registered rules (size,
T3 decomposition, density and the count floor, the T1 direction readout) are not re-decided with it.
"""

from __future__ import annotations

import math

SEEDS = tuple(range(10))
MIN_GAIN = 0.3  # board points; the same number as approved rule 7 (vec.t2_build.COUNT_RULE_MIN_GAIN)
MIN_AHEAD = 8  # of the 10 seeds; as rule 7 (vec.t2_build.COUNT_RULE_MIN_AHEAD)
MAX_LOSS = 0.3  # board points: a larger mean loss on any pseudo-task vetoes the change


def keep_change(paired: dict[str, dict[int, tuple[float, float]]]) -> dict:
    """Apply the adoption rule. `paired` maps each pseudo-task to {seed: (current candidate's score, the change's
    score)}, both scored on the same seed. Every pseudo-task must have a finite pair for each of SEEDS; otherwise the
    change is not kept and the record says why (as rule 7 keeps max_cells when its rerun is incomplete)."""
    if not paired:
        return {"keep": False, "reason": "no pseudo-tasks", "tasks": {}, "need": None, "passed": 0}
    tasks, incomplete = {}, []
    for name, by_seed in paired.items():
        missing = [s for s in SEEDS if s not in by_seed]
        bad = [s for s in SEEDS if s in by_seed and not all(math.isfinite(float(v)) for v in by_seed[s])]
        if missing or bad:
            incomplete.append(name)
            tasks[name] = {"complete": False, "missing_seeds": missing, "non_finite_seeds": bad}
            continue
        diffs = [float(by_seed[s][1]) - float(by_seed[s][0]) for s in SEEDS]
        mean = sum(diffs) / len(diffs)
        ahead = sum(d > 0 for d in diffs)
        tasks[name] = {"complete": True, "mean_gain": mean, "seeds_ahead": ahead,
                       "passes": mean > MIN_GAIN and ahead >= MIN_AHEAD, "vetoes": mean < -MAX_LOSS}
    need = 2 if len(paired) >= 3 else len(paired)
    if incomplete:
        return {"keep": False, "reason": f"incomplete pseudo-tasks: {sorted(incomplete)}", "tasks": tasks,
                "need": need, "passed": sum(bool(t.get("passes")) for t in tasks.values())}
    passed = sum(t["passes"] for t in tasks.values())
    vetoes = sorted(n for n, t in tasks.items() if t["vetoes"])
    if vetoes:
        reason = f"loses more than {MAX_LOSS:g} points on average on {vetoes}"
    elif passed < need:
        reason = (f"passes rule 7's test (mean gain > {MIN_GAIN:g}, >= {MIN_AHEAD} of {len(SEEDS)} seeds ahead) on "
                  f"{passed} of {len(paired)} pseudo-tasks; needs {need}")
    else:
        reason = f"passes on {passed} of {len(paired)} pseudo-tasks (needs {need}) and loses on none"
    return {"keep": not vetoes and passed >= need, "reason": reason, "tasks": tasks, "need": need, "passed": passed}
