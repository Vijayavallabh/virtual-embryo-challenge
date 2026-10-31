# Adoption rule for model changes (approved by the captain on 2026-11-01)

Raised by a red-team agent on 2026-10-31 and by a side agent on 2026-11-01; checked on 2026-11-01, before any
release-data result has been used to adopt a change (no decision under this rule exists yet).

## The current rule and its problem

ASSUMPTIONS.md A17, `docs/PLAN.md` section 4.1 (the T1 table's "Kept if" column) and feat-015 keep a change if it
"improves" / "beats" the previous candidate, "on at least 2 of 3 pseudo-tasks" for M2. No margin is written anywhere,
and feat-008 pairs the comparison over 10 seeds. A change with no effect is ahead on each pseudo-task with probability
one half whatever the number of seeds, so it passes "2 of 3" half the time (both of two: a quarter).

## Measured noise and simulated adoption rates (`gate_analysis.py`, results in `results/`)

Seed-to-seed sd of a paired score difference in the team's committed results: median 0.95 board points (IQR 0.57-1.52;
384 pairs from the density runs on the validation targets, the T1 pseudo-task on released stages, and the calibration).

| Pseudo-tasks | Rule | No effect kept | -0.5 kept | +0.5 kept | +1 kept | +2 kept |
|---|---|---|---|---|---|---|
| 3 | current | 0.50 | 0.01-0.10 | 0.90-0.99 | 0.99-1.00 | 1.00 |
| 3 | proposed | 0.00-0.04 | 0.00 | 0.32-0.34 | 0.79-0.92 | 1.00 |
| 2 | current | 0.25 | 0.00-0.04 | 0.65-0.91 | 0.92-1.00 | 1.00 |
| 2 | proposed | 0.00-0.03 | 0.00 | 0.14-0.17 | 0.52-0.69 | 0.97-1.00 |

(sd 0.95; ranges cover transitions that agree (tau 0) or disagree (tau 0.5 points) on the true gain.)

## The rule (the approved count rule 7's test, applied per pseudo-task; approved 2026-11-01, "yes for all")

A change is kept only if, on its board's pseudo-tasks scored with the calibrated emulator over the same 10 seeds (and,
where the code allows, the same cell draws for both candidates):
1. it passes rule 7's test, a paired mean gain above 0.3 points and ahead in at least 8 of the 10 seeds, on at least 2
   pseudo-tasks, or on every pseudo-task when the board has fewer than 3; and
2. on no pseudo-task does it lose more than 0.3 points on average.
Otherwise the simpler candidate stays.

It replaces the test in A17, the T1 table's "Kept if" entries (B0 against copy_last, M1 against B0, M2 against M1,
the stretch model) and "M2 replaces s2 only if it wins the backtests", and feat-015's condition; the pseudo-task
harness (feat-008) applies it in code through `vec.adoption.keep_change` (tests: `tests/test_adoption.py`). It changes none of the decisions already committed under their own
pre-registered rules (size, T3 decomposition, density and the count floor, the T1 direction readout), the T2 heart
shape rule (its own margin of 0.15 d2 skill on both transitions) or D3 (the 3-point Generality gap).

Cost: real gains under about half a point are mostly not kept (they are also too small to trust across transitions,
A17); with two pseudo-tasks, a real 1-point gain is kept a little over half the time. Scoring both candidates on the
same cell draws lowers the noise and raises those odds.
