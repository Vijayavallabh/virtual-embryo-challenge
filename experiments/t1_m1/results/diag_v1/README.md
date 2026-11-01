# M1 diagnostic, diag_v1 (2026-11-01)

**Not an adoption result.** The harness refuses Feng on every T1 pseudo-task: Feng's E9.5 cells lie inside each
task's version of the test's external window ((E8.5, E11.5] for E9.5 -> E10.5), though outside the test's own
((E9.5, E13.5]). This is the geometry of the pre-registered T1 direction readout (`experiments/t1_direction/`), whose
"build" verdict set M1's direction source to Feng (rule approved by the captain on 2026-10-21). The readout chose Feng
on this same transition, so the run is also in-sample for that choice: under the adoption rule's tuned-on reading
neither E9.5 -> E10.5 task would count for Feng's direction, and E8.5 -> E9.5 cannot use it, so the harness has no task
left to judge it on.

What ran (`experiments/t1_m1/diagnose.py`, commit 401e3eb, clean tree; DGX): B0 exactly as the harness builds it (its
scores here equal the harness's seed for seed) against B0 with Feng's E9.5 -> E14.5 direction on the cells of Feng's
groups (CM, Endocardium, Endothelium, Epicardium: 4,985 of 5,118 cells), scaled by the step's share of Feng's bracket
(1/5 here; 2/5 at test), held-flat genes left alone; the other cells (blood, mostly) keep B0's damped momentum. Every
setting was fixed before the run; nothing was chosen on its result. 10 seeds, bake protocol, scorer build main,
shared anchors per seed.

| Truth | B0 | B0 + Feng (M1 direction) | M1 - B0 | ahead | M1 against the floor |
|---|---|---|---|---|---|
| E10.5 chamber cut (the test's tissue mix) | 49.68 | 50.97 | +1.29 | 10 of 10 | +0.97, ahead in 10 of 10 |
| whole E10.5 | 44.12 | 45.21 | +1.09 | 10 of 10 | (scope differs; not compared) |

Per-metric skills on the chamber cut (floor 0.5): de_score 0.450 -> 0.480 (raw overlap 0.279 -> 0.347 against the
metric's null at 0.388), de_direction 0.608 -> 0.620, mmd_u 0.493 -> 0.497, variogram 0.423 -> 0.428.

Feng's group sizes behind the directions (E9.5 / E14.5 cells): CM 705 / 1,777; endocardial 70 / 107; endocardial plus
vascular 70 / 158; epicardial 13 / 143 (the readout's minimum is 10). 29,715 of 32,285 T1 genes map to Feng symbols.
