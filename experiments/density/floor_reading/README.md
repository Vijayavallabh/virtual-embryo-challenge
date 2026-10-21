# Count floor rule: how to read "next" (2026-10-21)

Written by the AI assistant (Claude) on 2026-10-21 from the four agents' reports. The user left the reading of "next"
to the AI assistant, and the code applies the judge's verdict.

## Question

Approved rule 7 (`docs/RELEASE_RULES_SUMMARY.md`) can pick 0.25 x max_cells on a validation target. On a test board
with max_cells under 4 x truth_cells, that level gives fewer cells than our floor (truth_cells where published, else
min_cells). The user approved "the next passing count the board allows, or max_cells if none does". Asked how to read
"next", the user replied: "Anything that would improve my chances of winning". Two readings:

- **G, next best gain:** of the passing levels the board allows, the one with the largest measured gain.
- **C, next count up:** of the passing levels the board allows, the one with the fewest cells.

They differ in one case only: 0.25 x wins below the floor, 0.5 x and 0.75 x both pass, and 0.75 x measures at least as
high as 0.5 x. G then submits 0.75 x max_cells and C submits 0.5 x. Choosing by distance to truth_cells was ruled out:
truth_cells stays a bound, never a chooser, while the organisers' answer on using it to choose the count is pending.
C orders the allowed levels by count, not by distance; all of them sit at or above the floor, so truth_cells enters C
only through the same bound that G uses.

## Panel

A workflow of four agents (read-only; released pseudo-target data only; no board files, nothing redacted opened):

| role | verdict | finding |
|---|---|---|
| simulation (these scripts) | C, medium | in the disagreement case, with scores a function of n / truth_cells (x mapping), C is ahead by +0.26 (smooth curves) to +0.40 points (all x families, wiggly curves included), P(C better) 0.82-0.89; absolute-count stress test +0.35 (0.93); fraction-of-max_cells stress test -0.20 (0.38), where G wins through rugged or wiggly level effects; all mappings pooled +0.17 (0.67); the case arises in 0.07-0.35% of releases (x mapping), so the expected gain is about +0.0002 to +0.0014 points per release |
| empirical data | indifferent, medium | all 36 adjacent count steps of the offline imperfect-prediction curves go up (10-20 seeds), so a measured 0.75 x over 0.5 x reversal is noise around a near-tie; C minus G about -0.01 to +0.04 points under its central assumptions (2% prior on real dips, 0.7 chance the validation result carries over), range -0.35 to +0.33 |
| advocate for G | C, medium | C is ahead by about 0.2-0.35 points in the case; a deliberately G-friendly prior only reaches a tie; the cost of being wrong is similar either way (G wrong 0.3-0.6, C wrong 0.1-0.5) |
| judge | **C, medium** | on every curve with at most one peak on which 0.25 x truly wins, true gain(0.5 x) >= gain(0.75 x), so the pattern there is noise around a near-tie (curves that rise with count produce it only when 0.25 x's win is itself noise, and there G is ahead by 0.04-0.13); C's level is the allowed count nearest the validated winner; G is ahead only if the validation result does not carry over (then rule 7 itself would be wrong) or a real dip at 0.5 x carries over (no offline curve shows one) |

Decision: **C**, the judge's verdict, adopted by the AI assistant. In code: `vec.t2_build.COUNT_FLOOR_RULE`, applied by
`vec.t2_build.n_cells_from_decision`. The stake is small: the case is rare, and inside it the difference is a few
tenths of a point on one board.

## Files

- `calibrate.py`: per-seed noise of paired score differences, from `experiments/density/review/results/` and the
  5-seed rows in `experiments/density/results/` (writes `calibration.json`, `calibration.txt`).
- `sim_floor.py`: the Monte Carlo (writes `results.json` and `results_table.md`; about 2 minutes; `--quick` for a short
  run, which also overwrites both). The table also copies in `shift_scan_extended.txt` (an extended shift scan whose
  script was not kept) and `floor_composition.txt` (from `floor_composition.py`: what each reading picks across all
  floor cases) when they are present.
- `verify_rule.py`: checks the simulation's rule 7 and reading G against `vec.t2_build.count_override` and
  `_passing_ranking` on 4,000 random tables (`verify_rule.txt`, a reworded copy of its printed line: 0 mismatches). Its
  docstring describes `n_cells_from_decision` as it stood when the panel ran (reading G). The code now applies C,
  which `tests/test_t2_build.py` checks; a verifier also compared the simulation's reading C with the new code on the
  same 4,000 tables (0 mismatches over 1,038 floor cases; that check was not kept).

Rerun from the repo root, in this order: `uv run python -I experiments/density/floor_reading/calibrate.py`, then
`sim_floor.py`, then `verify_rule.py`. The empirical analyst's, the advocate's and the judge's scripts were not kept;
their numbers above are as they reported them.
