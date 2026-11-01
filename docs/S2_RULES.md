# s2 rules (pre-registration)

**Authorship.** The AI assistant writes these rules under the captain's delegation of 2026-11-02 (`docs/PLAN.md`
decision log). The captain can overrule any of them until the freeze.

**Timing.** Every rule is fixed before any test-board score exists: on 2026-11-02 nothing has been uploaded. At the
freeze (by 2026-11-23, `docs/PLAN.md` section 7), the file's sha256 goes into the decision log and the commit is
pushed. After that a rule changes only with the captain's say-so, logged with its reason.

**No score feedback.** The rules are unconditional: s2 is the same file whatever s1 scores. No test score computes
anything (Official Rules: scores choose, never compute).

**Uploads.** The captain marks the better-scoring upload as final (section 8). Every official upload needs his
explicit go-ahead.

Evidence for each rule: `experiments/s2_hedge/README.md`. These are sensitivity analyses on released transitions, not
adoption results.

## T1 (decided 2026-11-02)

s2 is M1 built exactly as s1, with one difference: **cardiomyocytes at 0.60 of heart cells** instead of E10.5's own
0.265. Blood stays at E10.5's own chamber share, as in s1.

`uv run python experiments/builds/build_test.py t1 --cm-share 0.60 --out <s2 file>` (the seed and count as s1's).

**Why.** On the 2-D grid (T1/E9.5->E10.5:chamber, oracle base, 10 seeds), the cardiomyocyte share is the lever:

| Prior | Gain of the 60% s2 over s1 alone |
|---|---|
| Uniform over 15 tissue mixes | +4.74 |
| A third on own | +3.39 |
| Blood own, cm from own to 0.75 | +6.40 |
| Cardiomyocytes at 0.45 or above | +7.90 |

A blood-share hedge adds 1.5 points at most, and wins only if the cardiomyocyte share is known to stay at E10.5's. This
supersedes `docs/PLAN.md` 4.1's expectation of a blood-share hedge.

## T2 heart (rule fixed 2026-11-02, before the second transition's result)

### The test

The test is `docs/PLAN.md` 4.3 step 5 for the radial trim, read with `experiments/s2_hedge/payoff.py step5` on both
released transitions:
- T2h/<=E9.5->E10.5 (`results/heart_scope/step5.txt`);
- T2h/E8.25_late->E8.5 (runs at dcf1779).

The operation passes if, on both transitions, the trim at full magnitude (`core_radial`: the surrounding, non-core
cells removed farthest in plane from the core centroid first, until the in-plane scope is 0) gains at least 0.15 in
mean d2_shape skill over hold and does not lower mean occupancy_dice skill.

On the first transition it passes (d2 +0.416, Dice +0.041, ahead in 10 of 10 seeds).

### If it passes on both

s2 is H2 built as a trim:
1. Take E10.5's cells and remove the surrounding cells farthest from the core first, until the in-plane scope is half
   of E10.5's. This is step 5's half magnitude, and H2's own definition: "shell halved".
2. Rescale about the centroid to H2's pre-registered full-cloud RMS.
3. Everything else as s1: expression step, count rule and seed rule.

In the payoffs, this candidate (`half_radial`) adds 2.03 points as s2 under a uniform prior over the five truths,
against 1.08 for H2 as built.

**Physical rationale** (step 5's third condition). In both released transitions the later sample is more compact
around the heart core than the earlier sample's whole cloud: E8.25_late includes yolk sac and non-cardiac cells, and
E9.5 is an embryo portion while E10.5 is a dissected heart. A radial trim removes the tissue farthest from the core
first, as a tighter dissection would.

### If it fails on either

s2 is H2 as built: s1's cells and expression at H2's size.

## T2 embryo (open)

The candidates are:
- FGW with an alternative clock;
- the nearest-stage variant: E8.0 tissue shrunk to the predicted size, with mixture expression;
- psb_size.

A hedge analysis on the two embryo pseudo-tasks decides, before the freeze.

## T3 (decided 2026-11-01 by the captain)

s2 is the refused branch's s1, uploaded first (`docs/PLAN.md` section 4.4 and its decision log;
`experiments/t3_decomp/rules.yaml`).
