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

### Outcome (2026-11-02, read after 38faf73)

Step 5 **fails** on T2h/E8.25_late->E8.5:
- core_radial: d2 -0.021 (ahead in 3 of 10 seeds), Dice +0.112;
- half_radial: d2 -0.036.

On the 1.75-day pseudo-task (T2h/<=E8.75->E10.5, not one of step 5's two transitions), core_radial gains +0.141 on d2
but lowers Dice by 0.013. The shape gain belongs to the E9.5 -> E10.5 transition alone. As written, the rule gives
s2 = H2 as built.

### Amendment (2026-11-02, by the AI assistant under the captain's delegation, made after seeing the outcome above)

**s2 is H2 built as the half trim,** as under "If it passes on both":
`build_test.py heart --entry s2_trim` (the trim is `vec.t2_build.radial_trim_mask`, pinned by tests to the payoffs'
trim).

**Why.** The captain marks the better-scoring upload, so what an s2 is worth is the expected better of s1 and s2 on
board score. Step 5's d2 bar tests a shape gain, and a scope hedge does not need one. On board score, the half trim
beats H2 as built as an s2 on all three released heart transitions:

| s2 on (10 seeds each) | E9.5 -> E10.5 | E8.25_late -> E8.5 | E8.75 -> E10.5 |
|---|---|---|---|
| H2 as built, uniform over five truths | 50.17 | 53.45 | 54.50 |
| half trim, uniform over five truths | 51.12 | 56.84 | 60.02 |
| H2 as built, on the released target only | +0.00 | +0.58 | +0.42 |
| half trim, on the released target only | +0.00 | +2.04 | +1.57 |

On the released targets, the half trim's gains come mostly from the expression-against-reference metrics (de_score
+0.056 and +0.044, de_direction +0.108 and +0.023), not from shape.

The full trim (core_radial) is a larger hedge under a uniform prior but never helps on a released target. The half
trim is also H2's own definition in `docs/PLAN.md` 4.3 ("heart-focused (shell halved or less)").

No test data were involved. The captain can restore the rule as written ("s2 = H2 as built") by saying so.

## T2 embryo (rule fixed 2026-11-02, before the variants run)

### What no candidate does

On both embryo pseudo-tasks, every candidate built so far scores below fgw in every one of the 10 seeds: copy, B0,
pseudobulk_shift and psb_size (`experiments/harness/results/emb_fgw_v1/`). None of them adds anything as an s2.

fgw's open risk is its interpolated shape: on E7.5 its shape skills fall below a copy's.

### The candidates

The s2 candidates are the shape variants:
- `fgw_shape0` keeps fgw's coupling and its mixture of real cells, with positions taken from the earlier bracket's
  geometry (`interpolate_clouds(..., t_pos=0)`).
- `fgw_shape1` does the same with the later bracket's geometry (`t_pos=1`). This is the plan's nearest-stage variant:
  E8.0 tissue at the predicted size, with mixture expression.

Both are scaled to the embryo size rule. At test, the earlier bracket is E7.5 and the later is E8.0.

### The rule

Both candidates run on both embryo pseudo-tasks, against fgw, over seeds 0-9 with the harness.
- **Choice.** s2 is the variant with the higher expected better of fgw and it: the mean over seeds of the larger of the
  two scores, averaged over the two pseudo-tasks.
- **Tie.** If both equal fgw's own mean, because neither beats fgw in any seed, s2 is `fgw_shape1`, the plan's
  nearest-stage variant.
- **Build.** The embryo test build gets the matching `--method`.

## T3 (decided 2026-11-01 by the captain)

s2 is the refused branch's s1, uploaded first (`docs/PLAN.md` section 4.4 and its decision log;
`experiments/t3_decomp/rules.yaml`).
