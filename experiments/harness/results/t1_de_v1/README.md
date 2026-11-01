# T1 de_score audit (2026-11-01)

The D2 candidates on the three T1 pseudo-tasks again (floor, copy, pseudobulk_shift, b0, b0_scope; 10 seeds; bake
protocol, scorer build main; commit ab072da, clean tree), now with each candidate's metric panel kept in the rows, so
de_score's parts can be read. Scores match d2_v1. `de_audit.txt` is the output of
`uv run python experiments/harness/de_audit.py tasks/*.json`.

What it shows (veckit de_score, pinned git cc3859a: score = (raw - chance) / (1 - chance); the floor, which predicts
no change, scores exactly 0):

| Pseudo-task | Truth DE genes (up / down) | chance | copy raw | b0 raw | copy, de_score points | b0, de_score points |
|---|---|---|---|---|---|---|
| E9.5 -> E10.5, chamber-cut truth | 132 / 426 | 0.388 | 0.106 | 0.279 | -2.78 | -1.25 |
| E9.5 -> E10.5, whole truth | 69 / 382 | 0.410 | 0.107 | 0.203 | -3.04 | -2.26 |
| E8.5 -> E9.5 | 31 / 30 | 0.021 | 0.046 | (n/a) | +0.21 | (n/a) |

- The handicap is a property of the E9.5 -> E10.5 change, not of fresh cells as such. There the truth's changed genes
  are mostly genes highly expressed in the reference going down (426 of 558 on the chamber cut), so ranking genes by
  the reference's own expression (the metric's null, sign chosen in hindsight) already names about 40% of them. A copy
  names about 11% (its change is sampling noise) and B0 about 28%: both score below 0, the floor's value, in all 10
  seeds. On E8.5 -> E9.5 the null names 2% and a copy beats it.
- E9.5 is an embryo portion and E10.5 a dissected heart, so this transition mixes development with a change in what
  was sampled. How strongly the E12.5 change follows expression level is unknown; the test transition (E10.5 -> E12.5,
  both from hearts) need not share it.
- A copy's de_direction is above 0 (+0.10 to +0.22): its change shares the reference working copy's sampling noise
  with the truth's change (both are measured against that copy).
