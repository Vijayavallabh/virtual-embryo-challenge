# Gate estimate: our estimated leaderboard total from committed results

`estimate.py` computes "our estimated score in leaderboard" for the captain's upload gate (`docs/PLAN.md` section 8,
`results/upload_gate.json`) from committed results only. **It is aggregation only, not a try.** It reads committed
practice-task results and never a target, a test-board file or a board endpoint; from what it reads it uses only
labels, seeds and board scores (no cell count, anchor or target field is used or printed). It runs no candidate,
builder or scorer and computes no new score: every number it prints is a mean, a maximum or a percentile of board scores
already committed. No test score and no leaderboard figure enters it. The gate's threshold is read from
`results/upload_gate.json` and named only as the captain's threshold.

```bash
uv run --frozen --no-sync python -I -B experiments/gate/estimate.py [--json PATH]
uv run --frozen --no-sync python -I -B -m pytest -q tests/test_gate_estimate.py -p no:cacheprovider
```

It writes nothing unless `--json PATH` is given. The JSON holds the readings, bands, steps, sensitivities, the 20
self-checks and the sha256 of every file read.

## What it reads

All reads go through one door (`Door`), which refuses `boards/`, `panels/`, `docs/research/`, `data/`, `runs/` and any
`index.json`, also when a path resolves elsewhere through a link. Missing files stop the run with a list of what is
missing; a missing per-seed part, candidate or field stops it with the file's path.

- **T1**: the oracle-base grid per seed (`experiments/s2_hedge/results/per_seed/t1-grid/`), checked against
  `t1_grid/payoff.json` (mean[truth][pred]) and `t1_grid/hedge_bown_cown.json` (all five priors); M1 against B0
  (`experiments/t1_m1/results/diag_v1/T1_E9.5_E10.5_chamber.json`), whose B0 is checked against
  `experiments/harness/results/d2_v1/tasks/T1_E9.5_E10.5_chamber.json` seed for seed.
- **Heart**: the heart-scope payoffs per seed (`per_seed/heart-scope/`, `per_seed/heart-scope-e875/`), checked against
  `heart_scope*/hedge_hold.json` (all five priors); the entry's `b0` on d2_v1's `T2h_E9.5_E10.5.json` and
  `T2h_E8.75_E10.5.json`.
- **Embryo**: `experiments/harness/results/emb_shape_v1/tasks/` (fgw, fgw_shape0), whose fgw is checked against
  emb_fgw_v1's seed for seed.
- **T3**: `experiments/t3_decomp/results/gata4/gata4_records.csv` (only the columns target, seed, variant and t3 are
  used), checked against `gata4_summary.csv`; K2a and our WT copy on its fold from
  `experiments/harness/results/t3_k2_v1/tasks/T3_Gata4_Mab21l2.json`.

The synthesis's estimator plan named `experiments/t1_m1/results/d2_v1`. That folder does not exist; d2_v1 lives in
`experiments/harness/results/d2_v1/`.

## The two readings

- **A**: each board at the better of its two uploads, because a board shows a team's best. This is the reading of
  decision log 2026-11-04 (c).
- **B**: the first round alone: T1 s1, heart s1, embryo s1 and T3's s2 (T3's first upload is s2).

Overall = T1 + (heart + embryo) / 2 + T3. A board's outcomes are its committed per-seed scores: one per scenario
(truth, transition, backtest, stand-in or replicate) and seed. A takes the better of s1 and s2 in the same seed, as
`experiments/s2_hedge/payoff.py` `hedge()` does. Each reading is the weighted mean of its board's outcomes. Its band is
the 10th and 90th percentiles: exact on each board's weighted outcomes, and for Overall from 200,000 draws (seed
20261105), one outcome per board, with the boards drawn independently. The band therefore covers the committed seed
spread, the stand-ins' spread and the spread over the truths of the uniform prior. The other stated priors are printed
as a sensitivity and do not enter the band.

## The procedure

Decision log (c) says the outlook rests "on our own evidence (committed backtests and the oracle-base grid, no test
score)". It also says T1 and T3 carry bets that the stand-ins lack. Section 8 lists the evidence as "backtests carried
to the test boards, the oracle-base grids, the s2 hedge payoffs". Beyond that, (c) states no steps. Each rule below
marked *chosen* is the simplest one that fits the committed files. None was fitted to 184 or to any other figure.

- **Prior** (*chosen*). Every payoff uses its first stated prior, uniform over its truths. This is the prior that
  `payoff.py` and `docs/S2_RULES.md` list first.
- **Band and agreement** (*chosen*). The band is the 10th-90th percentiles. "Reproduces the outlook" means within 1
  point on every board.
- **T1**. Both entries' bet (E10.5's chamber mix with the s2 hedge) is priced by the oracle-base grid, with
  s1 = `bown_cown` and s2 = `bown_c0.6`. Under the uniform prior the grid gives A 68.28 and B 63.54.
  - *Chosen* backtest-to-test step: the grid is B0 and the entries are M1, so M1's paired gain over B0 on the same
    backtest (+1.29, diag_v1) is added to every outcome. diag_v1's B0 is d2_v1's B0 seed for seed.
  - *Not discounted*: the grid's oracle base (every prediction starts from the truth's own chamber mix). On the same
    backtest the grid's own/own cell (`t1_grid/payoff.json` mean[bown_cown][bown_cown], 69.92) stands 20.24 above
    d2_v1's `b0` (49.68), the entry's construction at E9.5's mix, in 10 of 10 seeds; that gap includes the blood and
    cardiomyocyte shares, and no committed file separates the other shares' part. The +1.29 is in-sample for the
    readout that chose Feng (`experiments/harness/TRIES.md`:99) and is measured at 1/5 of Feng's bracket (2/5 at
    test). So T1 is the reading's least certain board, and an upper one.
- **Heart** (*chosen*). The heart-scope payoffs (s1 `hold`, s2 `half_radial`, the half trim) are used on the two
  transitions the harness counts. The payoffs scale each prediction to the truth's core size, which is an oracle size.
  The entry's pre-registered size is scored by no backtest (A8). So each transition is shifted by the entry's own
  score on that backtest (d2_v1 `b0`) minus `hold` on the released truth: +0.20 on <=E9.5 -> E10.5 and -3.63 on
  <=E8.75 -> E10.5. The size bet gets no credit, as in (c), which names bets only for T1 and T3. E8.25_late -> E8.5
  has no counted backtest and is left out.
- **Embryo**. The two counted backtests, with s1 fgw and s2 fgw_shape0, read as the frozen s2 rule reads them.
- **T3** (*chosen*). s1's two stand-ins on the Gata4 replicates get even weight: A `xrep:F+depth` (73.33), which
  carries s1's knockout-depth bet, and B `transfer:mab21l2_all` (60.77). s2's stand-in is `blind:G2_bcat_central_d50`
  (50.87). On each replicate, s1 and s2 are paired by seed.
- **K2a, included symmetrically** (*chosen*; decision log 2026-11-04 (d) and 2026-11-05 (b)). Stand-in B was "the one
  reading of a factor carried across knockouts" (Mab21l2 -> Gata4). K2a is the other direction: s1's own Gata4 factor
  scored on Mab21l2 (50.08).
  - The cross-knockout half of s1's weight is split evenly between the two directions: A 1/2, B 1/4, K2a 1/4. A K2a
    above B would have raised T3 by exactly the rule that now lowers it.
  - Each T3 scenario is one truth with an s1 and an s2 stand-in, since a board's better of two compares two uploads on
    the same truth. K2a's fold has no s2 stand-in, so our WT copy on that fold (`copy`, 50.74; T3's s2 is near-WT)
    stands in for s2 there, paired by seed (*chosen*, after the independent review; pairing K2a with s2's stand-in on
    the Gata4 replicates instead would compare scores on different truths and read T3 0.28 higher).
  - So K2a lowers reading A's T3, and moves reading B's T3 only by K2a's fold standing in for a quarter of s2's
    (-0.03).

## Result (inputs as committed at 1435acf, unchanged since)

| | T1 | heart | embryo | T3 | Overall | 10th-90th |
|---|---|---|---|---|---|---|
| A, with K2a | 69.57 | 53.85 | 60.87 | 64.58 | 191.51 | 176.8-204.0 |
| A as of (c), before K2a | 69.57 | 53.85 | 60.87 | 67.05 | 193.98 | 184.0-204.0 |
| B, first round, with K2a | 64.83 | 49.20 | 59.95 | 50.84 | 170.24 | 160.2-179.1 |
| B as of (c), before K2a | 64.83 | 49.20 | 59.95 | 50.87 | 170.27 | 160.1-179.7 |
| outlook of (c) (`docs/PLAN.md`:990) | 63 | 52 | 60 | 65 | 184 | 170-198 |

**The procedure does not reproduce the outlook within 1 point per board.** Reading A as of (c) differs by:

- **T1 +6.57.** It is the grid's better of two (68.28) plus M1's +1.29. (c) does not say how it reached 63. The grid's
  s1 alone under the same prior is 63.54.
- **Heart +1.85.** Under the uniform prior, four of the five scope truths are trims, and the half-trim s2 gains on
  those. Under the released-truth prior, A is 51.72.
- **Embryo +0.87.** The better of two is 60.87 and fgw alone is 59.95.
- **T3 +2.05.** (c) states no weights on the stand-ins. Even weights give 67.05.

With K2a, A's T3 falls by 2.47, to 64.58.

Reading B is 170.24, just above the top of the "roughly 160-170 on unreviewed arithmetic" in decision log 2026-11-05
(b). No reading, and no band's upper edge, reaches the captain's threshold. An independent agent recomputed A as of (c)
and B with its own code from the committed files and matched them to 0.01 on every board (wf_ae139ac6-f35, before the
review's K2a pairing was adopted).

Risks that no committed file measures sit outside every band:
- the grid's oracle base on T1's other chamber shares (+20.24 on the backtest, all shares together);
- the two-day horizons of T1 and the heart at test, against one or 1.75 days on the backtests;
- the heart's pre-registered size;
- E7.75 acquired at E7.5-like depth.

The stand-ins and K2a's fold are not the test.
