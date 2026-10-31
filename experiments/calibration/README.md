# Emulator calibration (feat-005, assumption A2)

Run on 2026-10-31 on the DGX, after the team unlocked the validation answers (2026-10-28). `calibrate.py` scores the
organisers' reference baselines on all five validation boards under `vec.evaluate`'s `protocol="bake"` (the
validation bake that runbook step 4.6 recovered) and the scorer build `"prefix"`, seeds 0-9, and compares them with
the published reference rows in `reference_rows.json` (captured from the public API on 2026-10-18; extracted from the
raw research notes by a screening agent). Verdict rules were fixed in the script before the comparison.

## Result

- **Anchors: all 66 published validation anchors reproduce exactly at seed 0** (T1:val 8, each T2 board 16, T3:gata4
  10; every ranked metric, expression and shape). The local scorer computes the server's metric values.
- **Reference rows, at seed 0 with the published anchors:** T1 `pseudobulk_shift` 48.58 (published 48.5751), T2
  embryo `pseudobulk_shift` 55.22 (55.22), T3 `shift_transfer` 42.48 (42.4843) and `gene_ko` 42.50 (42.5033): exact to
  the published precision, every skill within 0.04 of the published skill. The organisers built these baselines on
  the 10% working copies (construction `wc`). The floors (`copy_last`, `wt_identity`) score 50 on every board.
- **One row is not reproduced:** T2 heart extrapolation `pseudobulk_shift`, 50.32 at seed 0 against 47.76 (10-seed
  mean 50.58, +2.82; de_score skill 50.4 against 40.9). E8.75 and E9.5 use different annotation schemes (5 of their
  33 / 22 `celltype` labels are shared; `cm_celltype` shares only "Unknown"), so the per-type deltas depend on how the
  organisers matched cell types across the two stages. Ten constructions were tried in two rounds (`variants.py`:
  label choices, base cells, the earlier heart stage); none matches. On 2026-10-16 the organisers wrote on Slack that
  a reference table mapping the labels to a consistent annotation existed but was unpublished; it is the likely
  cause. Unexplained.
- **A2 verdict by the fixed rule: fail** (that row's best construction is 2.82 points off; every other row and every
  anchor passes). What it means: the scorer itself is exact (anchors and four rows), and the miss is in rebuilding
  one baseline whose cell-type matching is not published. Local scores of our own predictions on the validation boards
  (seed 0, bake protocol, published anchors) are the server's scores.

## Provenance

Board runs at 1c8bf05 (clean tree) except `T1:val`, which started at 1c8bf05 and records e7d8230 because the repo was
synced to the DGX while it ran (the old script recorded the commit at the end; fixed in 65a9c29). The two commits
differ only in `summarise()` and the variants script, not in the run path. Variants: round 1 at e7d8230, round 2 at
2df4feb. Files: `results/boards/*.json` (per seed and construction), `results/calibration_summary.{json,txt}`,
`results/variants/`.
