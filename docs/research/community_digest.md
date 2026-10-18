# Virtual Embryo Challenge: community method digest (public repos)

Compiled 2026-10-18 from a read of every public VEC repository I could find (GitHub search and the organisers' community register). Every number below comes from a repo file I downloaded and read (paths are given), or from the public leaderboard JSON the landscape agent saved. Nothing was cloned and no data files were downloaded.

Tags used throughout:
- **[V]**: verified in the cited file.
- **[L]**: inferred only from a public leaderboard model label.
- **[U]**: claimed but not checkable from the public repo (code missing, or the description is ambiguous).
- **[I]**: my inference, not the author's.

All skill numbers are on the board scale: floor 50, ceiling 100. "Headline" means the board total.

**Licences decide what we may reuse.**
- `tubager/vec-t2`, `wqty123/virtual-embryo-challenge-methods` and `divyanshv2001/virtualembryo2026` have **no licence: use ideas only, do not copy code.**
- `i-habib/virtual-embryo-community`, `AJK0921/vecbench`, `WaelYTz/vec-leaderboard-insights`, EmbryoMatch, Anky9972, Shashwat-srivastav and random-guy-05 are MIT.
- `xxx12e/vec-community-kit` is MIT (its LICENSE adds a data note, so GitHub shows "NOASSERTION").
- `cadentann/*` is Apache-2.0.
- Rules §10 allows published code with disclosure, and separately forbids code sharing across teams. i-habib is very likely a competitor (the #2 Human team is "Ivan & Arnav"; this is a name inference **[I]**). Re-implementing the trivial operators ourselves is cleaner than vendoring them.

---

## 0. Top 10 actionable recipes, ordered by expected value for the test boards

The test boards are T1 E12.5, T2 embryo E7.75 (interpolation), T2 heart E12.5 (extrapolation) and T3 β-catenin. "Gain" is the measured board change on the validation boards.

| # | Recipe (board) | Measured gain | Source (licence) |
|---|---|---|---|
| 1 | **T2 interpolation stack** (embryo; partly heart). Steps, in order:<br>(a) Place real bracket cells onto ONE anchor cloud: cluster-wise Hungarian matching in PCA-32; each cell's X comes from the left or right anchor (Bernoulli, weight w); its xyz is that of the matched left-anchor cell.<br>(b) Developmental "slice": keep only cells whose within-cluster progress along μR−μL is near the clock.<br>(c) Partial-mean: move the top-16 \|μR−μL\| genes 15% toward the clock mean, within each cluster.<br>(d) Hungarian re-pair of the same X multiset onto fixed xyz, using local neighbourhood fields.<br>(e) Per-cluster quadratic curvature, added only where expression is non-zero ("keep-zeros"). | Embryo: copy_last 51.2 → 69.67 (+18.4).<br>Heart-interp: 52.3 → 69.67 (+17.4).<br>Board #1 at the time was 72.1 / 73.5. **[V]** | tubager/vec-t2 `doc/t2_p2_history.md` (no licence) |
| 2 | **Never shift the global mean. Keep real cells and add per-cluster deltas only on non-zero entries** (T2 heart-extrap; principle holds for T1/T3) | Heart-extrap: real-cell chassis 54.66 → 56.45 with α=0.5 (de 45.5→52.0, mmd +1.1, NFS +1.8). Refinements reach 57.57. A global α=1 shift scored **42.5**. T3 global Δ: variogram 49.8 → 12.2. **[V]** | tubager §39.33–39.55; wqty123 `experiments_t3.en.md` |
| 3 | **T3: mix KO-state cells into the matched E8.75 WT carrier, flatten composition, dilute with real WT** | Global Δ transfer 40.98.<br>mix 20% → 70%: 61.5 → 64.6 (sev 79.5 → 92.3).<br>Cell-type^β resampling (β=0.55): 65.31.<br>Real-WT dilution k=0.65: 66.36.<br>Combined (kbb): 66.61. **[V]**<br>The KO cells were presumably the released Mab21l2-KO E9.5 cells **[U]**. | wqty123 `experiments_t3.en.md` (no licence) |
| 4 | **T3 magnitude calibration**: tune response size by mixing or dilution ratio, never by additive scaling. Severity is forgiving on the current board. | A response 2× too strong still scores **91/100**; 90 is reached at 2.2× (veckit docstring) **[V]**. Leaders score ≥97. Additive α-scaling of the true response wrecks the variogram (0.016 → 0.163 as α goes 0.25 → 2; intuition-lab) **[V]** | veckit 0.1.2 `common/core_metrics.py`; i-habib `intuition-lab/RESULTS.md` |
| 5 | **T2 geometry hygiene**: build the cloud, then rescale isotropically to the target RMS radius. Sample rows randomly, because files are stored in spatial order. Use n ≈ 5000 on interpolation boards. Trim over-spread voxels lightly. | Scale skill 71.1 → 100: **+2.4 headline** (heart-interp); 88.9 → 98.7: +0.8 (embryo).<br>Contiguous row block instead of random sample: **48.94 → 38.6** (heart-extrap).<br>n=10k/17.6k vs 5k: −0.7/−0.5 (NFS falls).<br>Light trim of 0.76% of cells: +0.34 (embryo). **[V]** | wqty123 T2 reports; Shashwat-srivastav `virtualembryo-data-gotchas`; tubager §2.19–2.21, §15.8 |
| 6 | **Composition**: resample real cells to a predicted composition, and never drop a population | Frozen E9.5 composition beat trend extrapolation on heart-extrap: G0 54.66 vs trend 51.49. Each 0.1 of L1 composition shift cost about −1.3 headline. But i-habib's heart pseudo-pairs found composition-trend candidates best (+3.5/+4.9 over copy_last). **The two conflict: decide with the released validation answers.** **[V]** | tubager §5.6; i-habib `val-to-test/README.md` |
| 7 | **T1 cheap levers** (mid-table evidence only):<br>- variance amplification: multiply the deviations of the top-K HVGs by c;<br>- maturity-weighted resampling of real E9.5 cells;<br>- exact CP10k library size;<br>- a program-level (not per-gene) pseudobulk shift. | Global shift 48.85 → variance amplification 50.52 → program shift + K=500, c=1.5: 52.18 → maturity weighting b=0.5: **53.67** (b=0 control 51.37). The T1 board top is 76. **[V]** | wqty123 `experiments_t1.en.md` (no licence) |
| 8 | **Ensembling operators** to combine a DE-strong model with a distribution-strong carrier:<br>- mean graft: X_B + (μ_A − μ_B);<br>- quantile graft;<br>- spatial transplant: Hungarian matching in expression PCA. | Controlled tests: quantile graft kept the better parent on 75% of T1 primary metrics and was never worse than both; mean graft kept 67% (2/6 pairs worse than both on some metric). Random coordinate transfer hurt in 5/6 T2 pairs. **[V]** (150-cell toy data) | i-habib `ensemble-lab/*.md` (MIT) |
| 9 | **P3 two-slot selection**: rank by validation and take the top 2. Treat differences inside the noise band on mmd, variogram, NFS and d2 as ties. Do not trust local d2/scale rankings. | Val→test rank Spearman **0.879** (CI 0.81–0.92), top-2 regret ≈ 0. Per-metric Spearman: d2 **−0.54**, scale −0.21. Row-order noise is about ±0.6 headline per submission; de/dir/sev/scale are noise-free. **[V]** | i-habib `val-to-test`; tubager §11.2; vecbench `results/` |
| 10 | **Legal late-stage external data for E12.5 bracketing** (T1, T2 heart): E14.5+ heart scRNA/spatial (Feng 2022 GSE193346; MOSTA E14.5; GSE282547 E14.5 Visium covers 500/500 panel genes; GSE150817 E16.5). Normalise over the panel first. Align assays using a same-stage (E9.5) pair. | Only label-level evidence:<br>- Druggable, T1 #10 Human, 71.65: "…+feng_b4_4k…" **[L]**<br>- APC_TDC, T1 65.94: Feng2022 E9.5+E14.5 **[L]**<br>- LongWz, T1 61.09 / heart-extrap 59.51: "same-stage-assay-alignment [Feng2022 WT E9.5…]" **[L]**<br>No ablation is published. | landscape `model_names.txt`; i-habib `external-data-catalog`, `panel-normalise` |

**Biggest gap in the public record:** nobody has published a T3 recipe for β-catenin, or a heart-extrapolation geometry recipe. The board top on heart-extrap has d2 = 100 and dice ≈ 58, while every public extrapolation geometry experiment lost points. Section 2.3 has details.

---

## 1. Repositories read

| Repo | What it is | Licence | Value to us |
|---|---|---|---|
| [wqty123/virtual-embryo-challenge-methods](https://github.com/wqty123/virtual-embryo-challenge-methods) | Full methods journal: every submission with per-metric board scores (T1 65 rows, T3 42 rows, all T2 boards). Mid-table team (as of 09-23: T1 53.67, embryo 65.31, heart-interp 65.01, heart-extrap 54.51, T3 66.61). | **none (ideas only)** | High for failure modes and T3/T1 levers |
| [tubager/vec-t2](https://github.com/tubager/vec-t2) | T2 code and a ticket-by-ticket ledger (~240 board submissions with per-metric skills and mechanism notes). Locked 69.67 / 69.67 / 57.57 on 09-21 (rank 19, 11 and 24 on the three boards at the time). | **none (ideas only)** | **Highest for T2**; also noise and scorer mechanics |
| [i-habib/virtual-embryo-community](https://github.com/i-habib/virtual-embryo-community) | scorer-versions, metric-lens, leaderboard-anatomy, val-to-test, history-rescore, count-recovery, panel-normalise, external-data-catalog, data-safety, submission-blender, ensemble-lab, intuition-lab, morph-to-truth, metric-movies, test-phase-switchover | MIT | High for tooling, transfer and noise |
| [xxx12e/vec-community-kit](https://github.com/xxx12e/vec-community-kit) | validator, community baselines, local scoring wrapper (multi-seed band), daily rules watch | MIT | Medium: rules changes, protocol |
| [AJK0921/vecbench](https://github.com/AJK0921/vecbench) | noise calibration, fidelity (local-vs-board Kendall tau), metric "map" | MIT (in `vecbench/LICENSE`) | Medium: noise numbers |
| [WaelYTz/vec-leaderboard-insights](https://github.com/WaelYTz/vec-leaderboard-insights) | daily leaderboard snapshots, headroom per metric | MIT | Low (covered in landscape) |
| [ShaoliZhao/virtual-embryo-EmbryoMatch](https://github.com/ShaoliZhao/virtual-embryo-EmbryoMatch) | curated external-dataset registry and stage-alignment notes, no measured gains | MIT | Low–medium (datasets, legality) |
| [Shashwat-srivastav/virtualembryo-data-gotchas](https://github.com/Shashwat-srivastav/virtualembryo-data-gotchas) | 7 data pitfalls; spatial row order measured on the board | MIT | Medium (one 10-point pitfall) |
| [Anky9972/celltype-harmonization](https://github.com/Anky9972/celltype-harmonization) | cross-stage and cross-modality cell-type maps computed on the real files | MIT | Medium for T1/T2 composition |
| [divyanshv2001/virtualembryo2026](https://github.com/divyanshv2001/virtualembryo2026) | agent-run T1 attempts with board screenshots (best 52.16) | none | Low (negative evidence) |
| random-guy-05/vec-* (~40 repos) | tooling (samplecheck, bootstrap, failure gallery, replicate-range…); no board results found | MIT | Low |
| finnick123/VE-*, ShaoliZhao/VE-* labs | tiny synthetic-only labs (1–7 KB of code), no board results | mixed | Negligible |
| Shashwat-srivastav/VirtualEmbryo-basline | per-cell-type pseudobulk interpolation baseline, no scores | MIT | Low |
| cadentann/vec-scorelens, cadentann/virtual-embryo-community-kit | diagnostics and walkthrough, synthetic examples only | Apache-2.0 | Low |
| aristoteleo/PantheonOS `skills/virtualembryo/SKILL.md` | organisers' agent skill for their knowledge graph and 3D atlas API | BSD-2 | Legality warning (§9) |

Not VEC (despite matching searches): washieuanan/HeisenGAN, qkrdmsghk/PertVCFM, VirtualEmbryo/* (cell-mechanics project), EmilyBlaster/VEC-2026.

ManchesterCityWhy/Statistical_Learning_Term_Project and VictorRagazzi/virtual-embryo are course projects with plans only and no scores.

---

## 2. Task 2 (3D MERFISH) recipes

### 2.1 Facts about the released T2 data that the recipes depend on [V]
Sources: tubager `doc/t2_technical_plan_v2.md` §1, §5.2; `configs/t2.yaml`; data-gotchas.

**Stage sizes and radii** (RMS radius in `spatial_3D` units):

| Stage | Cells × genes | Types | RMS | Principal-axis sd |
|---|---|---|---|---|
| E6.75 | 7,093 × 498 | 18 | 122.7 | 100.9 / 52.2 / 46.3 |
| E7.25 | 13,295 × 498 | 18 | 147.3 | 105.3 / 86.7 / 55.6 |
| E8.0 | 31,671 × 500 | 26 | 326.6 | 226.6 / 202.1 / 120.3 |
| E8.25 | 58,716 × 500 | 33 | 354.1 | 245.3 / 210.4 / 144.9 |
| E8.75 | 24,826 × 500 | 33 | 216.9 | 136.8 / 128.3 / 108.9 |
| E9.5 | 53,742 × 500 | 22 | 335.0 | 244.1 / 178.3 / 144.4 |

**Radius changes come from the imaged field of view, not growth.**
- E8.25 (354) is larger than E8.75 (217). tubager attributes this to the field of view narrowing, not the heart shrinking.
- [redacted 2026-10-19: A]
- [redacted 2026-10-19: A]

**Board cell bounds** (validation): embryo [redacted 2026-10-19: A]–5000, heart-interp 1000–17616, heart-extrap 1000–[redacted 2026-10-19: A], T3 1000–7449.

**Annotation vocabularies differ by stage.**
- E8.25 and E8.75 share 100% of type names.
- E9.5 shares only 5 names with E8.75; only 32.3% of E9.5 cells fall in same-named E8.75 types.
- E8.0 vs E7.25: 51.9% of cells are in same-named types.
- tubager collapses labels into ~14 "working clusters" per setting (table in `t2_technical_plan_v2.md` §6).
  - Heart clusters: CM_V, CM_IFT, CM_OFT (born at E9.5), endo, peri, ncc, phm, pam, se, gut, nt, extra (field-of-view loss, set to 0 for t ≥ 9.5), branch (born), st (born).
  - Embryo clusters add neural (born at E8.0).

**Rows are stored in spatial order.** |Spearman ρ| of row index against a principal axis is 0.82–1.00. The first 5,000 rows cover only 9–77% of the occupied voxels. Always sample randomly. Measured cost of not doing so:
- heart-extrap: contiguous block of E9.5 scored **38.6**; random sample of E9.5 scored **48.94**;
- embryo: an 8,000-row loader cap scored 55.0, all rows 57.61;
- on T3 the block scored *higher*, because a severity validity gate flipped.

Source: `virtualembryo-data-gotchas/README.md` §7.

**Raw validation anchors** (index.json, transcribed in tubager v2 §1.4):

| Board | Metric | Floor | Ceiling |
|---|---|---|---|
| Embryo | neighborhood_mmd | 0.211 | 0.011 |
| Embryo | occupancy_dice | 0.705 | 0.766 |
| Embryo | scale_log_ratio | [redacted 2026-10-19: A] | 0.005 |
| Heart-interp | neighborhood_mmd | 0.057 | 0.0043 |
| Heart-interp | occupancy_dice | 0.675 | 0.880 |
| Heart-interp | scale_log_ratio | [redacted 2026-10-19: A] | −0.012 |
| Heart-extrap | neighborhood_mmd | 0.072 | 0.0004 |
| Heart-extrap | occupancy_dice | 0.775 | 0.947 |
| Heart-extrap | scale_log_ratio | [redacted 2026-10-19: A] | 0.0074 |

### 2.2 Interpolation boards: measured chain (tubager; no licence, ideas only)

**Evidence:**
- [`doc/t2_p2_history.md`](https://github.com/tubager/vec-t2/blob/main/doc/t2_p2_history.md) (summary);
- [`doc/t2_p2_val_2026-09-05.md`](https://github.com/tubager/vec-t2/blob/main/doc/t2_p2_val_2026-09-05.md) (per-ticket ledger, cited by §);
- [`doc/t2_technical_plan_v2.md`](https://github.com/tubager/vec-t2/blob/main/doc/t2_technical_plan_v2.md) (recipes).

| Step | Embryo E7.5 | Heart E8.5 | What moved |
|---|---|---|---|
| copy_last (upload) | 51.2 | 52.3 | — |
| isotropic rescale to target RMS ("scale_copy") | 54.6 | 55.0 | only scale skill: 50.8 → 91.9 (embryo); 50.5 → 82.5 (heart) |
| global pseudobulk Δ ("shift_scale") | 60.2 | 58.1 | de +11.8/+10.0, dir +16.9/+10.8, NFS +7.3/+7.4, **variogram −9.9/−18.3** |
| composition plus per-cluster placement ("full") | 61.93 | 64.58 | — |
| paired real-cell mix, no Δ (two clouds overlaid) | 65.12 | 61.76 | variogram 41.5 → 66.9; but dice 58 → 40 (embryo), d2 62 → 38 (heart) |
| OT lerp of expression | 62.46 | — | lerp washes out co-expression (variogram 58.7, mmd 55.6) |
| **OT-place** (§2.9–2.11) | **66.18** (w=1/3) | 65.31 (w=0.5) → **68.33 (w=0.25)** | variogram/mmd recovered; NFS 72.1 at w=0.25 |
| + TPS deformation of the left cloud by the clock (embryo only) | 66.32 | TPS failed (SDD 63.6 → 43.4) | removing TPS later cost −1.36 (dice −12.5) |
| n = 5000 + de-duplicated rows | 66.41 | n=5k best (10k: 67.67, 17.6k: 67.85) | dice 56.3 → 59.8 |
| row-order "coordinate transplant" (heart; see §2.5) | — | 68.94 | d2 69.7 → 72.6, dice 43.5 → 46.8 (deterministic) |
| **developmental slice X** frozen onto the old xyz (§15) | 66.77 | not used: on the W2 proxy raw de fell 0.402 → 0.359 | **de 54.6 → 59.5**, NFS +1.6, mmd −2.7 |
| span-trim of over-spread voxels (§14.4, §15.8) | 67.11 | heart version −0.19 (d2 −6.6) | dice 59.8 → 63.9 |
| spatially stratified cloud (keep=0.4, 6³ bins) + band X (§18) | 67.58 | — | d2 +2.5, dice +2.1 |
| **partial-mean** k=16, w=0.15 (§25.5) | 68.26 | 69.10 | NFS +2.0, variogram +1.6 |
| **Hungarian re-pairing** of the same X onto fixed xyz (§39) | 68.89 (k=8) → 69.41 (k=32) → **69.64 (k=256)** | 69.13 (k=8) → **69.67 (k=32)**; k=256 69.33 | NFS +5.4 embryo / +1.7 heart; de unchanged |
| **per-cluster quadratic−linear curvature, keep-zeros** (§39.30) | **69.67** | heart, α=2.5: 69.53 (failed) | embryo de 59.5 → 60.2 |

**How to re-implement each step (ideas only):**

**OT-place** (`t2_technical_plan_v2.md` §5.1; code `src/t2/ot.py::ot_interpolate_alloc`, x_mode=pick, xyz_mode=left).
1. Map labels to working clusters.
2. Composition: π(t*) = (1−w)·π_L + w·π_R. This is *arithmetic*. Log-linear was better offline (−31% mmd) but cost **de −3.6** on the heart board (66.64, §11.1), so it was closed.
3. Inside each cluster, Hungarian-match left and right cells in a panel PCA-32 (squared Euclidean; ≤800 pairs per cluster).
4. Draw n_R ~ Binomial(n_k, w).
5. A right pick gets X = X_R[i_R] and xyz = xyz_L[i_L]; a left pick gets X_L and xyz_L.
6. Clusters born in the interval are copied from the right anchor.
7. No global Δ. Jitter = 0.15 × median kNN distance.
8. Rescale isotropically to the target RMS.

Weights: the embryo is best at the clock weight w = 1/3. The heart is best at **w = 0.25**, not its clock 0.5 (w=0.20: 67.77).

Never:
- lerp expression;
- overlay two clouds;
- use right-anchor geometry for heart (53.08: SDD 62 → 31, NFS 70 → 42);
- attach right-anchor X to random left xyz (NFS −4.5).

**Developmental slice** (`--ot-x slice --slice-mode band --slice-keep 0.25`).
- In each cluster, project cells on the axis (μ_R − μ_L). Progress 0 is the left mean and 1 is the right mean.
- Bernoulli-sample real cells only from the band nearest the clock w (keep fraction 0.25).
- In the native slice geometry, dice collapsed (59.8 → 42.3). The fix was to freeze the slice X onto the previous best xyz sequence.
- On the heart W2 proxy (labels not aligned) the slice *lowered* de, so it was not used for heart.

**Spatially stratified cloud.**
- Bin the left anchor's canonical coordinates into 6³ bins per cluster and keep 0.4 of each bin by the band rule.
- Template quotas follow whole-cloud occupancy.
- Use this as the coordinate cloud and Hungarian-assign the band X onto it.

**Partial-mean** (`scripts/30_cluster_moment_morph.py --partial-genes`).
- In each cluster, move only the k genes with the largest |μ_R − μ_L| a fraction w toward the clock mean μ*. All other genes stay untouched, which preserves the co-expression "shell".
- k=16, w=0.15 was best. k ≥ 128 or w ≥ 0.4 hurt NFS and the variogram.

**Hungarian re-pairing.**
- Keep the 5000 X vectors and the xyz fixed. Only re-assign which X sits at which position, by Hungarian matching against a "local field" built from pooled kNN anchor neighbourhoods (k = 8/32/256; w = 1/3 embryo, 0.5 heart).
- Expression and shape metrics are unchanged by construction; only NFS moves.
- The exact local-field construction is described only in ledger §39; I could not find the script in the repo **[U]**.

**Cluster curvature** (`scripts/75_cluster_quad_curv.py`).
- Per cluster, compute μ_quad(t*) − μ_lin(t*) from a 3-knot Lagrange fit (E6.75 / E7.25 / E8.0).
- Add it to each cell, **only on non-zero entries**; original zeros are restored. α = 1; mae 0.012; about 150 genes with |Δμ| > 0.05.
- The same operator applied globally on 29 genes left de unchanged and cost −2.2 on the variogram (ticket 201).

**Interpolation failures to avoid** (all measured on the board):
- panel OT-CFM / joint CFM / PCA-32 generative X: **variogram explodes 4–6×** on both truth windows (§13);
- residual VAE, copula, gene-program and count-generative families: anti-transfer or noise;
- PCA-OT interpolation of real cells: 65.12 (−4.5); bidirectional hop-CFM transport: 64.24 (−5.4). These produce cells that "look like E7.25, not E7.5";
- voxel/uniform coordinate resampling ("occ"): d2 63.6 → 36.3;
- field-of-view crop: **58.40** (−9.9: mmd and variogram drop to the floor);
- occupancy blending: d2 72.6 → 61.2.

**wqty123's T2 interpolation journal** ([`experiments_t2_embryo_interp.en.md`](https://github.com/wqty123/virtual-embryo-challenge-methods/blob/main/experiments_t2_embryo_interp.en.md), [`experiments_t2_heart_interp.en.md`](https://github.com/wqty123/virtual-embryo-challenge-methods/blob/main/experiments_t2_heart_interp.en.md)) independently confirms several points:
- scale calibration alone gave +0.82 (embryo, 88.9 → 98.7) and **+2.4** (heart, 71.1 → 100);
- a type-level global shift destroyed shape (d2 23.8, dice 17.6);
- transplanting a displacement field crushed NFS (57.8 → 41.5);
- variance amplification (top-8 genes × 1.25) helps T1 but **hurts T2** (−0.56);
- their best runs were "homologous pairing + positional interpolation pos_f=0.15" (65.31, embryo) and "mixture + bracket + full-mark scale" (65.01, heart).

### 2.3 Heart extrapolation (E10.5 validation; E12.5 test)

**Evidence:** tubager ledger §2.17–2.22, §5.6, §39.33–39.65, and [`scripts/78_extrap_g0_clust_delta.py`](https://github.com/tubager/vec-t2/blob/main/scripts/78_extrap_g0_clust_delta.py); wqty123 [`experiments_t2_heart_extrap.en.md`](https://github.com/wqty123/virtual-embryo-challenge-methods/blob/main/experiments_t2_heart_extrap.en.md).

| Step (tubager) | Headline | Moved |
|---|---|---|
| copy_last E9.5 | 49.2 | — |
| E9.5 isotropically rescaled to RMS [redacted 2026-10-19: A] | 53.3 → 54.15 (n=5000) | scale 49.9 → 100 |
| global α=1 Δ | **42.5** | mmd, variogram and NFS collapse |
| **G0**: real E9.5 cells resampled by the *frozen* E9.5 composition, n=[redacted 2026-10-19: A], zero for field-of-view-lost clusters, each cell keeps its own xyz (rescaled) | 54.66 | d2 49.5, mmd 51.4, NFS 52.3 |
| composition trend ×0.5 (instead of frozen) | 51.49 | mmd −5.8, d2 −8.3, NFS −3.2 (composition also moves the point cloud) |
| **per-cluster keep-zeros Δ**: X += α(μ_E9.5,k − μ_E8.75,k) on non-zero entries only, α=0.5 | **56.45** | **de 45.5 → 52.0**, dir +2.2, mmd +1.1, NFS +1.8; shape locked |
| α = 0.75 | 56.58 | de saturates at 52.0 |
| + born clusters (CM_OFT, branch, st) get μ_E9.5[born] − μ_E8.75[progenitor] | 56.89 | mmd +1.1 |
| drop Δ for "descending lineages" (gut, ncc, phm, se) | 57.09 | de 52.0 → 52.4 |
| also drop endo, then CM_V | 57.16 → 57.34 | de → 53.3 |
| carriers (CM_IFT, CM_OFT, peri, pam, branch, st) at α=0.95 | 57.51 | de 53.7 |
| removing peri / CM_IFT from carriers | 56.81 / 56.72 | de −1.7 / −2.9 (they *carry* the DE signal) |
| non-carriers receive ε=0.20 × carrier-weighted Δμ ("projection"); finally ncc only | 57.56 → **57.57** | de 54.2 |
| composition boost 1.10 / shrink toward cluster mean / hop-CFM continuation | 56.24 / 56.87 / **41.71** | d2 −8.7 / mmd −4.1 / variogram, mmd and NFS collapse |

**Geometry on extrapolation was frozen throughout.** wqty123 ran four morphology experiments:

| Family | Headline | d2 skill |
|---|---|---|
| r6b70a20 | 38.97 | 11.9 |
| comp_c100 | 32.8 | — |
| vmorph_s120 | 46.8 | — |
| type-level displacement (cpfx) | 52.12 | **5.7** |

Their best was **scale_only 54.51** (E9.5 with coordinates [redacted 2026-10-19: A]; scale skill 100; every other metric ≈ 50).

The current board top on heart-extrap is Omazili Manbo 67.31 (d2 100.0, dice 58.5, NFS 61.6). MorphoFlux111 also has d2 100 (landscape §1.5). So **some teams do get geometry right on extrapolation, with an unpublished method.** Public evidence only says that naive deformation, crop, voxel resampling and per-type displacement all hurt.

**Heart extrapolation: other measured facts.**
- Expression across modalities: the T1 scRNA E9.5 file gives no transferable next-step direction for T2. Cross-modal correlation of 0.38–0.40 comes only from overlapping windows and falls to 0.015 for non-overlapping ones (§8.2(d)).
- A dataset-specific per-gene offset dominates the pseudobulk error. On the W2 truth window (E8.25 + E9.5 → E8.75), even oracle composition leaves pseudobulk relative error 0.277, against 0.324 for the pipeline (§11.4).
- Oracle pseudobulk would be worth about +14 headline. Every *estimable* pseudobulk target (clock-linear, curvature, variance calibration) was worse than verbatim real cells (λ-grids are monotone down; §11.4, §12.2–12.3).

### 2.4 Metric mechanics that the recipes exploit [V]

**Pseudobulk-only vs cell-level metrics.** de, dir (and severity on T3) depend only on pb(pred) − pb(ref), a 500- or 32,285-dimensional vector. mmd and variogram depend on the cells. Shape metrics depend only on the coordinates. NFS depends on the joint (X, xyz) coupling.
- Consequence: X and xyz can be optimised separately, and moves that keep the X multiset fixed leave de/dir/mmd/variogram unchanged by construction.
- Sources: tubager §8.1; wqty123 `reference-links.en.md` §VII "Structure theorem"; i-habib intuition-lab.

**NFS gap is in X, not coordinates.** On the W2 window:

| Combination | NFS |
|---|---|
| True cells (ceiling) | 0.0030 |
| True X + our coordinates | 0.0064 |
| Our X + our coordinates | 0.0706 |
| Our X + true coordinates | 0.1048 |

Replacing duplicated rows with fresh same-cluster cells raised NFS 0.0706 → 0.098 (it breaks the native "this cell sits here" correspondence). Shuffling cells within a cluster costs 2.8× NFS. Source: tubager §11.5, §12.5–12.6.

**NFS sensitivity curves** (i-habib [`metric-movies/README.md`](https://github.com/i-habib/virtual-embryo-community/blob/main/metric-movies/README.md); E9.5 truth, heart-interp anchors):
- NFS keeps full points until expression swaps reach ~6% of the cloud diameter; half its points are lost at 25%;
- position jitter up to ~40 µm is free;
- d2 loses 2.1 points at a 5% volume-preserving stretch;
- removing a slab hurts dice (0.78 at 50% removed) far more than random thinning (0.93);
- isotropic scale changes only scale_log_ratio;
- after the issue #7 fix, rotations change nothing.

**occupancy_dice is partly a count.** Dice = 2|A∩B| / (|A|+|B|) over occupied 16³ voxels, so a cloud that occupies too many voxels loses dice regardless of shape.

| Board | Occupied voxels: ours | Truth (anchor estimate) |
|---|---|---|
| Heart-interp | 289 | ~248 |
| Embryo | 204 | ~181 |
| Heart-extrap | 251 | ~250 |

- The local proxy transfer rate for dice was 0.95 on one calibration and 0.49 when realised.
- For **d2 the local proxy is anti-correlated with the board: −1.30** (local better = board worse).
- Source: tubager §8.3, §14.2.

**scale_log_ratio window.** It is ranked on |value| with a two-sided clip window of T·e^(±|ceiling|) (embryo ±0.53%). Being "below the true radius" is **not** free: [redacted 2026-10-19: A] scored only 91.9 (tubager §14.1). The legality of how a target radius is chosen is discussed in §9.

**Rotation / handedness.** With the issue #7 fix (veckit `main` cc3859a), rotations are score-neutral and mirror images are penalised.
- Without the fix (PyPI 0.1.2 and older), 5/12 random rotations of a perfect E9.5 cloud dropped dice 1.0 → 0.83, worth up to 1.7 (heart-interp) or 3.4 (heart-extrap) points.
- The board's scorer version is unpublished.
- Source: [`scorer-versions/README.md`](https://github.com/i-habib/virtual-embryo-community/blob/main/scorer-versions/README.md), intuition-lab.

### 2.5 Scorer-sampling artefacts (low value in P3, listed for completeness)

**What the scorer subsamples.** In the version tubager read, mmd and NFS take the same `rng(0).choice(n, 2000)` rows, the variogram takes `choice(n, 1500)`, and d2 takes `integers(0, n, 200000)` index pairs. All of this depends on row order. de, dir and scale use all rows.

**Two row-order tricks, both small.**
- "Row-order denoising" gave +0.02.
- Transplanting the coordinate *sequence* of an earlier, better-scored file deterministically reproduced its d2 72.6 / dice 46.8 (+0.51).
- This is selection among already-scored files, which is allowed in principle, but it is useless with 2 test slots.
- The organisers told xxx12e that "the subsample seed will depend on each submission" ([vec-community-kit README](https://github.com/xxx12e/vec-community-kit/blob/main/README.md), 2026.10.0 notes). Expect this to stop working.

---

## 3. Task 1 (scRNA, E12.5 test) recipes

### 3.1 wqty123 T1 journal (no licence; mid-table, best 53.67)
Evidence: [`experiments_t1.en.md`](https://github.com/wqty123/virtual-embryo-challenge-methods/blob/main/experiments_t1.en.md), [`pitfalls.en.md`](https://github.com/wqty123/virtual-embryo-challenge-methods/blob/main/pitfalls.en.md).

**A team-uploaded copy_last does not score 50.** Uploaded copies of E9.5 scored 46.56–46.84 (de 37.7–38.5, dir 51.3–52.4, mmd 48.1–49.1, variogram 48.2–48.5).
- Mechanism: the floor row is built from the scorer's own 10% subsample of the reference, and your copy is a different sample. tubager §2.1 calls 49–52 for an uploaded copy "subsampling noise"; the xxx12e tutorial says the same. That the DE nulls specifically score this small noisy dp below zero is **[I]**.
- The same happens on T3: a wt_identity upload scored 45.81 (de 38.8, dir 46.4).
- **Implication [I]:** a prediction whose true change is small is penalised for its noise on de/dir. You need a real, coherent dp.

**Measured families:**

| Family | Exact operation | Best | Notes |
|---|---|---|---|
| pseudobulk shift | X̂ = clip(X_E9.5 + d·(μ_E9.5 − μ_E8.5), 0); d ∈ {0.25…2.5} | 48.85 (d=2.5: de 45.3, dir 56.4, mmd 55.4, **variogram 34.1**) | de caps near 45; asymmetric up/down variants worse; adding noise to the shift: 42.29 (variogram −19) |
| WOT/OT mixture | type-transition resampling from a WOT coupling | 47.52 (variogram 28.8) | — |
| Markov flux | apply the E8.5→E9.5 transition rule again | **44.19** | time-homogeneity rejected; a local proxy had said 0.964 |
| program mixtures | linear combinations of gene programs (with or without shift) | 45.22–47.34 (variogram 31.7) | convex-hull deadlock |
| low-rank multiplicative | X·exp(0.25·rank-8 factor) | 48.27 | — |
| hand-made selection (by "day", PC1, …) | resample subsets by rules | **30.21–47.02** | mmd and variogram collapse |
| mix 20% E8.5 cells | — | 46.37 | old-stage cells hurt |
| library-size normalisation ("mm") | renormalise so expm1(X) sums to 1e4 per cell; n=1249 | 49.93 | necessary fix |
| **variance amplification ("vs")** | for the K most-variable genes, multiply each cell's deviation from the gene mean by c (K=200, c=1.65) | **50.52** (mmd 55.7) | K=500, c=1.35: 50.01 |
| **knowledge-program shift + variance amplification ("kp3")** | shift the pseudobulk along 10 textbook programs (271 genes; a=0.25, weighted by program hits), then K=500, c=1.5; n=2000 | 52.18 (dir 55.7, mmd 56.7, variogram 54.2) | "raw" variants 47–48 (variogram 46.6) |
| **maturity-weighted resampling ("kp7")** | as kp3, plus resample real cells by a maturity/direction score with exponent b=0.5 | **53.67** (de 44.2, dir 57.6, mmd 57.6, variogram 54.7) | b=0 control 51.37 (de −4.1, dir −3.2, mmd −1.7). The score definition is not published **[U]** |
| state population dynamics | type-level birth/death momentum + program-modulated velocities | 47.93 | — |

Rejected ideas listed in `reference-links.en.md`:
- the organisers' KG (gene sets) has an "empty intersection with the true DE set" **[U]**;
- TF-Sapiens embeddings were anti-correlated (−0.69) as a T1 criterion;
- moscot/MIOFlow structure extrapolation, scVelo/Dynamo/CellRank: no spliced/unspliced layers;
- Harmony/scVI: would erase temporal change.

### 3.2 T1 facts and tools from other repos

**count-recovery** (i-habib, MIT; [`count-recovery/README.md`](https://github.com/i-habib/virtual-embryo-community/blob/main/count-recovery/README.md)).
- The T1 log1p(CP10k) files invert exactly to integer UMIs: divide expm1 by each cell's smallest non-zero value.
- All 17,057 E9.5 cells recover: library sizes min 731, median 11,080, max 94,591. E8.5 library sizes range 1,299–63,850.
- Normalisation was over all 32,285 genes.
- **Use [I]:** generate in count space (binomial thinning, Poisson resampling, NB models, scVI), then re-normalise. This yields realistic per-gene sparsity for the 32k-gene variogram.

**Cell-type harmonisation** (Anky9972, MIT; [`analysis/REPORT.md`](https://github.com/Anky9972/celltype-harmonization/blob/main/analysis/REPORT.md)), on the real files.
- Raw label turnover overstates emergence 2–3×. From E8.75 to E9.5, 85% of label names vanish, but only 8/22 (≈36%) of E9.5 populations lack a confident E8.75 antecedent (pseudobulk correlation ≥ 0.75; permutation p ≈ 0.01).
- Confident renames include Endo → Chamber-Endo (r 0.914) and IFT-CM → A-CM (0.885). The cm_celltype sub-labels renumber: vCM1 → {vCM3, vCM4}.
- Cross-modality map from T1 scRNA E9.5 to MERFISH heart E9.5: V-CM ↔ V-CM 0.54, BEC 0.60, Hepatocyte 0.66, Endocardium ↔ Chamber-Endo 0.60, IFT-CM ↔ A-CM 0.54, OFT/RV-CM ↔ OFT-CM 0.52.
- The organisers' own wording (verified in `site_text/challenge_tasks_temporal.txt`): about half of the next-stage cells have no same-named predecessor.

**Mean graft** (i-habib submission-blender, MIT).
- X_out = X_B + (μ_A − μ_B): A supplies the pseudobulk, B supplies the cell-level structure. The CLI requires `--clip-min 0` and reports the mean error introduced by clipping.
- **Combine with keep-zeros [I]:** tubager's on-board evidence is that adding Δ only on non-zero entries avoids the variogram collapse.

**Failed T1 agent attempts** (divyanshv2001; no licence). Extended Mouse Atlas (Development 2024 dev201867, 25,963-cell cohort) plus a CNF/kinetics model scored 51.87–52.16 (de 44.9–47.2). Negative evidence that atlas pre-training alone does not move T1.

**External late-stage data on the T1 board** (leaderboard model labels only **[L]**):

| Team | T1 | Label |
|---|---|---|
| Druggable (Human #10) | 71.65 | `t1-sel_t05h3k1fdc1000+feng_b4_4k_sprpb1.6` |
| APC_TDC | 65.94 | `t1mean_Feng2022GSE193346_E95E145_support_s1_k10fn` |
| LongWz | 61.09 | `[Feng2022 WT E9.5+E14.5; MOSTA E14.5]` |
| BTW | 58.61 | `s_feng_amp` |
| Leslie Lab | 52.71 | `GSE230531 E8.5+E14.5` |
| bingo / oxtium-ai | 55.16 / 55.02 | `GSE126128 WT/time` |

Source: `SCRATCH/landscape_work/api/lbh_T1_val.json`. These labels show that legal E14.5 data is used up to the T1 top-10. **The landscape's statement that teams using external data are only mid-table needs updating.**

**The organisers' framing for T1** (temporal task page): clearing the floor "plausibly needs supervision across the whole stage ladder, a non-autonomous field, and a treatment of populations appearing". No public T1 method above 55 discloses how it does this.

### 3.3 De_score levers distilled from all repos [V, with I where marked]
- de_score is a signed top-n overlap against a chance level built from ±pb(ref). Global or monotone transforms of dp do not change it. You score only by naming the right genes in the right direction.
- Global shifts cap T1 de near 45 skill (wqty123). In T2, global Δμ left the between-cluster relative DE unchanged, so de stayed fixed (tubager ticket 201).
- On the board, de moved only through operations that change the cell-type-specific means of the real population:
  - per-cluster deltas: heart-extrap +6.5;
  - developmental slicing: embryo +4.9;
  - per-cluster curvature: +0.7;
  - cell-type^β resampling: T3 +3.0;
  - maturity-weighted resampling: T1 +4.1.
- **Implication [I]:** de comes from *which cells are present and their within-type state*, more than from gene-wise extrapolation.
- A do-nothing file could once score de between 0.000 and 0.358 purely through normalisation choice (veckit docstring, `library_size_ratio`). The current expression-magnitude null was built against that attack. Always submit exact log1p(CP10k): expm1 row sums equal 1e4 (wqty123/vec-scalecheck; xxx12e/gh-dv-openclaw validators).

---

## 4. Task 3 (perturbation; β-catenin test) recipes

### 4.1 wqty123 T3 journal (no licence; best 66.61; board top 80.27)
Evidence: [`experiments_t3.en.md`](https://github.com/wqty123/virtual-embryo-challenge-methods/blob/main/experiments_t3.en.md) (the Chinese version, `experiments_t3.md`, has the same wording).

| Family | Operation | Headline (de / dir / sev / mmd / variogram) |
|---|---|---|
| wt_identity upload | copy E8.75 WT | 45.81 (38.8 / 46.4 / 50 / 50.5 / 49.8) |
| global Δ transfer ×1.5 | add pb(Mab21l2 KO E9.5) − pb(WT E9.5) to the E8.75 carrier | **40.98** (38.8 / 48.4 / 50 / 31.2 / **12.2**) |
| cardiac-restricted transfer | Δ applied only to cardiac cells | 44.36–44.85 (de 37.4–37.7) |
| literature prior (Mab21l2 targets) / "Gata4 state program" | — | 45.78–45.80 (zero gain) |
| **mix** p% "real KO cells" + (100−p)% WT carrier | p = 20/30/40/50/60/70 | 61.5 / 62.37 / 63.58 / 64.37 / 64.51 / **64.60** (sev 79.5 → 92.3; de 51.3–53.3; mmd 52.3–58.8) |
| mix 40% with an early carrier | — | 50.09 |
| per-gene transforms on named genes (pk26, ctlp26, kodir, koq, mx50amp) | — | 38.9–58.0; **variogram 10.8–27.1 in 5/5 cases** |
| **cmp**: resample the KO pool by cell type with probability ∝ count^β | β = 1.0 / 0.70 / **0.55** (n=6000) | 63.35 / 64.93 / **65.31** (de 56.3, sev 97.1; mmd 44.2, variogram 47.1) |
| cmpw: 70% flattened KO pool + 30% E8.75 WT carrier | — | 65.49 (de 54.8, sev 92.4, **mmd 55.2, variogram 53.9**) |
| official population_reweight + GSE208162 prior | — | 54.49 (sev 70.6) |
| **kb**: KO pool diluted with real WT, k = 0.65 | — | **66.36** (56.3 / 62.2 / 97.9 / 45.0 / 50.5); k=0.75: 65.32 |
| **kbb**: kb (k=0.65) × cmp (β=0.55) | — | **66.61** (56.3 / 62.2 / 98.0 / 46.1 / 51.8) |

**Reading the T3 results.**
- The "real KO cells" are not named. In P2 the only released KO was Mab21l2 at E9.5, so they are presumably those cells **[U]**.
- The carrier must be the matched E8.75 WT. The authors' local A/B "falsification" of kb had wrongly used the E9.5 WT as carrier, which is why it disagreed with the board.
- **[I]:** [redacted 2026-10-19: B]
- With Gata4 answers released at P3 this can be tested directly. Its relevance to β-catenin is unknown.
- Priority (pitfalls §17): de (0.30) + severity (0.25) outweigh mmd + variogram (0.20), so hitting the effect beats distribution fidelity on T3. Dilution with real WT then buys back mmd and the variogram.
- "Displacement + zero-inflation mask": `np.where(X>0, X+d, X)` instead of `X+d` (`reference-links.en.md` §VII.4). This is the same keep-zeros idea as tubager's.

### 4.2 Magnitude calibration (severity_slope) [V]
- **Current board metric:** log of the zero-intercept OLS slope of dp on dt over the truth DE genes. A worst value of log(1e-3) applies if β ≤ 0 or R² ≤ 0.01. veckit's own docstring (0.1.2 `common/core_metrics.py::severity_ratio`): measured on the T3 Gata4 board, a response **2× too strong still scores 91/100**, and 90 is only reached at 2.2×. Overshooting is cheap; an inverted or unrelated response hits the floor.
- **intuition-lab, Mab21l2 mini-sample** ([`intuition-lab/RESULTS.md`](https://github.com/i-habib/virtual-embryo-community/blob/main/intuition-lab/RESULTS.md)). Scaling the true response by α gives severity −1.49 / −0.80 / **−0.106** / +0.30 / +0.59 at α = 0.25 / 0.5 / 1 / 1.5 / 2. de and dir are flat for any α > 0. The variogram goes 0.0158 / 0.0285 / 0.0668 / 0.1128 / 0.1629 because the additive global response damages co-expression. A reversed response gives de −0.13 and dir −0.95.
- **On the board** (landscape §1.5): leaders have sev ≥ 97 (|log β| ≈ 0.04); the spread is in de (top 74.1) and dir (top 74.0).
- **Hedge [I]:** PR #6 (shipped in PyPI 0.1.2, not on the board) would replace this with min(β, 1/β). Then 2× too strong becomes 0.5 and magnitude would matter much more. Aim for β ≈ 1, not "anything ≥ 1".
- **Replicates:** the T3 board scores a submission against both KO replicates and publishes the score difference ([data page](https://virtualembryo.ai/challenge/data), verified in `site_text/challenge_data.txt` line 241). This gives a direct biology-noise reading for the β-catenin test.

### 4.3 β-catenin-specific community content
None found. No repo has a Ctnnb1/Wnt recipe.
- The landscape and brief note that Ctnnb1 and classic Wnt targets are off-panel.
- [redacted 2026-10-19: B]
- i-habib catalogues the Mouse GO GAF (CC BY 4.0; 1,111 GO edges for Gata4/Ctnnb1/Mab21l2) as a non-measured prior.

---

## 5. Cross-cutting: what NOT to do (measured)

| Idea | Board result | Source |
|---|---|---|
| OT-CFM / joint flows / VAEs / diffusion generating panel X | variogram 4–6× worse; "never stably exceed copy-last" | tubager §13, §25; wqty123 pitfalls §14 |
| Adding noise or scatter to predictions | variogram −19, mmd down | wqty123 T1 |
| Lerping expression between stages | variogram and mmd lost | tubager §2.7 |
| Overlaying two anchor clouds | dice 58 → 40 (embryo), SDD 62 → 38 (heart) | tubager v2 §3 |
| Right-anchor geometry on heart interp | 53.08 (SDD 31, NFS 42) | tubager v2 §3 |
| Voxel/uniform resampling, field-of-view crop, slab removal | d2 −27 to −33 | tubager §2.17, §5.2 |
| Any deformation on extrapolation | 32.8–52.12 | wqty123 heart-extrap |
| Global α=1 shift (extrapolation) | 42.5 | tubager §2.3 |
| Composition trend on heart-extrap | 51.49 (−3.17) | tubager §5.6 |
| Per-gene transforms on T3 | variogram 10.8–27.1 | wqty123 T3 |
| Markov / time-homogeneous transition reuse (T1) | 44.19 | wqty123 T1 |
| Contiguous row slices | −10.3 | data-gotchas |
| Re-uploading identical files | burns quota (happened 6×) | wqty123 pitfalls II.1 |
| Uploading 500 genes to the embryo board or 32,285 to T2 | rejected | wqty123 pitfalls II.2 |

---

## 6. Noise estimates per metric

| Metric | Estimate | Context / source |
|---|---|---|
| de_score, de_direction, severity_slope, scale_log_ratio | **sd = 0 across scorer seeds** (full-population metrics) | vecbench `results/noise_scorer.txt`; tubager §11.2; wqty123 pitfalls §11 |
| de_score (construction) | sd **1.03 skill** across 5 builder seeds of a resampling method | vecbench `noise_construction.txt` (embryo E6.75→E7.25 holdout, blend baseline) |
| mmd_u | scorer seed sd **1.22 skill** (range 2.88); row-order sd **0.85–1.6** (raw 0.0005–0.001); a "champion file, 6 scorer seeds" run gave **±5.75** (board and statistic not stated, likely T3; **[U]**) | vecbench; tubager §11.2 (heart, pseudo-truth E8.75); wqty123 pitfalls §11 |
| variogram | scorer sd 0.94; row-order sd 1.1–1.6; ±2.20 (wqty123) | same |
| neighborhood_mmd | row-order sd 1.1–1.3 (heart); 0 in the vecbench run (n ≤ 2000 there, so no subsampling) | tubager §11.2; vecbench |
| d2_shape | scorer sd 0.17 (embryo); row-order sd 1.2–2.2 (heart); index-pair sampling sd ≈ 0.0009 raw ≈ 1.7 skill ≈ 0.43 headline | vecbench; tubager §14.2 |
| occupancy_dice | heart sd 0 (row-order independent); **embryo dice sd 0.0075** (near-degenerate principal axes flip the SVD; worth 0.4–0.65 headline) | tubager §11.3 |
| T2 headline | **±0.6 (1 sd) per submission** from row order alone | tubager §11.2 |
| Local pseudo-split headline | copy_last on a synthetic heart pair: 53.11 ± 3.54 across 5 target-subsample seeds | xxx12e README (shows how much the 10% *target* subsample moves scores) |
| Local candidate headline | 3-seed sd 0.04–0.30 (4000-cell subsamples) | i-habib `val-to-test/results/analysis.json` |
| Unweighted composite resolution | do not act on differences below ~0.39–0.57 points | vecbench |

**Rules of thumb used by the teams:**
- treat a difference smaller than 2× the band as a tie (wqty123);
- only de/dir/scale (and heart dice) are noise-free, so a single-ticket conclusion on mmd, variogram, NFS or d2 needs more than 2 skill points (tubager §11.2).

---

## 7. How teams estimate validation-to-test transfer

**Pseudo-pair study** (i-habib [`val-to-test`](https://github.com/i-habib/virtual-embryo-community/blob/main/val-to-test/README.md), MIT).
- Design: heart MERFISH, val pair E8.0 + E8.25 → E8.75, test pair E8.25 + E8.75 → E9.5; 16 candidates; veckit cc3859a; 3 seeds.
- Rank agreement: **Spearman 0.879** (CI 0.81–0.92), Kendall 0.717.
- Two-slot regret: top-2 by validation ≈ **0.00**; "top-1 plus best of another family" 1.42.
- Per-metric Spearman:

| de | dir | mmd | variogram | d2 | dice | scale | NFS |
|---|---|---|---|---|---|---|---|
| 0.85 | 0.58 | 0.68 | 0.89 | **−0.54** | 0.79 | **−0.21** | 0.80 |

- The winning candidates were composition resampling (s = 0.5 / 1.0): val 53.54, test 54.90. Mean-drift candidates were worse than copy_last on both pairs (val 47.1 → 43.1 and test 45.6 → 38.8 as a goes 0.5 → 1.5).
- `configs/oct20_template.json` lists the T1, heart and KO pairs to rerun on 20 October.

**History rescore** (i-habib `history-rescore`). After 20 October, rescore every past submission against the released validation answers using a board protocol (10% subsample, split half, multi-seed) and a full protocol. Then measure local-vs-board difference and Spearman per board, and which local metric tracked the board.

**Honest truth windows** (tubager).
- W2: E8.25 + E9.5 → E8.75. W3: E8.75 copy → E9.5. Embryo: E6.75 + E8.0 → E7.25.
- Lessons:
  1. local d2 is anti-correlated with the board (slope −1.30);
  2. local dice transfers at 0.4–1.4;
  3. **generative or transport X families that win big on hold-out windows lose on the board.** Interpolated cells resemble the endpoint, not the midpoint (tickets 206, 209, 238);
  4. a proxy leaked oracle composition (`--target 8.75` read π(8.75)), which explains several proxy-passed / board-failed tickets (§10.3).

**wqty123** (pitfalls §10, §13).
- "Local proxies can only exclude bad candidates."
- A proxy reading must fall inside its calibration interval: an extrapolated reading was off by +11.24.
- A 500-gene out-of-fold check rejected kp7, which then gained +2.3 on the 32,285-gene board.

**vecbench** `map`.
- de, dir and severity are "truth-set" metrics; a different target has a different DE set, so local movement there is not evidence.
- scale_log_ratio and count_log_ratio are target scalars and cannot be ranked locally.
- Only "direct" distribution and shape metrics transfer.
- `fidelity` computes Kendall tau of local vs board on your own past files; `stability` checks orderability across target redraws.

**xxx12e** `vec_local_score --pred A --pred B`. A paired, same-seed difference B − A; only a min..max band that excludes 0 counts as a result.

---

## 8. External datasets in community use

| Dataset | Stages / content | Who / evidence | Legality (rules §10) | Batch handling noted |
|---|---|---|---|---|
| **GSE193346** (Feng et al. 2022, Nat Commun; heart scRNA) | E9.5 (887 cells in i-habib's normalise example), E14.5 | T1 labels: Druggable 71.65, APC_TDC 65.94, LongWz 61.09, BTW 58.61 **[L]**; heart-extrap: APC_TDC 56.90 ("extsame_w0.5mix"), LongWz 59.51 ("same-stage-assay-alignment") **[L]**; i-habib panel-normalise | E9.5 allowed; **E14.5 allowed** (after E13.5) with disclosure; nothing in (E9.5, E13.5] | normalise over panel genes first. Normalising over all genes then subsetting makes values median 2.3 too low (panel genes hold a median 5.3% of counts) [V]. Same-stage E9.5 assay alignment **[L]** |
| **MOSTA E14.5** (Stereo-seq) | whole embryo E14.5 | LongWz label **[L]** | allowed after E13.5 with disclosure; **MOSTA E10.5–E13.5 sections forbidden** for T1/T2-heart | — |
| **GSE282547** (developing heart Visium / Slide-seq series) | validated sample GSM8645799: E14.5 CD1 WT, 456 spots × 32,285 genes, **500/500 panel genes**, real 2-D coordinates | i-habib external-data-catalog (validated); wqty123 listed it | E14.5+ samples allowed for T1/T2-heart with disclosure; filter out any (E9.5, E13.5] sample; article is CC BY-NC-ND (data licence separate) | 2-D only; no fake z |
| GSE247450 (endothelial MERFISH) | E9.5 WT (315 genes, **68/500** overlap), E15.5 WT | i-habib (validated, narrow) | E9.5 and E15.5 allowed | x/y only |
| GSE150817 | E16.5 heart scRNA, 4,200 cells, 468/500 panel genes | i-habib panel-normalise example | allowed (after E13.5) | panel-first CP10k |
| Extended Mouse Atlas (Imaz-Rosshandler et al.) | E6.5–E9.5, 430,339 cells, 25 GB | divyanshv2001 (T1 ≈ 52); wqty123 audited (leak-test Pearson 0.39–0.50) | T1 allowed; **T2 embryo: remove (E7.25, E8.0)**; **T2 heart: remove (E8.25, E8.75)** | — |
| sc3D / GSE197353 (Slide-seq 3D) | E8.5, E9.0, partial E9.5 | wqty123 closed it (contains Tbx6 mutants; 38% panel overlap) | E8.5 is inside the heart-interp window; WT only for T3 | — |
| GSE208162 | E12.5 ventricle, Gata4 KO/Het | wqty123 used it as a T3 prior (no gain) | **E12.5 is a held-out T1/T2-heart stage, so forbidden for those boards**; same gene as the Gata4 board; irrelevant to β-catenin | — |
| GSE126128 (de Soysa et al. 2019 cardiac atlas; identity not verified **[U]**) | early heart time course | T1 labels bingo 55.16 / oxtium-ai 55.02 **[L]** | check stages vs (E8.25, E8.75) for heart-interp | — |
| GSE230531 | "E8.5 + E14.5" per label | Leslie Lab T1 52.71 **[L]** | E14.5 allowed | — |
| GSE278603 (whole embryo) | — | wqty123: only 15,669 shared genes; stage direction cosine −0.383 with the global direction | stage-check | — |
| Cross-species: rabbit PRJNA943458, pig GSE236766, macaque GSE193007 / CNP0007017, human CS7 E-MTAB-9388, human heart GSE283967 / EGAS…, chicken GSE149457 | various | EmbryoMatch registry (no measured gains) | **rules give no cross-species stage-equivalence policy; EmbryoMatch recommends asking the organisers** | orthology plus unmeasured-gene mask |
| Mouse GO GAF | annotations | i-habib | non-measured prior; CC BY 4.0 | — |

**How teams handle batch effects.**
- wqty123 deliberately avoids Harmony/scVI, because integration would treat true temporal change as batch. They use only gene-wise rank or program-shape alignment.
- tubager measured that a dataset-level per-gene offset dominates the cross-stage pseudobulk error even between the challenge's own MERFISH stages (§11.4).
- Label evidence suggests the stronger external-data teams calibrate on a same-stage pair (external E9.5 vs released E9.5) before using the later stage **[L]**.

**Pretrained models** (wqty123 `reference-links.en.md` §III).
- TranscriptFormer (non-Sapiens) and UCE are flagged high-risk: their corpora include mouse embryo time courses that could contain held-out stages.
- TF-Sapiens, scGPT (human) and Geneformer are judged usable.

---

## 9. Legality red flags seen in public repos (do not copy these practices)

1. **Solving for target properties from scores or anchors.** tubager derived the true target RMS radii from board read-backs and the published floor anchors:
   - embryo [redacted 2026-10-19: A] from two embryo board readings ([redacted 2026-10-19: A]);
   - heart E8.5 [redacted 2026-10-19: A] from the floor anchor ([redacted 2026-10-19: A]), refined to [redacted 2026-10-19: A] from a board reading;
   - heart E10.5 [redacted 2026-10-19: A], the floor anchor inverted.

   They then set their predictions to exactly those radii (`t2_p2_val_2026-09-05.md` §2.3, §14.1; config comment "官方 copy_last TSR ≈ [redacted 2026-10-19: A], β = −TSR"). Rules §10 (brief): "Score returns may be used to choose among predictions, not to compute them (no inverting the skill map / differencing scores to solve for target properties)". vecbench's README states the same boundary. PR #6 shows the organisers see the scale exploit ("rescaling the floor's coordinates took the old size score from 50 to 100"). scale_log_ratio = 100 on 85, 109 and 41 of the 200 listed rows on the three T2 boards (landscape §0), which suggests the practice is widespread. **For the test, choose the radius from a model of released stages, and document how.**
2. **The organisers' own code publishes held-out radii.** The veckit **0.1.2** docstring (PyPI wheel / PR #6, `common/shape_metrics.py`) lists heart RMS radius and nearest-neighbour distance by stage, including held-out stages and the E12.5 test target. **[Values deliberately redacted from this digest: they are measurements of held-out stages and must not inform any prediction.]** Do not open that docstring while building size models; [redacted 2026-10-19: A]
3. **The organisers' public atlas and knowledge graph** (`kg.virtualembryo.ai`, documented in aristoteleo/PantheonOS `skills/virtualembryo/SKILL.md`). It lists 3D spatial datasets including `digiembryo_e7_5_rep1`, MOSTA sections and Spateo E9.5/E11.5 reconstructions. Anything in (E7.25, E8.0) is forbidden for T2 embryo, and anything in (E9.5, E13.5] for T1/T2-heart. Do not pull these for those boards.
4. **"Authorized cross-team reuse"** label (REVISE, embryo 69.31 **[L]**) and several identical metric vectors across teams (landscape §1.9). Code sharing across teams is prohibited; avoid any appearance of it.
5. **GSE208162** is inside the T1/T2-heart excluded window (E12.5).
6. **Cross-species data** analogous to protected mouse windows: unresolved; ask.

---

## 10. Rule and scoring changes recorded publicly

**From the xxx12e rules watch** ([`rules-watch/CHANGES.md`](https://github.com/xxx12e/vec-community-kit/blob/main/rules-watch/CHANGES.md); only 3 runs recorded a change):
- **2026-09-30 (baseline):**
  - 5 boards; T1:val 32,285 genes, cells 1000–5118; embryo 498 genes, [redacted 2026-10-19: A]–5000; heart-extrap 500 genes, 1000–[redacted 2026-10-19: A]; heart-interp 1000–17616; **T3:gata4 1000–7449**;
  - phase p2, daily quota 8, nominations closed;
  - veckit `main` 46d41e6 (0.1.1); branches main and feat/metric-hardening-0.2.0.
- **2026-10-03:** veckit `main` moved 46d41e6 → **cc3859a** (PR #9, right-handed PCA frame).
- **2026-10-16:**
  - FAQ and Data pages: the outright extrapolation window now explicitly covers "(Task 1, and Task 2 in the heart setting)";
  - Evaluation submissions pages for T1/T2/T3: "only the 1,000-cell minimum is enforced" became "only the board's own minimum";
  - new veckit branch `docs/library-size-ratio-no-renormalisation` (f56a6e3; PR #12: the scorer never renormalises, so always submit log1p CP10k).

**Other dated notes in the kit README:**
- Evaluation pages were corrected after the organisers reviewed the kit (2026-09-22): `max_cells` is authoritative.
- 2026-09-30: the organisers said the subsample seed "will depend on each submission".

**EmbryoMatch rules notes.** The 2026-09-02 transcription of the T2 interpolation windows was wrong and was corrected on 2026-09-10 to the strict open intervals (E7.25, E8.0) and (E8.25, E8.75). This matches the rules §10 restatement of 26 Aug.

**Scorer versions** (i-habib `scorer-versions`).
- PyPI 0.1.2 equals the unmerged PR #6. It reports severity_ratio and size_fidelity, which the board does not rank, and it lacks the rotation fix.
- `main` cc3859a is 0.1.1 plus the rotation fix only.
- T1 metrics are identical across all variants. Install with `pip install "git+https://github.com/aristoteleo/veckit@cc3859ab2a2305efc4544815d3f44f660efa77d4"`.

**Final-phase questions the official text does not answer** (i-habib `test-phase-switchover`):
- how the two files per board are designated;
- whether P3 returns per-metric breakdowns;
- whether the §10 held-out list still applies to released validation stages;
- the number of uploads per day.

---

## 11. Implications for our P3 plan (inference from the above) [I]

- **T2 embryo E7.75 (test).** At P3, E7.5 truth becomes training input, so the brackets become E7.5 and E8.0, with clock w = (7.75 − 7.5)/(8.0 − 7.5) = 0.5. If E7.5 is not used, the brackets are E7.25 and E8.0 with w = 2/3. tubager's "geometry from the nearer-left anchor" rule was learned at w ≤ 0.5; at w = 2/3 the nearest anchor is E8.0, so validate the anchor choice on the released E7.5 answers.
- **T2 heart E12.5.** E10.5 truth becomes the last observed stage. The recipe-2 chassis (real E10.5 cells, frozen composition, keep-zeros per-cluster Δ from E9.5→E10.5 with α < 1, carrier selection) is the strongest public template, but it was tuned on a 1-day step and E12.5 is 2 days away. Geometry is the open problem: top teams reach d2 ≈ 100 on extrapolation by unknown means. The scale must come from a defensible model, [redacted 2026-10-19: A].
- **T1 E12.5.** Combine:
  - real-cell carriers from E10.5 (released at P3);
  - composition and emergence handling (harmonised labels);
  - keep-zeros per-type deltas;
  - variance and maturity reweighting;
  - legally bracketing E14.5 data (Feng 2022 / MOSTA E14.5) with same-stage assay alignment.
- **T3 β-catenin.** Gata4 (2 replicates, at the test stage) plus Mab21l2 become training data. The mix + flatten + dilute framework is the only measured T3 recipe. Test with Gata4 as pseudo-test before using it on β-catenin. Keep severity β near 1.
- **Selection.** Run history-rescore and the val-to-test pseudo-pairs on 20 October, pick the top 2 by validation per board, and ignore local d2/scale rankings.

---

## 12. Things I could not verify
- The exact definition of the wqty123 "maturity" score (kp7) and the identity of the "real KO cells" in its T3 mix families.
- The exact local-field construction for tubager's Hungarian re-pairing (tickets 192–200); the script is not in the repo.
- tubager's claim of "T1 #1" (`t2_p2_history.md` line 27, "T1 已是第一"); the repo contains no T1 board scores.
- wqty123's "mmd ±5.75 over 6 scorer seeds": board and statistic not stated.
- Whether "Ivan & Arnav" (Human #2) is i-habib's team (name match only).
- The identity of GSE126128 and GSE230531 and whether their stages are legal for each board.
- Methods of the top-5 Human teams: no public write-ups exist.

---

## References (URLs)

**Methods journals and code:**
- wqty123 methods journal: https://github.com/wqty123/virtual-embryo-challenge-methods (README.en.md, experiments_t1.en.md, experiments_t2_embryo_interp.en.md, experiments_t2_heart_interp.en.md, experiments_t2_heart_extrap.en.md, experiments_t3.en.md, pitfalls.en.md, tools.en.md, reference-links.en.md)
- tubager T2: https://github.com/tubager/vec-t2 (doc/t2_p2_history.md, doc/t2_p2_val_2026-09-05.md, doc/t2_technical_plan_v2.md, doc/t2_p2_plan_first_2026-09-13.md, doc/t2_p2_weekly_2026-09-16_22.md, doc/t2_flow_eval.md, configs/t2.yaml, scripts/75_cluster_quad_curv.py, scripts/78_extrap_g0_clust_delta.py, scripts/30_cluster_moment_morph.py, src/t2/ot.py)

**i-habib community (all MIT):**
- repo root: https://github.com/i-habib/virtual-embryo-community
- scorer-versions: https://github.com/i-habib/virtual-embryo-community/tree/main/scorer-versions
- val-to-test: https://github.com/i-habib/virtual-embryo-community/tree/main/val-to-test
- history-rescore: https://github.com/i-habib/virtual-embryo-community/tree/main/history-rescore
- count-recovery: https://github.com/i-habib/virtual-embryo-community/tree/main/count-recovery
- panel-normalise: https://github.com/i-habib/virtual-embryo-community/tree/main/panel-normalise
- external-data-catalog: https://github.com/i-habib/virtual-embryo-community/tree/main/external-data-catalog (standalone copy: https://github.com/i-habib/external-data-catalog)
- submission-blender: https://github.com/i-habib/virtual-embryo-community/tree/main/submission-blender
- ensemble-lab: https://github.com/i-habib/virtual-embryo-community/tree/main/ensemble-lab
- intuition-lab: https://github.com/i-habib/virtual-embryo-community/tree/main/intuition-lab
- metric-movies: https://github.com/i-habib/virtual-embryo-community/tree/main/metric-movies
- morph-to-truth: https://github.com/i-habib/virtual-embryo-community/tree/main/morph-to-truth
- leaderboard-anatomy: https://github.com/i-habib/virtual-embryo-community/tree/main/leaderboard-anatomy
- test-phase-switchover: https://github.com/i-habib/virtual-embryo-community/tree/main/test-phase-switchover
- data-safety: https://github.com/i-habib/virtual-embryo-community/tree/main/data-safety

**Other community tools:**
- xxx12e kit and rules watch: https://github.com/xxx12e/vec-community-kit (rules-watch/CHANGES.md, docs/tutorial_en.md)
- vecbench: https://github.com/AJK0921/vecbench (vecbench/README.md, vecbench/results/*.txt)
- Leaderboard insights: https://github.com/WaelYTz/vec-leaderboard-insights
- EmbryoMatch: https://github.com/ShaoliZhao/virtual-embryo-EmbryoMatch
- Data gotchas: https://github.com/Shashwat-srivastav/virtualembryo-data-gotchas
- Interpolation baseline: https://github.com/Shashwat-srivastav/VirtualEmbryo-basline
- Cell-type harmonisation: https://github.com/Anky9972/celltype-harmonization (analysis/REPORT.md)
- Agent T1 attempts: https://github.com/divyanshv2001/virtualembryo2026 (outputs/t1_submissions/*/official_result.json)
- Scale check: https://github.com/wqty123/vec-scalecheck
- Submit check: https://github.com/gh-dv-openclaw/vec-submit-check
- random-guy-05 tools: https://github.com/random-guy-05/vec-samplecheck , https://github.com/random-guy-05/vec-replicate-range , https://github.com/random-guy-05/vec-failure-gallery
- cadentann: https://github.com/cadentann/vec-scorelens , https://github.com/cadentann/virtual-embryo-community-kit

**Organisers:**
- PantheonOS virtual embryo skill: https://github.com/aristoteleo/PantheonOS/blob/main/pantheon/factory/templates/skills/virtualembryo/SKILL.md
- veckit: https://github.com/aristoteleo/veckit (main cc3859a); PyPI 0.1.2: https://pypi.org/project/veckit/0.1.2/ (docstrings in common/shape_metrics.py, common/core_metrics.py)

**Datasets:**
- GSE193346: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193346
- GSE282547: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282547
- GSE247450: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247450
- GSE150817: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150817
- GSE197353: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197353
- GSE208162: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208162
- Extended Mouse Atlas: https://marionilab.github.io/ExtendedMouseAtlas/
- Mouse GO GAF: https://current.geneontology.org/annotations/gaf/MOUSE-mod.gaf.gz

**Challenge pages:**
- Rules: https://virtualembryo.ai/challenge/rules
- Data: https://virtualembryo.ai/challenge/data
- Board index: https://virtualembryo.ai/challenge/panels/index.json

**Local copies (this run):**
- Downloaded repo files: `SCRATCH/community/repos/`
- Trees and metadata: `.../scratchpad/community/trees/`
- Searches and community register list: `.../scratchpad/community/search/`
