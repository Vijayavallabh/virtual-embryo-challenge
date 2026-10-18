# Task 2 methods: spatiotemporal interpolation and extrapolation of 3D MERFISH

Virtual Embryo Challenge (NeurIPS 2026), Human Team track. Written 2026-10-18. The test phase (P3) opens 2026-10-20; final submissions are due 2026-12-02.

Scope: the two Task 2 test boards. These are **embryo E7.75** (interpolation between E7.5 and E8.0) and **heart E12.5** (extrapolation from E8.25–E10.5).

Confidence labels:
- **[High]**: verified from the scorer code, the challenge site, the live leaderboard API, data files or a primary paper.
- **[Medium]**: reasoned from verified facts, or from a secondary source.
- **[Low]**: judgement or extrapolation.
- **[Unverified]**: I could not confirm it.

Scratch artefacts used here:
- `SCRATCH/scripts/inspect_samples.py`: inspects the sample h5ads.
- `SCRATCH/scripts/nbhd_density_sim.py`: synthetic test of how cell count affects neighbourhood MMD.
- `SCRATCH/dl_js/resp_board_T2_*.json`: live leaderboard snapshots, pulled 2026-10-18 at about 16:14 UTC.

---

## 0. Executive summary

1. **What a T2 board score is [High].** Per board:
   - S = 25·[E + D + g·M + N], where
   - E = (s_de + s_dir)/2,
   - D = 0.6·s_mmd + 0.4·s_vario,
   - M = (s_d2 + s_dice + s_scale)/3,
   - N = s_nbhd.

   The structure gate is **g = min(1, N/0.5)**. This reproduces all 208 published T2 headlines to within 0.000 (§1.2). A plain product g = N is rejected. So the gate only bites when the neighbourhood term is worse than copy_last. Above that, `neighborhood_mmd` is simply a 25% metric.

2. **The live board is public.** It shows what P2 methods reach on the validation boards (§1.4).

   | Board | Top headline | Top-20 median |
   |---|---|---|
   | Embryo E7.5 interpolation | 77.7 | 75.7 |
   | Heart E8.5 interpolation | 75.4 | 71.0 |
   | Heart E10.5 extrapolation | 67.8 | 64.2 |

   Almost every top team scores 99–100 on `scale_log_ratio`. This is most likely not a modelling feat: the public `index.json` floor anchor gives the target RMS directly (e.g. [redacted 2026-10-19: A]), and score-based bisection would also work. [Medium] On the test boards the size must actually be predicted, unless test anchors are published (see the rule risks). So test scores should fall, and **genuine size, shape and neighbourhood prediction is where to gain an edge.**

3. **Coordinates are in micrometres [High for "µm-scale physical"; Medium for the exact section pitch].** In the 150-cell samples:
   - z takes discrete values exactly 15 units apart (section pitch);
   - xy extents are 0.6–1.1 mm;
   - RMS radius is 309.6 ± 7.7 at E9.25 and 333.5 ± 9.5 at E9.5.

   The organisers' own table (PR #6, released stages) gives:

   | Stage | RMS radius (µm) | Median NN distance (µm) |
   |---|---|---|
   | E8.25 | [redacted 2026-10-19: A] | [redacted 2026-10-19: A] |
   | E8.75 | [redacted 2026-10-19: A] | [redacted 2026-10-19: A] |
   | E9.5 | [redacted 2026-10-19: A] | [redacted 2026-10-19: A] |

   The RMS radius is **not monotone in time**. It tracks the imaged field of view, not heart size alone. [redacted 2026-10-19: A] A single log-RMS growth curve across the whole heart ladder is therefore invalid (§2).

4. **Integrity issue [High].** The open, unmerged PR #6 in `aristoteleo/veckit` (the organisers' branch, 2026-08-25) contains a docstring that tabulates RMS radius and NN distance for **held-out stages, including the E12.5 test target**. Treat those entries as held-out data:
   - do not use them;
   - do not copy them into notebooks;
   - notify the organisers.

   This report does not reproduce them (§1.5). PR #6 would also replace `scale_log_ratio` with `size_fidelity`. That metric scores the tissue radius **and** the median NN spacing, which would make cell count and spacing part of the ranked size question.

5. **Cell count matters even though it is "not scored" [High for the mechanism, Medium for the size of the effect].** `neighborhood_mmd` averages expression over the 15 nearest neighbours among the *submitted* cells. The truth half is only 5% of the target's cells, so its 15-NN neighbourhoods are spatially large. In a synthetic test that mirrors the scorer:

   | Submitted count / truth count | Neighbourhood MMD |
   |---|---|
   | 1× | 0.0013 |
   | 2× | 0.0016 |
   | 4× | 0.0074 |
   | 0.5× | 0.0057 |

   A uniform rescale of the coordinates changed nothing. **Submit about as many cells as the truth half** (`truth_cells` in `index.json`). This also protects against a `size_fidelity` switch.

6. **Recommended pipelines.** Both are built on "transport real cells", not on synthesising expression, because the MERFISH matrix is very sparse: the median cell has 69 counts and 35 detected genes at E9.5.
   - **Embryo E7.75:**
     1. Fused Gromov-Wasserstein (FGW) coupling of canonicalised E7.5 and E8.0.
     2. Proper-rotation Procrustes alignment.
     3. McCann displacement interpolation of positions at t ≈ 0.5.
     4. Expression as a mixture of coupled partners, using count-space thinning rather than additive log shifts.
     5. Log-RMS interpolated with a monotone curve through E6.75, E7.25, E7.5 and E8.0.
     6. Cell count matched to the truth.
   - **Heart E12.5:**
     1. Start from real E10.5 tissue.
     2. Apply a damped, scope-aware composition and per-type expression extrapolation. Its direction can be anchored on *allowed* post-E13.5 external data, after calibrating scRNA to MERFISH on the competition's paired E8.5/E9.5/E10.5 stages.
     3. Insert emergent populations (epicardium-derived cells, coronary endothelium, valve mesenchyme, smooth muscle) into their anatomical niches.
     4. Grow the coordinates with a size model fixed in advance from training stages.
     5. Spend the two official test submissions on **two scope/size hypotheses**; the better one counts.

---

## 1. Scorer facts that change the method choice

### 1.1 Metric mechanics, per T2 board [High; from veckit code at commit cc3859a]

- **Reference stage.** `ref` is the preceding stage. On the validation boards, `ref_cells` = 10% of E7.25 (1,330) for embryo, E8.25 (5,872) for heart interpolation and E9.5 (5,374) for heart extrapolation.
  - **For the test boards, check `ref_cells` in the test `index.json`.** It shows whether the embryo reference is E7.5 ([redacted 2026-10-19: A]) or E7.25 (1,330). The heart reference should be E10.5 ([redacted 2026-10-19: A]). [Medium]
- **Subsampling.** Every stage is subsampled to 10%, and the target is then halved.
  - The truth halves are [redacted 2026-10-19: A] cells (embryo E7.5), [redacted 2026-10-19: A] (heart E8.5) and [redacted 2026-10-19: A] (heart E10.5).
  - Implied full stage sizes are about [redacted 2026-10-19: A] (E7.5), 13.3k (E7.25), [redacted 2026-10-19: A] (E8.5), 58.7k (E8.25), 53.7k (E9.5) and [redacted 2026-10-19: A] (E10.5).
  - These counts reflect MERFISH sectioning and sampling, not biology: [redacted 2026-10-19: A]. **Do not use them as a growth signal.** [High]
- **DE metrics.** `de_score` and `de_direction` are rank-based on dp = pb(pred) − pb(ref). Any positive rescaling of dp is a no-op, so only the **direction** of the predicted pseudobulk change matters. `de_direction` is a rank partial correlation over **all** 498 or 500 genes, controlling for the rank of pb(ref).
- **mmd_u.** Computed in a 30-PC space fitted on the truth, on 2,000 cells, with RBF bandwidths from 0.25× to 4× the truth median. The variogram uses 20,000 gene pairs, 1,500 cells and |x_i − x_j|^0.5. Both are sensitive to sparsity and to the discrete level structure of log1p(1e4·c/Σc) values. For context, a single count in a 69-count cell maps to log1p(145) ≈ 5.
- **Shape terms.**
  - `d2_shape` compares within-cloud distances divided by RMS. It is blind to reflection.
  - `occupancy_dice` canonicalises (centre, PCA, divide by RMS, right-handed frame) and counts occupancy on a 16³ grid over [−3, 3] RMS (voxel 0.375 RMS). It maximises over the four proper sign flips only. **It does not permute axes and does not reflect.** Both clouds are subsampled to min(n_pred, n_true) first.
  - `scale_log_ratio` = log(RMS_pred / RMS_true), scored on its absolute value.
- **neighborhood_mmd.** Each cell's 15-NN mean expression (self included) is computed on the full submitted set. Then `mmd_u` is run on 2,000 of these profiles, with the PCA fitted on the truth's profiles. **The kNN graph is invariant to uniform scaling, so absolute size does not enter this term; cell density does** (§1.3).

### 1.2 The structure gate, re-derived from the live board [High]

I pulled `GET https://api.virtualembryo.ai/challenge/leaderboard?task=T2&board=<key>` for all three boards. These are 200 rows each plus reference rows. I then recomputed each headline from its eight per-metric skills.

- g = min(1, N/0.5) gives **max |error| = 0.000 on 208/208 rows**, including the only row with N < 50: the heart-extrap `pseudobulk_shift` reference, 47.76.
- No gate is off by 0.27 on that row. g = N (a plain product) is off by up to 10 points.
- A sibling report (`landscape.md`) independently confirms this on 600 rows across two snapshots.

**Implication.** Copy_last scores 50 because its gate is 1. A model whose neighbourhoods beat copy_last gets full shape credit, and every gain in `neighborhood_mmd` above the floor is a plain 25%-weight gain. A model whose positions are structureless loses shape credit in proportion; that is how `ctrl_random_cube` was killed.

### 1.3 Cell count and density: a scored lever hidden behind "count is not scored" [High mechanism, Medium magnitude]

**Why.** The truth half holds about 5% of the stage's cells, so its 15-NN neighbourhoods span many tissue layers. If a submission has k times the truth's density, its neighbourhoods are k^(1/3) times smaller in radius. Their averaged profiles are then less mixed, and the neighbourhood-profile distribution shifts.

**Test.** I re-implemented the scorer's kNN-mean and unbiased-MMD logic in my own script; no downloaded code was run. The synthetic tissue had 200k cells, 8 spatially patterned types at scales of 80–400 µm and 100 genes. The truth was 4,000 cells. Result: NN-MMD at 1×/2×/4×/8× the truth count = 0.0013/0.0016/0.0074/0.0138, and at 0.5×/0.25× = 0.0057/0.0169. Coordinates ×1.7 gave an identical value (0.00110). Random positions in the bounding box gave 0.0893.

**Size of the effect on real boards.** Converted with the real anchors, an extra 0.006 of NN-MMD costs:
- about 0.03 skill on embryo (floor 0.211, ceiling 0.011);
- about 0.08 on heart E10.5 (floor 0.072, ceiling 0.0004).

So expect **roughly 1–2 board points** from count mismatches of 4× or more. Real tissue mixes cell types at finer scales than this simulation, so the effect could be larger.

**Rule.**
- Submit n_pred ≈ `truth_cells` ([redacted 2026-10-19: A]).
- If the test index does not publish it, aim within 2× of a best estimate.
- Never go below about 4,000 on heart E12.5. Under the proposed `size_fidelity`, NN distance is measured on min(n, 4000) cells per cloud.
- Avoid exact duplicate coordinates, which cause degenerate kNN; veckit issue #11 covers this for Moran's I.

### 1.4 What P2 teams achieve per metric [High; snapshot 2026-10-18 16:14 UTC]

Medians of the top-20 rows, skill × 100:

| Board | Headline | de_score | de_dir | mmd_u | variogram | d2 | dice | scale | nbhd |
|---|---|---|---|---|---|---|---|---|---|
| Embryo E7.5 interp | 75.7 (top 77.7) | 58.2 | 71.1 | 77.2 | 65.9 | 80.6 | 100.0 | 99.4 | 74.3 |
| Heart E8.5 interp | 71.0 (top 75.4) | 63.0 | 69.2 | 69.4 | 58.7 | 92.4 | 56.5 | 100.0 | 71.7 |
| Heart E10.5 extrap | 64.2 (top 67.8) | 59.1 | 61.2 | 63.4 | 59.0 | 74.3 | 53.0 | 100.0 | 64.5 |

Reading:
1. **Scale is saturated on validation** because the anchors leak the target RMS. Expect a drop on the test.
2. **Embryo dice saturates at 100.** With only [redacted 2026-10-19: A] truth cells, the half-vs-half ceiling dice (0.766) is limited by sampling. A prediction that spreads its points evenly over the right canonical region occupies more of the truth's voxels than a random half does, so it beats the ceiling (clipped at 100). A mild uniformisation of point placement (Lloyd or blue-noise relaxation within the predicted shape) is a cheap test on E7.75, **if** local validation shows it does not hurt d2 or nbhd. [Medium]
3. **Heart dice is hard** (about 50–60 even for leaders). The canonical frame has to match: PCA axis order and handedness (§3.4).
4. **Human-track leaders** to calibrate against:
   - Embryo: Omazili Manbo (77.7), GISL (76.1), Rapunzel (76.0).
   - Heart E8.5: Omazili Manbo (75.4), Ivan & Arnav (73.7), EmbryoForge (73.7).
   - Heart E10.5: Omazili Manbo (67.3), John Roos (67.0), Ivan & Arnav (66.5).
5. `pseudobulk_shift` (constant velocity) scores 55.2 on embryo and **47.8 (below the floor)** on heart extrapolation. This is consistent with the organisers' note that constant-velocity extrapolation loses to copy_last.

### 1.5 PR #6: a possible metric change and a leak [High]

`aristoteleo/veckit` PR #6 ("feat(metrics)!: PSS/TSR scoring change, sync with vec_baselines, bump to 0.2.0") is open and unmerged. It was opened on 2026-08-25 from the branch `aristoteleo:feat/metric-hardening-0.2.0`.

**The proposed change.** It replaces the ranked `scale_log_ratio` with
`size_fidelity = ½·[min(r, 1/r) + min(q, 1/q)]`, where
- r = RMS ratio;
- q = ratio of median NN distances, each computed on min(n, 4000) cells.

The live `index.json` anchors still use `scale_log_ratio` and `severity_slope`, so the board has **not** adopted it yet. However:
- PyPI's latest `veckit` is **0.1.2** (2026-08-25, the same day as the PR). A sibling report says it ships the PR #6 metrics. **Score locally from the GitHub `main` checkout**, which keeps `scale_log_ratio` and the #7 right-handed fix, not from `pip install veckit`.
- If `size_fidelity` is adopted, the remedy is the same as in §1.3: real-cell spacing (about [redacted 2026-10-19: A] in the organisers' units) plus a matched count.

**The leak.** The `size_fidelity` docstring contains a per-stage table of RMS and NN distance for E8.25, E8.5, E8.75, E9.5, E10.5 **and E12.5**.
- The E8.5 and E10.5 values match what the public anchors imply, and those stages are released on 2026-10-20.
- **The E12.5 entries describe the hidden test target.**
- Under Official Rules §10 ("measured data from a held-out stage … by any route"), and the clause against solving for target properties from scores, using them would be a violation.

Actions:
- (a) Do not use or copy them.
- (b) Email the organisers about the leak.
- (c) Log that our E12.5 size model was fixed from training stages [redacted 2026-10-19: A] only.

Disclosure: I saw the value while reading the PR. The size recommendations in §5.3 were derived from training-stage data and anchor-derived ratios. Still, use the pre-registered procedure there, not any single number, and hedge across the two submissions.

### 1.6 Two more scorer properties worth exploiting legitimately [Medium]

- **Handedness.** `occupancy_dice` will not credit a mirror image. All coordinate transforms must be proper (det = +1): Procrustes without reflection, and PCA bases with det fixed to +1.
  - Also check whether released stages share one handedness convention. Run FGW-coupled Procrustes with and without reflection between consecutive stages.
  - If some pairs fit much better mirrored, the acquisition frames differ. The E12.5 frame is then a coin flip, which is a reason to hedge across the two submissions.
- **The PCA axis order defines the canonical frame.** If two principal SDs are close in the truth, small shape errors swap axes and dice collapses. When predicting shape, match the predicted principal-SD **ratios** as well as the RMS.

---

## 2. Coordinate units and what absolute RMS growth implies (Q4)

### 2.1 Inspection of the veckit samples [High]

Files: `ext/veckit/data/sample_heart_9.25.h5ad` and `sample_heart_9.5.h5ad`. Both have 150 cells × 500 genes. Read with anndata via `uv run --no-project`.

| | sample_heart_9.25 | sample_heart_9.5 |
|---|---|---|
| obs | slide (`M118`), section (40 unique), `z_coord`, celltype (14: Paraxial Mesoderm, Foregut, vCM, Branchial Arch, Endothelium, Endocardium, LPM, DPM, OFT/RV-CM, SV-CM, Surface Ectoderm, Neural Tube, Pericardium, Proepicardium), Level1/Level3 | celltype (20: V-CM, A-CM, Peri, Lateral-PAM, BEC, Chamber-Endo, aPHM, SV-CM, …), cm_celltype |
| obsm | X_pca, X_pca_harmony, X_umap, spatial (slide coords, ±6–7 mm), spatial_2D, **spatial_3D**, spatial_nonrigid_aligned_3D | spatial_2D (slide coords, x/y 0–6 mm, z ≈ 0), **spatial_3D** |
| X | not log-normalised (max 63; `counts` layer is integer) | log1p(1e4·c/Σc): every cell has expm1(X).sum = 10,000; `counts` layer median 69 counts/cell, 35 genes/cell |
| spatial_3D extents (x, y, z) | 586 × 912 × 645 | 632 × 1,065 × 585 |
| z values | 40 unique, **step exactly 15.0** | 39 unique, **step exactly 15.0** |
| RMS radius (bootstrap SE) | 309.6 (7.7) | 333.5 (9.5) |
| principal SDs | 216 / 170 / 142 | 247 / 170 / 146 |
| median NN distance (150 cells) | 60 | 65 |

**Interpretation.** Units are micrometres:
- z is a section index × 15 µm;
- xy extents of 0.6–1.1 mm match a heart-region cross-section at E9.25–E9.5;
- the E9.5 sample's RMS (333.5 ± 9.5) [redacted 2026-10-19: A], so the sample shares the challenge frame and units.

The same lab's human MERFISH+ heart work used 16 µm cryosections and 15 µm of z imaging per section (Kern et al. 2025 preprint). This supports a pitch of about 15 µm. [Medium]

The E9.25 sample is **not** challenge data: no E9.25 stage is released for T2, and its X is not log-normalised. Use it only as a structural example.

**Scope.** The E9.25 sample includes neural tube, paraxial mesoderm and surface ectoderm. The E9.5 vocabulary has no neural tube or somites, but does include foregut, hepatocytes, septum transversum, pharyngeal mesoderm and NCC. The "heart scope" is therefore a heart-centred region whose extent varies by stage. [Medium]

### 2.2 RMS radius across stages [High for the numbers, Medium for the interpretation]

Heart, released stages (PR #6 table plus our sample): [redacted 2026-10-19: A].

From the anchors: [redacted 2026-10-19: A].

Log rates per day:

| Interval | Log rate per day | Note |
|---|---|---|
| E8.25 → E8.5 | [redacted 2026-10-19: A] | over 0.25 d |
| E8.5 → E8.75 | [redacted 2026-10-19: A] | |
| E8.75 → E9.5 | [redacted 2026-10-19: A] | |
| E9.5 → E10.5 | [redacted 2026-10-19: A] | |

The heart does not shrink. **The field of view changes**: E8.25_late has the largest cell count, 58.7k, and a wide field. [redacted 2026-10-19: A]

Embryo: the anchor gives RMS(E7.5)/RMS(E7.25) = [redacted 2026-10-19: A], i.e. [redacted 2026-10-19: A]. [redacted 2026-10-19: A]
- Isotropic proliferation with a [redacted 2026-10-19: A] generation time (Snow 1977) would give about [redacted 2026-10-19: A] linear growth per 6 h at most.
- [redacted 2026-10-19: A]

[Medium/Low]

**What this means for modelling.**
1. Absolute size matters only for `scale_log_ratio` (and for `size_fidelity` if adopted). `d2_shape` and `occupancy_dice` divide by RMS, and `neighborhood_mmd` is scale-invariant.
2. **Embryo:** a monotone curve for log RMS(t) through the four stages should be accurate. PCHIP or quadratic in t, checked by leave-one-out on E7.25 and E7.5.
3. **Heart:** do **not** fit one curve through E8.25–E10.5. Model a scope-invariant core (cardiac cells only: CM, endocardium, epicardium, cushion) separately from the "shell" (pharyngeal, foregut, liver, pericardial body wall), and decide the E12.5 scope explicitly (§5.3).

Organisers' warning about copy_last on interpolation boards: the site says copy_last is "roughly [redacted 2026-10-19: A] too small by volume" there. This is inconsistent with the current anchors. [redacted 2026-10-19: A] It probably comes from the retired bounding-box volume metric. [Medium]

---

## 3. Methods survey and code status (Q1)

Column legend for the main table:
- **Rel.** = relevance to our boards, H/M/L.
- **Code** = repository, licence and last push, checked 2026-10-18 with `gh api`.
- **PyPI** = latest version and upload date, from pypi.org JSON.

### 3.1 Spatiotemporal OT and alignment

| Method | What it does | Rel. | Code / PyPI | Reference |
|---|---|---|---|---|
| **Spateo** (Qiu et al., Cell 2024; same lab as the organisers; Qingquan Zhang is a co-author) | 3D reconstruction of serial sections; `morpho_align` (partial non-rigid alignment, modes SN-S/SN-N, GPU, KL or Euclidean expression dissimilarity); `tdr.morphofield_sparsevfc`/GP "morphometric vector fields" learned from aligned stage pairs, with curl, divergence, acceleration and so on; mesh and voxel models; interpolation (sparseVFC, GP, deep) | H (alignment and vector fields between our stage pairs; non-rigid correspondences for displacement interpolation) | aristoteleo/spateo-release (BSD-2, 345★, pushed 2026-02-25); PyPI `spateo-release` 1.1.1 (2025-03-25) | [10.1016/j.cell.2024.10.011](https://doi.org/10.1016/j.cell.2024.10.011); erratum [10.1016/j.cell.2025.02.026](https://doi.org/10.1016/j.cell.2025.02.026). Its mouse heart analysis used whole-embryo Stereo-seq atlases at E9.5 and E11.5. **E11.5 is inside the heart extrapolation exclusion window; E9.5 is allowed.** |
| **moscot** `SpatioTemporalProblem` (Klein et al., Nature 2025) | FGW between consecutive time points: GW term on within-stage spatial distances plus a linear term on expression (`joint_attr`). Unbalanced via `tau_a`/`tau_b`, low rank via `rank`, `epsilon` defaults to 1e-3, `alpha` defaults to 0.5. Growth-rate priors through proliferation/apoptosis scores. ott-jax backend on GPU. Used on mouse embryogenesis (MOSTA). | H (ready-made FGW coupling at scale) | theislab/moscot (BSD-3, pushed 2026-10-14); PyPI `moscot` 0.5.2 (2026-08-25, py ≥ 3.11) | [10.1038/s41586-024-08453-2](https://doi.org/10.1038/s41586-024-08453-2) |
| **DeST-OT** (Halmos et al., Cell Systems 2025) | Semi-relaxed OT between developmental slices; models growth, death and differentiation; growth-distortion and migration metrics | M (semi-relaxed coupling when types appear; growth-rate estimates) | raphael-group/DeST_OT (BSD-3, pushed 2025-03-19); no PyPI package found | [10.1016/j.cels.2024.12.001](https://doi.org/10.1016/j.cels.2024.12.001) |
| **STORIES** (Huizing et al., Nat Methods 2025) | Learns a spatially informed potential (Wasserstein gradient flow) with FGW across time points; benchmarked on Stereo-seq atlases | M (a learned potential could forecast; small number of time points is a concern) | cantinilab/stories (**GPL-3.0**, pushed 2025-10-21); PyPI `stories-jax` 0.1.1 (2025-10-21) | [10.1038/s41592-025-02855-4](https://doi.org/10.1038/s41592-025-02855-4) |
| **PASTE / PASTE2** (Zeira et al. 2022; Liu et al. 2023) | FGW slice alignment and centre-slice integration; PASTE2 allows partial overlap | M (partial FGW is useful when the field of view differs, as in the heart) | raphael-group/paste, paste2 (BSD-3); PyPI `paste-bio` 1.4.0, `paste2` 1.0.1 | [10.1038/s41592-022-01459-6](https://doi.org/10.1038/s41592-022-01459-6), [10.1101/gr.277670.123](https://doi.org/10.1101/gr.277670.123) |
| **SLAT** (Xia et al., Nat Commun 2023) | Graph-adversarial matching of heterogeneous slices across technologies; spatio-temporal matching demo | L–M (alternative correspondence engine) | gao-lab/SLAT (MIT); PyPI `scSLAT` 0.3.0 (2024-12-20) | [10.1038/s41467-023-43105-5](https://doi.org/10.1038/s41467-023-43105-5) |
| **spaTrack** (Shen et al.; Cell Systems 2025 per its docs) | OT trajectories in space; weights expression against distance by spatial autocorrelation | L | yzf072/spaTrack (MIT) | preprint [10.1101/2023.09.04.556175](https://doi.org/10.1101/2023.09.04.556175); journal details not checked on the publisher page [Medium] |
| **FGW** (Vayer et al., ICML 2019) and **FGW barycenters** | Joint feature-and-structure OT. POT ships `fused_gromov_wasserstein[2]`, `entropic_*`, `semirelaxed_*`, `partial_*`, `quantized_fused_gromov_wasserstein*` (scales to large N) and `fgw_barycenters` | H (core coupling) | PythonOT/POT (MIT); PyPI `pot` 0.9.7.post1 (2026-07-29); ott-jax 0.6.0 (Apache-2.0) for GPU and low-rank GW (Scetbon et al., ICML 2022) | [arXiv:1805.09114](https://arxiv.org/abs/1805.09114); [arXiv:2106.01128](https://arxiv.org/abs/2106.01128) |
| **Displacement interpolation** (McCann 1997) | μ_t = push-forward of a coupling π under (x, y) ↦ (1−t)x + t·y | H (the E7.75 interpolant) | (formula) | [10.1006/aima.1997.1634](https://doi.org/10.1006/aima.1997.1634) |
| **novoSpaRc** (Nitzan et al., Nature 2019) | Places cells in a target shape via (F)GW, with smoothness of expression over space | M (assigning positions to cells inside a predicted shape) | rajewsky-lab/novosparc (BSD-3); PyPI 0.4.4 (2022) | [10.1038/s41586-019-1773-3](https://doi.org/10.1038/s41586-019-1773-3) |

FGW barycenters for **extrapolation** (weights outside [0, 1]) are not well defined. Forward McCann extrapolation x + t(T(x) − x) with t > 1 is only sensible for small overshoot. For a two-day heart extrapolation, use a *damped* version calibrated on pseudo-targets. [Medium]

### 3.2 Navigo: the probable `dynode_flow` reference

Navigo (Fan et al. 2026, bioRxiv [10.64898/2026.06.18.733286](https://doi.org/10.64898/2026.06.18.733286); code aristoteleo/Navigo-release, BSD-2, pushed 2026-06-10) is **expression-only**.
- It is rectified flow matching between consecutive time points, with an `MLPTimeGRN` whose velocity follows RNA kinetics (spliced/unspliced halves; α, β, γ).
- Per the preprint abstract, it is trained on a mouse atlas of about 12.4M cells and 43 time points. This is consistent with the Qiu et al. 2024 sci-RNA-seq3 time-lapse (12.4M nuclei, E8 to P0) [High for the architecture from code; Medium for the identity of the atlas].

Two consequences:
- **Pretrained Navigo checkpoints include stages inside the heart exclusion window (E9.5 < t ≤ E13.5). They are forbidden for heart E12.5 and for T1.** Retraining on allowed data is fine.
- MERFISH has no spliced/unspliced split. How `dynode_flow` adapts Navigo is [Unverified]. Its T2 coordinates are presumably copied.

### 3.3 Generative models for point clouds and spatial transcriptomics

| Model | Notes | Rel. |
|---|---|---|
| PointFlow (Yang et al., ICCV 2019) [arXiv:1906.12320](https://arxiv.org/abs/1906.12320); stevenygd/PointFlow (MIT) | Two-level CNF: shape latent, then points | L (we need one shape per stage, not a shape family) |
| Luo & Hu, CVPR 2021 [arXiv:2103.01458](https://arxiv.org/abs/2103.01458); luost26/diffusion-point-cloud (MIT, archived) | Point-cloud diffusion | L |
| PVD (Zhou et al., ICCV 2021) [arXiv:2104.03670](https://arxiv.org/abs/2104.03670); alexzhou907/PVD (MIT) | Point-voxel diffusion, shape completion | L |
| LION (Zeng et al., NeurIPS 2022) [arXiv:2210.06978](https://arxiv.org/abs/2210.06978); nv-tlabs/LION (NVIDIA licence) | Latent point diffusion | L (licence check needed) |
| OT-CFM / torchcfm (Tong et al., [arXiv:2302.00482](https://arxiv.org/abs/2302.00482)); PyPI `torchcfm` 1.0.7 | Minibatch-OT conditional flow matching; multi-marginal time-conditioned variants | M (optional learned joint (position, expression) flow for the embryo) |
| **Wasserstein Flow Matching** (Haviv, Pooladian, Pe'er, Amos; ICML 2025, PMLR 267:22238–22258) [arXiv:2411.00698](https://arxiv.org/abs/2411.00698); DoronHav/WassersteinFlowMatching (MIT) | Flow matching over *distributions* (point clouds) with attention; demonstrated on spatial-transcriptomics niches | M (niche-level generator) |
| **NicheFlow** (Sakalyan, Palma, Guerranti, Theis, Günnemann; NeurIPS 2025, OpenReview 5ofJyjgrth); kristiyansakalyan/nicheflow (MIT) | Variational flow matching of **microenvironments over time**, paired by OT on pooled coordinates and expression; includes a mouse embryonic development dataset; 2D (dx/dy grid) | M (closest published analogue of a "local-neighbourhood-preserving" temporal generator; needs a 3D port) |
| Mimyr (bioRxiv [10.1101/2025.11.24.690239](https://doi.org/10.1101/2025.11.24.690239); ICLR 2026) | Plane-conditioned diffusion of cell positions, then type classifier, then transformer for expression; fills missing slices | L–M (design pattern: positions first, then expression conditioned on geometry) |
| STITCH (bioRxiv [10.64898/2026.06.03.729557](https://doi.org/10.64898/2026.06.03.729557)) | OT-conditioned flow matching for 3D cross-slice gaps plus latent conditional FM for expression; no code link seen | L |
| stDiffusion (ICLR 2025 MLGenX workshop) | DDPM interpolating unseen slices | L |
| stDiff (Brief Bioinform 2024, [10.1093/bib/bbae171](https://doi.org/10.1093/bib/bbae171)); fdu-wangfeilab/stDiff (MIT) | scRNA-to-ST gene imputation | L |
| SpaDiT (Brief Bioinform 2024, [10.1093/bib/bbae571](https://doi.org/10.1093/bib/bbae571)); wenwenmin/SpaDiT (no licence file) | scRNA-conditioned DiT for ST gene prediction | L |
| CeLEry (Nat Commun 2023, [10.1038/s41467-023-39895-3](https://doi.org/10.1038/s41467-023-39895-3)); QihuangZhang/CeLEry (MIT) | Supervised location recovery for scRNA cells | M (place external or generated cells into our coordinate frame) |
| "scSpaceGen" | **Not found** on GitHub, PubMed or the web; may not exist under that name [Unverified] | — |

**Verdict.** No published model jointly generates realistic 3D tissue (positions plus 500-gene expression) for a held-out developmental stage at this quality level. Generic point-cloud generators solve a different problem (families of shapes). The best-supported building blocks are:
- (F)GW couplings (moscot or POT);
- non-rigid alignment (Spateo);
- displacement interpolation;
- real-cell resampling in count space.

Learned flow-matching models are an optional upgrade, to be justified only by local pseudo-validation.

### 3.4 Growth and morphogenesis models

| Model | Use for us |
|---|---|
| Le Garrec et al., eLife 2017 [10.7554/eLife.28951](https://doi.org/10.7554/eLife.28951): looping as buckling between fixed poles, with HREM 3D shapes | Explains the shape change from E8.5 to E9.5; qualitative only |
| Esteban et al., Nat Cardiovasc Res 2022 [10.1038/s44161-022-00065-1](https://doi.org/10.1038/s44161-022-00065-1): 3D+t atlas of heart tube formation | E7.75–E8.75 morphometrics. Respect the windows: heart data strictly between E8.25 and E8.75 is excluded, and the embryo window excludes E7.25 < t < E8.0 |
| de Boer et al., Dev Biol 2012 [10.1016/j.ydbio.2012.05.001](https://doi.org/10.1016/j.ydbio.2012.05.001): 3D atlas of proliferation and growth from heart formation to late gestation | Regional growth rates: proliferation drops on recruitment into the heart and then rises locally at sites of chamber formation. Use only values at ≤ E9.5 or > E13.5 for calibration |
| Spateo morphometric vector fields; CPD/BCPD-type non-rigid registration | Learn displacement fields between consecutive stages; extrapolate a **damped** field |
| Simple statistical shape model | With only 3–5 stages, PCA over a shape population is impossible. Model the trajectory of (log RMS, log principal-SD ratios) and apply an affine stretch in the canonical frame |

---

## 4. Biology needed for shape, size and composition (Q3)

### 4.1 Heart morphogenesis from E8.25 to E12.5 [High unless noted]

Sources: Savolainen et al. 2009 (atlas from E11.5 to E18.5, H&E, with a review of earlier stages), Krishnan et al. 2014 (EFIC from E9.5 to birth), Red-Horse et al. 2010 and de Boer et al. 2012.

- **E8.0–E8.5 (tube and onset of looping).**
  - The bilateral fields fuse into one beating heart tube by E8.0, with an arterial pole (aortic sac) and a venous pole (sinus horns). The heartbeat is regular by E9.0.
  - Rightward looping runs from E8.5 to E10.5. At E8.5 the bulbus cordis (future RV), the primitive LV and the common atrium are distinguishable.
  - The second heart field adds cells to both poles. [High]
- **E9.5–E10.5 (looped heart, early chambers).**
  - Endocardial cushions appear in the AV canal and OFT from E9.5.
  - The proepicardium invests the heart surface from E9.0 to E11.0; the **epicardium covers the heart by about E11.0**.
  - Septum primum starts at about E10.0.
  - The first and second pharyngeal arch arteries regress by E10.5; the third, fourth and sixth arch arteries form at E9.5–E10.5. [High]
- **E11–E11.5.**
  - **Coronary plexus begins at E11.5 on the dorsal surface, sprouting from the sinus venosus** (apelin-nlacZ). Intact E10.5 hearts have epicardium but no coronary sprouts yet.
  - Large ventricular trabeculations are present from E11.0.
  - The primary pulmonary vein canalises at about E11.
  - The muscular interventricular septum (IVS) begins at E11.5.
  - OFT cushions spiral at E11.5, and neural-crest-derived tissue divides the distal OFT. [High]
- **E12.5 (TS21).**
  - **The distal OFT is fully septated, including the semilunar valve primordia (thick, plate-like).** The proximal conal septum closes from distal to proximal, zipper-like.
  - Krishnan reports two completely distinct OFT channels by E12.5.
  - Separate right and left AV valves are forming (thick leaflets). The mesenchymal cap, vestibular spine and AV cushions fuse at about E12.0. Septum secundum starts at E12.5.
  - **The ventricles still communicate at E12.5**; the interventricular communication closes at E13.5–E14.0.
  - Coronary vessels begin invading the myocardium at E12.5–E14.5 and switch on arterial markers (Dll4, Hey1, ephrinB2) at E12.5.
  - [redacted 2026-10-19: A] Compaction follows at about E13–E14. [High]
- **Absolute dimensions.** I could not find a verified table of heart length or volume per stage for E9.5–E14.5 [Unverified]:
  - de Boer 2012 has quantitative volumes in an interactive supplement I could not access.
  - DMDD HREM at E14.5 (Geyer et al. 2017, [10.1111/joa.12663](https://doi.org/10.1111/joa.12663)) reports compact-wall thickness of [redacted 2026-10-19: A] across substages S21–S23, and interventricular communication in 97% of S21 embryos and 2% of S23 embryos.

  **Rule caution:** heart measurements at E10–E13.5 count as external data from the excluded window. Use only values at ≤ E9.5 or > E13.5.

**What a heart-scope MERFISH at E12.5 plausibly contains** [Medium; checked against the panel genes].

These panel genes **are** present: Myh7, Myl4, Myl7, Nppa, Nppb, Ttn, Hey2, Tbx5, Hand1/2, Gata4/6, Mef2c, Isl1, Bmp10, Cited1, Hcn4, Shox2, Tbx2/3, Gja1/5, Irx3, Pitx2, Tbx18, Wt1, Tcf21, Upk3b, Msln, Aldh1a2, Pdgfrb, Postn, Col1a1/1a2/3a1, Lum, Twist1, Msx2, Tgfb2, Vcan, Cdh5, Pecam1, Npr3, Nfatc1, Apln, Aplnr, Cldn5, Emcn, Gja4, Nr2f2, Plvap, Lyve1, Acta2, Myh11, Eln, Sox10, Foxd3, Plp1, Phox2b, Tfap2a/b, Afp, Apoa1, Sox2, Foxa2, Wnt2, Fgf10, Tbx4, Foxf1, Top2a, Cdk1.

These are **absent**: all globins, Myh6, Myl2, Nkx2-5, Tnnt2, Tbx20, Sox9, Kdr, Pdgfra, Nkx2-1, Sftpc.

Expected E12.5 content by population:
1. Ventricular CM split into **trabecular** (Nppa, Gja5, Irx3, Bmp10 high) and **compact** (Hey2, Top2a/Cdk1 proliferative) layers, plus LV/RV identity (Tbx5, Hand1 in the LV).
2. Atrial CM (Myl7, Nppa), and SAN and venous-pole cells (Shox2, Hcn4, Tbx3, Isl1).
3. AV canal and conduction cells (Tbx2/3, Bmp2, Gja5 in trabeculae).
4. Endocardium, with Npr3 in the chamber endocardium and Nfatc1 in the valve endocardium.
5. **Cushion and valve mesenchyme** (Postn, Twist1, Msx2, Tgfb2, Vcan).
6. **Epicardium** (Wt1, Tbx18, Upk3b, Msln, Aldh1a2) and **subepicardial EPDCs and fibroblast precursors** (Tcf21, Pdgfrb, Col1a1/3a1, Postn, Lum).
7. **Coronary endothelium** (Apln/Aplnr, Cdh5, Pecam1, Plvap, Emcn), concentrated on the dorsal surface and the AV groove.
8. Great-artery smooth muscle and NCC-derived OFT mesenchyme (Acta2, Myh11, Eln, Sox10, Foxd3, Plp1, Phox2b for cardiac ganglia).
9. Blood (no globin probes, so detected only weakly, via Gata1).
10. Depending on scope: liver (Afp, Apoa1), lung bud and foregut (Sox2, Foxa2, Wnt2, Fgf10, Tbx4, Foxf1) and pericardium or body wall.

Relative to E10.5, the robust *directions* are: populations 5, 6, 7, 8 and the compact/trabecular split go up, while early progenitor and pharyngeal-mesoderm signatures (Isl1+ SHF, Hand1 in the inner curvature and so on) go down.

### 4.2 Embryo gastrulation from E7.25 to E8.0 [Medium; staging approximate]

- **Staging.** Pijuan-Sala et al. 2019 describe their E6.5–E8.5 series as TS9–TS12: pre/early streak; mid/late streak; neural plate; headfold to somitogenesis. In Downs & Davies (1993) terms:
  - about E7.25 = late streak (LS), with no or early allantoic bud;
  - about E7.5 = early/late bud (EB/LB), with node and head process present;
  - about E7.75 = early headfold (EHF);
  - about E8.0 = late headfold, with 1–4 somites.
  - Per UNSW (secondary): the node is visible at TS10c; at TS11 the amniotic cavity is sealed, the neural plate and notochordal plate form and the headfolds enlarge; the foregut pocket begins at TS11d; at TS12a (E8.0) there are 1–4 somites, the foregut pocket is clear and the cardiogenic plate forms.
- **Events across the E7.5 → E8.0 window that change shape and composition:**
  - headfold elevation (anterior neural plate);
  - foregut pocket (anterior intestinal portal) invagination;
  - definitive endoderm completes intercalation into the visceral endoderm (by E7.5);
  - node and notochordal plate;
  - allantois elongation towards the chorion (contact at 5–7 somites, about E8.25);
  - yolk-sac blood islands (first primitive erythrocytes at about E7.5);
  - cardiac crescent moving towards the midline;
  - the first somites (paraxial mesoderm: Meox1, Tbx6 on the panel);
  - the posterior growth zone (Nkx1-2, Wnt3a, Fgf8 on the panel).
- **Growth.** Snow (1977): before the streak forms, the mean generation time is about 9 h; [redacted 2026-10-19: A] ([10.1242/dev.42.1.293](https://doi.org/10.1242/dev.42.1.293)).

  Our data show RMS [redacted 2026-10-19: A] from E7.25 to E7.5. [redacted 2026-10-19: A] The released E6.75, E7.25, E7.5 and E8.0 data will show the exact curve.
- **External data usable for the embryo board.** Exclude anything strictly between E7.25 and E8.0.
  - Pijuan-Sala atlas: E6.5–E7.25 and E8.0–E8.5 only. Drop E7.5, E7.75 and the "mixed gastrulation" samples.
  - Tyser et al. 2021 cardiac-crescent data, E7.5–E8.25: mostly excluded.
  - Lohoff et al. 2022 seqFISH at 8–12 somites: later than E8.0, allowed; but this is beyond the bracket.

---

## 5. Concrete pipelines (Q2)

### 5.0 Shared infrastructure (week 1)

1. **Local scorer replica.**
   - Use veckit from the GitHub `main` checkout (cc3859a), not PyPI 0.1.2.
   - Re-create the board protocol: 10% subsample of each stage; split the target in half; floor = copy_last at 10%; ceiling = other half; skill() as in `core_metrics.skill`; Table-1 weights; **gate = min(1, 2·s_nbhd)**.
   - Repeat over several seeds; the official seeds are unknown.
   - Validate the replica: the three validation boards' floor and ceiling values should approach the `index.json` anchors once the answers are released on 2026-10-20.
2. **Pseudo-targets** for model selection. There are only two official test submissions per board, and only the better counts.
   - Embryo:
     - E7.5 from {E6.75, E7.25, E8.0}, t = 1/3;
     - E7.25 from {E6.75, E7.5}, t = 1/2;
     - E7.5 from {E7.25, E8.0} without E6.75, for robustness.
   - Heart:
     - E10.5 from {≤ E9.5}, a 1-day extrapolation;
     - E9.5 from {≤ E8.75}, a 0.75-day extrapolation;
     - **E10.5 from {≤ E8.75}, a 1.75-day extrapolation and the closest analogue of the 2-day E12.5 task.**
3. **Count-space expression edits.** These preserve the sparsity that the variogram and mmd see.
   - Recover integer counts. Use the `counts` layer if present (the E9.5 sample has one). Otherwise invert log1p(1e4·c/Σc): the smallest non-zero value in a cell gives Σc.
   - To lower gene g by a factor f < 1 in a cell: binomial thinning, Bin(c, f).
   - To raise it: add Poisson(λ·(f−1)), where λ is the cell's count, or a type-specific rate when c = 0, so that detection rates also rise.
   - Re-normalise to 1e4 and apply log1p.
   - **Never add Δ to log values the way `pseudobulk_shift` does.** That fills zeros with small positive values, and is one reason that baseline sits at or below the floor.
4. **Cell-type harmonisation** across stages. The vocabularies differ by stage. Build a joint embedding on the 500 genes (Harmony or scVI) and transfer labels. Also define "core cardiac" versus "shell" groups for the heart (§2.2).
5. **Count rule**: n_pred = `truth_cells` from the test `index.json` if published, otherwise a best estimate (§1.3). No duplicate coordinates.

### 5.1 Embryo E7.75: FGW displacement interpolation with mixture expression [Medium/High expected value]

Inputs: E7.5 (released at P3) as A, E8.0 as B, plus E6.75 and E7.25 for curvature. Working sets of 3–10k cells per stage are enough, since the truth half is about [redacted 2026-10-19: A] cells.

1. **Canonicalise each stage.** Centre, apply a PCA rotation with det = +1, and divide by RMS. Store RMS_A and RMS_B, and the principal-SD ratios.
2. **Couple A and B with FGW.**
   - Feature cost: squared Euclidean distance in a joint 30–50-PC expression space fitted on A∪B.
   - Structure cost: within-stage canonical distance matrices.
   - α ≈ 0.3–0.6, small entropic ε, mild unbalancing (τ ≈ 0.95–0.99) so that proportions can change.
   - Tools: moscot `SpatioTemporalProblem(adata).prepare(time_key, spatial_key, joint_attr)` then `.solve(alpha, epsilon, tau_a, tau_b)`; or POT `entropic_fused_gromov_wasserstein` / `semirelaxed_fused_gromov_wasserstein` on GPU.
   - Optionally refine local correspondences with Spateo `morpho_align` (mode SN-N).
3. **Align frames.** Run weighted Procrustes under π with **proper rotations only**: R = argmin Σπ_ij‖R a_i + b − b_j‖², det R = +1.
4. **Interpolate positions.** Sample pairs (i, j) ~ π and set z = (1 − t)(R a_i + b) + t·b_j.
   - Nominal t = 0.5.
   - Calibrate a "developmental clock" on the E7.5 pseudo-target, where nominal t = 1/3: if t* ≠ 1/3, development is non-linear in clock time, so apply the same warp to 7.75.
5. **Interpolate expression.** Take cell i's profile with probability (1 − t), otherwise partner j's (a mixture of real cells). dp is then proportional to pb(B) − pb(A). Because the DE metrics are rank-based, this is right whenever gene trajectories are monotone.
   - **Optional curvature correction:** fit per-gene pseudobulk trajectories through E6.75, E7.25, E7.5 and E8.0 (quadratic, or PCHIP per type then mixed by composition). Then reweight the sampled pairs by exponential tilting, w ∝ exp(λᵀ·PC(x)), so that the weighted pseudobulk hits the curved target.
6. **Set the scale.** RMS(7.75) = exp of a PCHIP curve of log RMS at 7.75, through the four stages. Validate by leave-one-out at E7.25 and E7.5.
   - If principal-SD ratios differ between A and B, interpolate them in log space and apply the corresponding stretch in the canonical frame, before the final RMS rescale.
7. **Optional dice uniformisation.** A few Lloyd or repulsion iterations restricted to the occupied region, capped at about 0.3 × median NN distance. Keep only if pseudo-targets show d2 and nbhd unchanged.
8. **Export.** n_pred ≈ truth_cells; `.X` on the 498-gene test panel in the published order; `obsm["spatial_3D"]` in micrometres.

**Expected outcome.** On the E7.5 validation analogue, the leaders already sit at about 76–78 with dice and scale nearly free. Our edge on E7.75 should come from:
- a genuinely correct size (no anchor);
- better neighbourhoods (coupled real tissue at matched density);
- correct DE direction.

The bracket is shorter at test time (E7.5–E8.0) than for validation (E7.25–E8.0). [Low; validate]

**Second submission (hedge).** The same pipeline with the alternative clock (t* from calibration versus 0.5), or a "nearest real stage" variant (E8.0 tissue shrunk to RMS(7.75) with mixture expression). Pick by pseudo-target scores.

### 5.2 Heart E12.5: what makes it hard

- 2-day extrapolation (twice the validation gap).
- New populations (§4.1), with no stage at E10.5–E13.5 to bracket.
- An unknown imaged field of view. This drives both size and composition, through the "shell" of non-cardiac tissue.
- Constant-velocity extrapolation is below the floor on the validation analogue (47.8).

### 5.3 Heart E12.5: recommended pipeline [Medium expected value; highest-variance board]

**Size model.** Pre-register it from training stages before looking at anything else.

1. Compute RMS of the **core cardiac** cells (CM, endocardium, epicardium, cushion) and of the full cloud for E8.75, E9.5 and E10.5. Also compute the shell fraction.
2. Fit three candidate core models. These numbers are illustrative, from the full-cloud anchor-derived rates; refit on the core:
   - (i) **Log-linear at the latest rate** ([redacted 2026-10-19: A]): Δlog RMS(10.5 → 12.5) ≈ [redacted 2026-10-19: A].
   - (ii) **Exponentially decaying rate** fitted to the E8.75→E9.5 and E9.5→E10.5 rates ([redacted 2026-10-19: A]): Δ ≈ [redacted 2026-10-19: A].
   - (iii) **"Half-rate" rule.** [redacted 2026-10-19: A] Applied once over two days at E9.5→E10.5's rate: Δ ≈ [redacted 2026-10-19: A].
3. Choose by pseudo-target error, i.e. predict E10.5 from ≤ E8.75 (and E9.5 from ≤ E8.75).
4. Optionally add a legal external bracket: heart size at > E13.5 (for example from MOSTA E14.5–E16.5 sections or HREM at E14.5), interpolated with E10.5.
5. Combine with the scope decision below: whole cloud = core plus shell.

**Composition.**
- Harmonised types; centred log-ratio (CLR) trajectories over E8.75, E9.5 and E10.5, extrapolated with **damping**. Choose the damping factor on the E10.5-from-≤E8.75 pseudo-target. Expect it in 0.3–0.6.
- Explicit priors for emergent or expanding populations: epicardium/EPDC, coronary EC (dorsal surface first), valve mesenchyme, compact versus trabecular CM, great-artery SMC.
- **Scope hypothesis H1** ("continuation"): shell fraction follows the E9.5 → E10.5 trend.
- **Scope hypothesis H2** ("heart-focused"): shell reduced by about half or more. This is consistent with a larger heart filling the imaged field.

**Expression direction.**
1. Per-type shifts in count space (§5.0.3), damped.
2. Better: fix the direction from **allowed** later data. External E14.5–E16.5 heart cells are allowed for the extrapolation target because they are > E13.5 (with disclosure). Candidates: Qiu et al. 2024 sci-RNA-seq3 E8–P0, keeping only t ≤ E9.5 or t > E13.5; MOSTA Stereo-seq STDS0000058, E14.5–E16.5 sections; DeLaughter et al. 2016 E14.5/E18.5 [Unverified accession].
3. Map them into MERFISH space with gene-wise, per-type calibration learned on the competition's own paired modalities: **T1 scRNA E8.5, E9.5 and E10.5 (E10.5 released at P3) against T2 MERFISH E8.5, E9.5 and E10.5**. Same lab, same stages.
4. Then set the E12.5 per-type means and proportions to a time-weighted interpolation between E10.5 (own data) and E14.5+ (mapped), for example t = 0.5 in linear time.
5. Since DE is rank-based, **getting the sign and order of the E10.5 → later change right is what pays.** Magnitude does not.

**Tissue generation.**
1. Start from real E10.5 cells: about [redacted 2026-10-19: A], so sample without replacement up to n_pred.
2. For types that increase, add extra real cells of that type from the E10.5 pool, placed next to existing same-type cells: small offsets of about 1 NN distance, with no duplicates.
3. For new populations, place cells into anatomical niches defined on E10.5 geometry:
   - epicardium: outer surface shell of the myocardium;
   - EPDC: just beneath it;
   - coronary EC: subepicardial, dorsal side and AV groove first;
   - valve mesenchyme: inside the cushions;
   - SMC: OFT and great-artery walls.

   Take their expression from mapped external cells, or from E10.5 precursor cells edited in count space.
4. Delete cells of declining types.
5. Grow all coordinates with the chosen size model. Optionally apply a *damped* extrapolation of the E9.5→E10.5 non-rigid deformation field (Spateo morphofield) and of the principal-SD ratio trend.

**Submissions.**
- Submission 1 = H1 scope with its size.
- Submission 2 = H2 scope with its size.

The best of the two counts; this hedges the single largest uncertainty. Also check the handedness consistency of the released stages (§1.6). If it is ambiguous, let the two submissions differ in handedness as well. Dice is the only ranked term that is not reflection-invariant.

**Optional learned upgrade.** A time-conditioned microenvironment flow (NicheFlow-style, ported to 3D), trained on E8.25 → E8.75 → E9.5 → E10.5 and rolled forward two days. Accept it only if it beats the structured pipeline on the E10.5-from-≤E8.75 pseudo-target. [Low prior]

### 5.4 Per-metric checklist

| Metric | What maximises it here |
|---|---|
| de_score / de_direction | Correct direction of pb change across all genes. Mixture interpolation for E7.75; external-bracket direction for E12.5. Never negate; keep std(dp) ≥ 0.01·std(dt) |
| mmd_u / variogram | Real cells; count-space edits; correct composition; correct per-gene sparsity |
| d2_shape | Interpolated or extrapolated canonical form (displacement interpolation; damped deformation) |
| occupancy_dice | Proper rotations only; correct principal-axis order and ratios; on sparse boards, even point placement |
| scale_log_ratio | Pre-registered size model (embryo: monotone log-RMS curve; heart: core plus scope); never score-based bisection |
| neighborhood_mmd (25%, plus the gate) | Coupled real tissue; matched density (n ≈ truth_cells); niche-aware insertion; no duplicates |

### 5.5 Six-week plan (2–4 H100s; everything light except optional flows)

| Week | Embryo | Heart |
|---|---|---|
| Oct 19–25 | Register and download (downloadable heart T2 files total about [redacted 2026-10-19: A]: E8.25 225, E8.5 [redacted 2026-10-19: A], E8.75 84, E9.5 164, E10.5 [redacted 2026-10-19: A] MB; embryo file sizes not published); local scorer replica; copy and scaled-copy baselines; RMS curves | Same; type harmonisation; core/shell RMS table; handedness check |
| Oct 26–Nov 1 | FGW coupling (moscot/POT), displacement interpolation, mixture expression; pseudo-targets | Damped composition and shift baseline in count space; pseudo-target E10.5 from ≤ E8.75 |
| Nov 2–8 | Clock calibration, curvature tilting, dice uniformisation ablations | scRNA ↔ MERFISH calibration (T1 E8.5/E9.5/E10.5 against T2); external > E13.5 ingestion with exclusion log |
| Nov 9–15 | Freeze embryo submission 1 | Niche insertion; scope hypotheses H1/H2; size model selection |
| Nov 16–22 | Embryo submission 2 | Ablations; optional learned flow |
| Nov 23–Dec 2 | Final submissions (2 per board); method summary with data disclosures | — |

Compute: FGW on 5k×5k runs in seconds to minutes on one H100; quantized or low-rank FGW for 50k+ cells; flow-matching models are hours on one GPU. Run on the DGX, per the team's compute policy.

---

## 6. Rule risks (Human track)

1. **The PR #6 leak.** It contains E12.5 test-stage geometry statistics. Using them is held-out-data use. Do not use them, and notify the organisers.
2. **Test-board anchors.** If the test `index.json` publishes floor anchors, the `scale_log_ratio` floor reveals RMS(E10.5)/RMS(E12.5). Deriving our scale from it is "solving for target properties" from scores. Ask the organisers before using it; default to not using it.
3. **Score-based search.** Bisection or differencing on official test scores to set scale (or anything else) is forbidden. Choosing the better of two hedged submissions is allowed and automatic.
4. **Pretrained models.** Navigo checkpoints (and similar atlas-trained models) contain E9.5 < t ≤ E13.5 data, so they are forbidden for heart E12.5. Check every pretrained model's training stages.
5. **External data windows.**
   - Heart: exclude E9.5 < t ≤ E13.5. E9.5 exactly is fine; > E13.5 is fine with disclosure.
   - Embryo: exclude E7.25 < t < E8.0.
   - Heart interpolation window: E8.25 < t < E8.75.
   - Spateo's E11.5 and MOSTA E10.5–E13.5 sections are excluded for the heart.
   - Resources spanning a window need in-window cells removed and the removal documented.
   - Literature *measurements* at E10–E13.5 (e.g. heart volumes) are arguably external data from the window; do not use them to set predictions.
6. **Metadata.** Using `truth_cells` and `ref_cells` from `index.json` to match counts is published board metadata, not a score. It is probably fine, but ask publicly on Slack or Discord. Inferring cell counts from published file sizes is unreliable: T3 KO files store about 8 KB per cell versus 3.4 KB for E8.75 WT.
7. **Licences.** STORIES is GPL-3.0 (fine to use; disclose it). LION carries an NVIDIA licence. SpaDiT has no licence file. Keep code private; code sharing across teams is prohibited.

---

## 7. Open questions

- Is the test embryo reference E7.5 or E7.25? Check `ref_cells`.
- Will test boards publish `truth_cells` and the anchors? Will PR #6 (`size_fidelity`, `severity_ratio`) be adopted before or during P3?
- Is the heart imaging field at E12.5 heart-only or heart plus thoracic shell? This is the biggest driver of size and composition. It is unknowable from allowed data, hence the two-hypothesis hedge.
- Do released stages share one handedness convention?
- Does a calibrated scRNA → MERFISH mapping (T1 ↔ T2 at E8.5/E9.5/E10.5) transfer to external E14.5+ platforms such as sci-RNA-seq3 or 10x?
- How large is the count/density effect on real boards? Re-run §1.3 on E9.5 and E10.5 once they are downloaded.
- Does dice uniformisation generalise from sparse validation truths to the test truths?

---

## 8. References

Competition sources (all fetched 2026-10-18):
- https://virtualembryo.ai/challenge (tasks/spatial, data, evaluation, baselines, faq, rules)
- Live API: https://api.virtualembryo.ai/challenge/leaderboard?task=T2&board=T2:heart:val_extrap (and `:val_interp`, `T2:embryo:val_interp`); https://api.virtualembryo.ai/challenge/phase ; https://api.virtualembryo.ai/challenge/boards
- veckit: https://github.com/aristoteleo/veckit (commit cc3859a); PR #6: https://github.com/aristoteleo/veckit/pull/6 ; issues: https://github.com/aristoteleo/veckit/issues ; PyPI: https://pypi.org/project/veckit/
- Proposal and Nature Methods comment (local copies in `site_assets/`)

Methods:
1. Qiu X et al. Spatiotemporal modeling of molecular holograms. Cell 187:7351–7373 (2024). https://doi.org/10.1016/j.cell.2024.10.011 ; code https://github.com/aristoteleo/spateo-release
2. Klein D et al. Mapping cells through time and space with moscot. Nature 638:1065–1075 (2025). https://doi.org/10.1038/s41586-024-08453-2 ; https://github.com/theislab/moscot
3. Halmos P et al. DeST-OT. Cell Syst 16:101160 (2025). https://doi.org/10.1016/j.cels.2024.12.001 ; https://github.com/raphael-group/DeST_OT
4. Huizing GJ et al. STORIES. Nat Methods (2025). https://doi.org/10.1038/s41592-025-02855-4 ; https://github.com/cantinilab/stories
5. Zeira R et al. PASTE. Nat Methods 19:567–575 (2022). https://doi.org/10.1038/s41592-022-01459-6
6. Liu X, Zeira R, Raphael BJ. PASTE2. Genome Res 33:1124–1132 (2023). https://doi.org/10.1101/gr.277670.123
7. Xia CR et al. SLAT. Nat Commun 14:7236 (2023). https://doi.org/10.1038/s41467-023-43105-5
8. Vayer T et al. Optimal transport for structured data with application on graphs (FGW). ICML 2019. https://arxiv.org/abs/1805.09114
9. Scetbon M, Peyré G, Cuturi M. Linear-time Gromov-Wasserstein distances using low rank couplings and costs. ICML 2022. https://arxiv.org/abs/2106.01128
10. McCann RJ. A convexity principle for interacting gases. Adv Math 128:153–179 (1997). https://doi.org/10.1006/aima.1997.1634
11. POT. https://github.com/PythonOT/POT ; ott-jax https://github.com/ott-jax/ott
12. Nitzan M et al. Gene expression cartography (novoSpaRc). Nature 576:132–137 (2019). https://doi.org/10.1038/s41586-019-1773-3
13. Fan Y et al. Generative modeling of mouse embryogenesis for fate and disease prediction (Navigo). bioRxiv (2026). https://doi.org/10.64898/2026.06.18.733286 ; https://github.com/aristoteleo/Navigo-release
14. Yang G et al. PointFlow. ICCV 2019. https://arxiv.org/abs/1906.12320
15. Luo S, Hu W. Diffusion probabilistic models for 3D point cloud generation. CVPR 2021. https://arxiv.org/abs/2103.01458
16. Zhou L, Du Y, Wu J. Point-voxel diffusion. ICCV 2021. https://arxiv.org/abs/2104.03670
17. Zeng X et al. LION. NeurIPS 2022. https://arxiv.org/abs/2210.06978
18. Tong A et al. Minibatch OT conditional flow matching. https://arxiv.org/abs/2302.00482 ; https://github.com/atong01/conditional-flow-matching
19. Haviv D et al. Wasserstein Flow Matching. ICML 2025. https://arxiv.org/abs/2411.00698 ; https://github.com/DoronHav/WassersteinFlowMatching
20. Sakalyan K et al. NicheFlow. NeurIPS 2025. https://openreview.net/forum?id=5ofJyjgrth ; https://github.com/kristiyansakalyan/nicheflow
21. Mimyr. bioRxiv 2025. https://doi.org/10.1101/2025.11.24.690239
22. STITCH. bioRxiv 2026. https://doi.org/10.64898/2026.06.03.729557
23. Li K et al. stDiff. Brief Bioinform 25:bbae171 (2024). https://doi.org/10.1093/bib/bbae171
24. Li X et al. SpaDiT. Brief Bioinform 25:bbae571 (2024). https://doi.org/10.1093/bib/bbae571
25. Zhang Q et al. CeLEry. Nat Commun 14:4050 (2023). https://doi.org/10.1038/s41467-023-39895-3
26. Shen X et al. spaTrack. bioRxiv https://doi.org/10.1101/2023.09.04.556175 ; https://github.com/yzf072/spaTrack
27. Kern C et al. MERFISH+ … 3D human developing heart. bioRxiv (2025). https://doi.org/10.1101/2025.11.02.686137
28. Farah EN et al. Spatially organized cellular communities form the developing human heart. Nature 627:854–864 (2024). https://doi.org/10.1038/s41586-024-07171-z

Biology:
29. Savolainen SM, Foley JF, Elmore SA. Histology atlas of the developing mouse heart (E11.5–E18.5). Toxicol Pathol 37:395–414 (2009). https://doi.org/10.1177/0192623309335060
30. Krishnan A et al. A detailed comparison of mouse and human cardiac development. Pediatr Res 76:500–507 (2014). https://doi.org/10.1038/pr.2014.128
31. Red-Horse K et al. Coronary arteries form by developmental reprogramming of venous cells. Nature 464:549–553 (2010). https://doi.org/10.1038/nature08873
32. de Boer BA et al. Growth of the developing mouse heart: an interactive qualitative and quantitative 3D atlas. Dev Biol 368:203–213 (2012). https://doi.org/10.1016/j.ydbio.2012.05.001
33. Le Garrec JF et al. A predictive model of asymmetric morphogenesis from 3D reconstructions of mouse heart looping dynamics. eLife 6:e28951 (2017). https://doi.org/10.7554/eLife.28951
34. Esteban I et al. Pseudodynamic analysis of heart tube formation in the mouse. Nat Cardiovasc Res 1:504–517 (2022). https://doi.org/10.1038/s44161-022-00065-1
35. Geyer SH et al. Morphology, topology and dimensions of the heart and arteries … stages S21–S23 (E14.5). J Anat (2017). https://doi.org/10.1111/joa.12663
36. Tyser RCV et al. Characterization of a common progenitor pool of the epicardium and myocardium. Science 371:eabb2986 (2021). https://doi.org/10.1126/science.abb2986
37. Snow MHL. Gastrulation in the mouse: growth and regionalization of the epiblast. J Embryol Exp Morphol 42:293–303 (1977). https://doi.org/10.1242/dev.42.1.293
38. Downs KM, Davies T. Staging of gastrulating mouse embryos by morphological landmarks. Development 118:1255–1266 (1993). https://doi.org/10.1242/dev.118.4.1255
39. Kojima Y, Tam OH, Tam PPL. Timing of developmental events in the early mouse embryo. Semin Cell Dev Biol 34:65–75 (2014). https://doi.org/10.1016/j.semcdb.2014.06.010
40. McDole K et al. In toto imaging and reconstruction of post-implantation mouse development at the single-cell level. Cell 175:859–876 (2018). https://doi.org/10.1016/j.cell.2018.09.031
41. Pijuan-Sala B et al. A single-cell molecular map of mouse gastrulation and early organogenesis. Nature 566:490–495 (2019). https://doi.org/10.1038/s41586-019-0933-9
42. Lohoff T et al. Integration of spatial and single-cell transcriptomic data elucidates mouse organogenesis. Nat Biotechnol 40:74–85 (2022). https://doi.org/10.1038/s41587-021-01006-2
43. Chen A et al. MOSTA (Stereo-seq). Cell 185:1777–1792 (2022). https://doi.org/10.1016/j.cell.2022.04.003 ; data https://db.cngb.org/stomics/datasets/STDS0000058/
44. Qiu C et al. A single-cell time-lapse of mouse prenatal development from gastrula to birth. Nature 626:1084–1093 (2024). https://doi.org/10.1038/s41586-024-07069-w ; https://omg.gs.washington.edu/jax/public/about.html
45. UNSW Embryology, Mouse Timeline Detailed (secondary). https://embryology.med.unsw.edu.au/embryology/index.php/Mouse_Timeline_Detailed

Literature metadata for items 1–7, 10, 12, 23–25, 28–37, 41–42 and 44 was retrieved from PubMed with the DOIs listed above.
