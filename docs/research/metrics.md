# VEC metric playbook: how the official scorer behaves and what moves each ranked metric

Virtual Embryo Challenge (NeurIPS 2026), Human Team track. Compiled 2026-10-18 (test phase P3 opens 2026-10-20; final submissions due 2026-12-02).

Confidence labels: **[High]** = read in the scorer code, or verified numerically against published numbers. **[Medium]** = derived from code plus an assumption, or shown in a calibrated simulation. **[Low]** = plausible inference that has not been verified. Where something could not be verified, the text says so.

---

## 0. Executive summary

1. **The aggregation is now pinned down exactly, including the T2 structure gate [High].** I fitted every public row of the live leaderboard API (608 rows on 3 T2 boards; 203 on T1; 204 on T3; snapshot 2026-10-18 15:44 UTC). The fit holds to at most 1e-4 headline points:
   - T1 = 25 de_score + 25 de_direction + 30 mmd_u + 20 variogram.
   - T3 = 30 de_score + 25 de_direction + 25 |severity_slope| + 20 (0.6 mmd_u + 0.4 variogram).
   - T2 = 25 (½ de_score + ½ de_direction) + 25 (0.6 mmd_u + 0.4 variogram) + **25 · mean(d2_shape, occupancy_dice, |scale_log_ratio|) · min(1, 2·skill(neighborhood_mmd))** + 25 neighborhood_mmd.
   - Every coefficient above multiplies a skill in [0, 1].
   - The gate only bites when neighborhood_mmd is worse than the floor (skill < 0.5). The 17 rows on heart-extrap with skill below 0.5 rule out the no-gate and pure-multiplication alternatives.
   - A task's score is the mean over its boards. The overall standing is T1 + T2 + T3.
2. **Hyperbolic map, one-line intuition [High].**
   - skill = 1/(1 + x), where x is the fraction of the floor-to-ceiling gap still remaining. So: floor 0.5; halfway 0.667; a quarter of the gap left 0.8; ceiling 1.0. Twice as far from the ceiling as the floor gives 0.333.
   - For a distance metric, halving the floor's value is worth about +0.17 skill. For a correlation metric (floor 0), reaching half the ceiling is worth +0.17.
   - Each 0.1 of skill is worth 100 × weight points: 3.0 points for T1 mmd_u, 2.5 for T1 de_score or de_direction, 2.5 for neighborhood_mmd on each T2 board.
3. **Biggest point pools, measured against the current human leaders.**
   - T1: de_score headroom about 9 points (leader skill 0.64, raw about 0.37), de_direction about 6.5, mmd_u about 4.6, variogram about 3.9.
   - T2, per board: neighborhood_mmd is the single largest pool (about 5 to 8.5 board points of headroom). On heart-extrap, occupancy_dice is far from its ceiling for everyone (best raw dice about 0.84 against a ceiling of 0.9465).
   - T3: de_score and de_direction about 7 to 9 points each. Severity is nearly saturated (median team skill 0.92).
4. **Cheap points that are easy to miss.**
   - (a) **T3 |severity_slope| is almost free once the response direction is right [High].** The worst case is pinned at log(1e-3), which is exactly the floor value. A wrong-sign prediction therefore scores 0.5 (no loss). Any β > 0 with R² > 0.01 scores 0.6 to 1.0: β = 0.37 or 2.7 gives 0.87; β = 0.5 to 2 gives at least 0.91. That is up to +12.5 points.
   - (b) **T2 |scale_log_ratio|** is a single scalar, the RMS radius in absolute units. [redacted 2026-10-19: A]
   - (c) **T2 occupancy_dice** rewards even ("blue-noise") point placement. In my test, a jittered lattice scored dice 0.924 against 0.904 for random sampling of the same shape. That plausibly explains why 25 teams exceed the embryo dice ceiling (0.766).
   - (d) **Pseudobulk hygiene for de_direction (see item 5).**
5. **de_direction gets a large share of its value from reference-subsample noise [Medium, simulation + algebra].**
   - The server's reference is a 10% subsample. Both the true change (dt) and the predicted change (dp) contain the same term −ε_ref, provided the prediction is built from the full released reference.
   - In a simulation calibrated to the published T1 ceiling (0.79), a prediction with zero biological information still scored de_direction about 0.20.
   - Adding random noise on genes you know nothing about cost about 0.08. Leaving uncalibrated sampling noise in the submission's pseudobulk cost about 0.04.
   - A perfect, calibrated mean model reached 0.84, which beats the half-vs-half ceiling.
   - Recommendation: build deltas relative to the full released reference, make the submission's pseudobulk exactly the intended vector, and set Δ = 0 for unknown genes (no noise).
6. **de_score is a contest against an expression-level null.**
   - Only the order of dp at the two tails matters (top n_up and bottom n_dn, sizes set by the truth).
   - Rescaling is a no-op. A negated or noise-driven tail scores below 0. In the simulation, a no-change prediction scored −0.37 against a null of 0.35.
   - Only genes with reference level ≥ 0.25 can ever be down-hits.
   - If the model nominates too few confident genes, filling the tail with the expression-level prior in the orientation favoured by the training transitions helped (+0.08 to +0.13). The wrong orientation hurt (a change of −0.06 to −0.10).
7. **On 32,285 genes the variogram is mostly a per-gene sparsity statistic [High on samples].**
   - About 50% of random pairs have one gene silent. Those pairs carry about 55% of the score and equal exactly the difference in per-gene E[√x].
   - Shuffling every gene independently across cells (all co-expression destroyed, marginals kept) raised the variogram only 2% (32K genes) and 6% (500 genes), while mmd_u rose 2.4 to 2.9×.
   - kNN smoothing multiplied the variogram 7.4×. Adding a constant delta to every cell doubled it.
   - **Generate real-looking sparse cells, never smoothed or decoded means.**
8. **The structure gate and neighborhood_mmd reward keeping expression attached to positions.**
   - Shuffling expression over a correct-shape cloud raised neighborhood_mmd 18× (to random-cube level) while leaving every shape metric unchanged.
   - Jitter of 0.15 RMS tripled neighborhood_mmd and cut dice from 0.94 to 0.84.
9. **Evaluator.** A server-faithful local evaluator is specified in section 6 and implemented in `research/vec_local_eval.py` (smoke-tested on the sample files).
   - Unknowns: the exact 10% subsampling rule (RNG, seed, rounding), which half is the truth, and bake.py's seeds. The evaluator therefore averages several subsample/split replicates.
10. **Risk to watch: veckit PR #6 (open, not merged) would replace both "free" terms [High that it exists; Low that it goes live].**
    - severity_slope would become min(β, 1/β): magnitude would start to matter a lot.
    - scale_log_ratio would become size_fidelity: tissue radius and median nearest-neighbour spacing, so submitted density would start to matter.
    - The live board and index.json still use the 0.1.1 definitions.

---

## 1. Sources and what was verified

| Source | What it gave | Status |
|---|---|---|
| veckit 0.1.1 at `cc3859a` (local copy `ext/veckit`; github.com/aristoteleo/veckit, MIT): `common/core_metrics.py`, `common/shape_metrics.py`, `T{1,2,3}/metrics_v2.py`, `T2/metrics.py`, `score_h5ad.py` | Exact metric code, guards, seeds | Read in full **[High]**. GitHub main is still `cc3859a` (last code commit 2026-10-02). |
| Live leaderboard API (`api.virtualembryo.ai/challenge/leaderboard?task=…&board=…`, `/challenge/overall`, `/challenge/boards`) | Per-metric skill ×100 and headline for 200 rows per board, plus reference rows | Fetched 2026-10-18. Used to verify the weights and infer the gate **[High]**. |
| Site text, plus the site's own JS bundle (`/_nuxt/BPXbFRVP.js`, which holds the Evaluation "Scoring" tab) | Metric derivations, NaN rule, abs rule, the gate sentence, data/count rules | The derivation page never states the gate formula. The only mention is on the Reference-rows page ("shape credit multiplied by neighbourhood-structure skill"). The exact form below comes from the API fit. |
| Proposal (`site_assets/proposal.txt`) | Table 1 weights; eq. (1) to (10) | Matches the code. Exceptions: the proposal says S_overall = ⅓ Σ S_t, but the Rules and the API both use the sum; and the proposal calls severity "undefined" for β ≤ 0 or R² ≤ 0.01, while the code returns log(1e-3) (the index.json floor of −6.9078 confirms the code). |
| `site_assets/index.json` (= `/challenge/panels/index.json`) | Anchors, cell bounds, truth and reference sizes | Used for all point computations. |
| veckit PRs/issues #3, #6, #12 | #12: the server never renormalises; #6: proposed PSS/TSR change; #3: T1 is the 10x Multiome RNA release (~32k genes) | Read 2026-10-18 **[High]**. |
| My probes (section 5) on veckit's 150-cell samples plus synthetic data | Empirical behaviour | Each run took under about 2 min. Scripts are in `scratchpad/metrics_agent/work/`. |

---

## 2. The scoring pipeline, step by step

1. **Validation** (unlimited format checks) **[High]**.
   - `var_names` must equal the board panel element by element (T1 32,285; T2 embryo 498; T2 heart, T3 500).
   - X must be finite and non-negative, and is cast to float32. Sparse or dense is fine.
   - Cell count must be in [min_cells, max_cells], inclusive.
   - T2/T3 need `obsm["spatial_3D"]`, finite, with at least 3 columns (only the first 3 are read).
   - `obs["celltype"]` is ignored.
   - **No renormalisation:** a raw-count matrix is scored as if it were log1p-normalised (veckit PR #12). Submit log(1 + 1e4·c/Σc). In the released files, expm1(X).sum(1) = 10,000 exactly.
2. **Working copies** **[High for the fact; unknown for the exact RNG]**.
   - Every released stage (reference and target) is subsampled to 10%. Examples: E9.5 T1 has 17,057 cells, giving 1,706 (rounding, or ceil); the E9.5 heart reference is 5,374; the E8.75 WT is 2,483.
   - The target's 10% is split in half with `split_half` at a fixed seed: one half is the truth, the other half is the "ceiling" prediction. So **the truth is about 5% of the target's cells**: T1 val 1,249; T2 embryo [redacted 2026-10-19: A]; heart-interp [redacted 2026-10-19: A]; heart-extrap [redacted 2026-10-19: A]; T3 Gata4 2,320.
   - The board cell cap is max_cells = 3 × max(truth, ref), with a minimum cap of 5,000. Examples: 5,118 = 3×1,706; [redacted 2026-10-19: A]; 3,990 is raised to 5,000.
   - min_cells appears to be min(1000, truth_cells), e.g. [redacted 2026-10-19: A] on the embryo board **[Medium]**.
   - There is no sign that the submission itself is subsampled before scoring. The code scores it as submitted, and the data page says only the heavy metrics draw their own subsamples **[Medium]**.
3. **Metric panel**: veckit `score_task{1,2,3}_v2` functions, seed 0 **[High for the code; Medium that bake.py uses seed 0]**.
4. **Anchors per board**: floor = copy_last (T1/T2: the reference subsample itself, including its coordinates) or wt_identity (T3); ceiling = the other half. Both are scored by the same functions against the same truth half. Because the floor is the reference subsample itself, its de_score and de_direction are exactly 0 (the guard fires).
5. **Skill** = d_floor/(d_floor + d), clipped at 1, where d_floor = |ceiling − floor| and d = signed distance from the ceiling.
   - scale_log_ratio and severity_slope are scored on |value| as lower-is-better. I assume abs() is also applied to the anchors; the published anchors are signed and the site says "absolute value first" **[Medium]**.
   - A metric that is undefined or NaN gets skill 0 within its group (FAQ, leaderboard page) **[High]**.
6. **Groups and weights**, including the gate: section 0, item 1 **[High, fitted to the API]**.
7. **Board → task → overall**: task = mean over its boards; overall = T1 + T2 + T3; an unattempted task counts 0 **[High, Rules + API]**.
   - In P3 there are 2 scored submissions per board for the whole phase, and the better one counts.
   - Test boards: T1 E12.5; T2 embryo E7.75 (interpolation) plus heart E12.5 (extrapolation), **each worth half of T2**; T3 β-catenin. This makes one T2 test-board point worth 0.5 overall points.

---

## 3. The hyperbolic map and where the points are (question 6)

### 3.1 Generic shape
skill = 1/(1 + x), where x = d/d_floor is the remaining fraction of the floor-to-ceiling gap. In practice the floor and ceiling are usually far apart, so x is just the "fraction of the floor's error you still have".

| x (remaining fraction of the gap) | 3 | 2 | 1.5 | 1 (floor) | 0.75 | 0.5 | 0.25 | 0.1 | 0 (ceiling) |
|---|---|---|---|---|---|---|---|---|---|
| skill | 0.250 | 0.333 | 0.400 | **0.500** | 0.571 | 0.667 | 0.800 | 0.909 | **1.000** |

The map is convex towards the ceiling: the last 25% of the gap is worth as much as the first 50%. It is gentle below the floor (2× the floor's error costs only −0.17).

### 3.2 Board points for realistic raw values (published validation anchors)
"Points" means change in that board's 0 to 100 score relative to the floor. T1 and T3 points are overall points. **T2 points must be divided by the number of boards** (2 in the test phase) to get overall points. T2 shape terms assume the gate is open. Computed by `metrics_agent/work/points.py`.

**T1 (E10.5 validation; truth 1,249, ref 1,706)**

| raw value → | | | | | | |
|---|---|---|---|---|---|---|
| de_score (ceiling 0.846) | −0.2: −1.3 | 0: 0 | 0.2: +1.7 | 0.4: +3.9 | 0.6: +6.9 | 0.8: +11.2 |
| de_direction (ceiling 0.790) | −0.2: −1.4 | 0: 0 | 0.2: +1.8 | 0.4: +4.2 | 0.5: +5.8 | 0.6: +7.7 |
| mmd_u (floor 0.0836, ceiling 0.0041) | 2×floor: −5.2 | floor: 0 | ½: +5.3 | ¼ (0.021): +9.8 | ⅛: +12.8 | ceiling: +15.0 |
| variogram (floor 0.00522, ceiling 0.00016) | 2×floor: −3.4 | floor: 0 | ½: +3.5 | ¼: +6.3 | ⅛: +8.2 | ceiling: +10.0 |

**T3 (Gata4 validation; truth 2,320, WT 2,483)**

| | | | | | | |
|---|---|---|---|---|---|---|
| de_score (ceiling 0.870) | −0.2: −1.5 | 0.2: +1.9 | 0.4: +4.5 | 0.5: +6.1 | 0.6: +7.9 | 0.8: +12.8 |
| de_direction (ceiling 0.925) | −0.2: −1.2 | 0.2: +1.5 | 0.4: +3.4 | 0.5: +4.6 | 0.6: +6.0 | 0.8: +9.5 |
| \|severity\| = \|log β\| (floor 6.908 = wrong sign or no response) | 6.9 or wrong sign: 0 | 4.6 (β=0.01): +2.5 | 3: +4.9 | 2: +6.9 | 1: +9.4 | 0.5: +10.8 / 0: +12.5 |
| mmd_u (floor 0.0314) | 2×: −2.0 | ½: +2.0 | ¼: +3.7 | ⅛: +4.8 | | |
| variogram (floor 0.0322) | 2×: −1.4 | ½: +1.4 | ¼: +2.5 | ⅛: +3.2 | | |

**T2, per board** (divide by 2 for overall points in P3)

| metric | embryo interp (E7.5; truth [redacted 2026-10-19: A]) | heart extrap (E10.5; truth [redacted 2026-10-19: A]) | heart interp (E8.5; truth [redacted 2026-10-19: A]) |
|---|---|---|---|
| de_score 0.2 / 0.4 / 0.6 | +0.8 / +1.9 / +3.5 | +0.7 / +1.7 / +2.9 | +0.9 / +2.0 / +3.6 |
| de_direction 0.2 / 0.4 / 0.6 | +0.8 / +1.7 / +3.0 | +0.7 / +1.6 / +2.7 | +0.7 / +1.7 / +2.9 |
| mmd_u ½ / ¼ floor | +2.6 / +4.8 | +2.5 / +4.5 | +2.6 / +4.8 |
| variogram ½ / ¼ floor | +1.7 / +3.1 | +1.7 / +3.1 | +1.8 / +3.2 |
| d2_shape ½ / ¼ floor | +1.5 / +2.7 | +2.7 / +4.2 (gap 0.0093→0.0034 is narrow) | +1.7 / +3.2 |
| occupancy_dice: halfway / ceiling | +1.4 / +4.2 (floor 0.705, ceiling 0.766) | +1.4 / +4.2 (0.775 / 0.9465) | +1.4 / +4.2 (0.675 / 0.880) |
| \|scale\| = 0.10 / 0.05 / 0.02 | [redacted 2026-10-19: A] | [redacted 2026-10-19: A] | [redacted 2026-10-19: A] |
| neighborhood_mmd ½ / ¼ floor | +4.5 / +8.2 | +4.2 / +7.6 | +4.6 / +8.5 |

**Gate arithmetic.** With neighborhood_mmd skill s_nb < 0.5, the shape group is multiplied by 2·s_nb. For example, s_nb = 0.40 with mean shape skill 0.90 gives a shape group of 0.72 instead of 0.90, which is −4.5 board points on top of the neighborhood_mmd loss itself. Above the floor (s_nb ≥ 0.5) the gate has no effect.

### 3.3 What the field is reaching now (context only)
Raw levels implied by public skills (`inverse_lb.py`), for P2 validation boards:

| board | best human headline | current max raw (approx.) |
|---|---|---|
| T1 | 76.0 | de_score ≈ 0.37; de_direction ≈ 0.57; mmd_u ≈ 0.018; variogram ≈ 0.0011 |
| T2 embryo | 77.7 | de ≈ 0.33; dir ≈ 0.58; mmd ≈ 0.022; vario ≈ 0.024; d2 ≈ 0.0055; dice and scale at or above ceiling for many |
| T2 heart extrap | 67.3 | de ≈ 0.41; dir ≈ 0.54; mmd ≈ 0.009; vario ≈ 0.0175; dice ≈ 0.84 (ceiling 0.9465); nbhd ≈ 0.032 |
| T2 heart interp | 75.4 | de ≈ 0.44; dir ≈ 0.58; mmd ≈ 0.0059; dice ≈ 0.80; nbhd ≈ 0.018 |
| T3 Gata4 | 80.3 | de ≈ 0.57; dir ≈ 0.60; severity median \|log β\| ≈ 0.64; mmd ≈ 0.008; vario ≈ 0.0037 |

- Overall, the best human total is 224.7 (EmbryoForge), followed by 224.5. 171 human teams entered all three tasks; 150 of them clear the 150 floor.
- Published references: pseudobulk_shift is **below** the T1 floor (48.6). On T3, shift_transfer and gene_ko score 42.5: de_score raw ≈ −0.41 and variogram skill 0.155, because adding a delta to every cell is penalised hard.
- These inversions describe other teams' metric levels, not target properties. The Rules forbid using inverted scores to compute submission content (section 8).

### 3.4 Where the points are (leader headroom × weight, then adjusted for cost)
- **T1:** de_score (≈9 pts of headroom) > de_direction (≈6.5) > mmd_u (≈4.6) > variogram (≈3.9).
  - mmd_u is the most reliable large-delta lever for a team starting from the floor: ½ floor gives +5.3, ¼ floor gives +9.8.
- **T2, per board:** neighborhood_mmd (≈6 to 8.5) is largest everywhere. Then de_score (≈5) and, on heart extrapolation, dice (≈3.6) and mmd_u (≈5).
  - scale and d2 are mostly saturated by the leaders.
- **T3:** de_score (≈9) and de_direction (≈7.5) are the open pools. Severity is a must-get (+10 for nearly everyone). mmd/variogram are small (≤2.4).

---

## 4. Metric-by-metric behaviour (questions 1 to 5)

Notation: pb = mean over cells of log-normalised X. ref_sub = the server's 10% reference subsample. dp = pb(pred) − pb(ref_sub). dt = pb(truth) − pb(ref_sub).

### 4.1 de_score (T1 .25; T2 .125/board; T3 .30)
- **Construction [High].**
  - Truth DE genes are found by a two-sided Mann-Whitney U test with BH at 0.05 and |dt| ≥ 0.25. The docstring says BH removed zero genes on T1/T3, so in practice DE = {|dt| ≥ 0.25}.
  - The genes are split into n_up and n_dn.
  - The prediction contributes its top-n_up and bottom-n_dn genes ranked by dp (`np.argsort`, unstable among ties). Overlap is computed sign for sign.
  - Score = (ov − chance)/(1 − chance), where chance = max(overlap when ranking by +pb(ref_sub), overlap when ranking by −pb(ref_sub)).
  - Guard: std(dp) < 0.01·std(dt) gives score 0.
- **What moves it.** Only the ordering of dp at the two tails. Magnitude is irrelevant (positive rescaling is a no-op), and so is everything in the middle of the ranking.
- **Structural facts to exploit [High, algebra]:**
  - A gene can be a down-hit only if pb(ref_sub) ≥ 0.25, since log-expression cannot fall below 0. Never let low-expressed genes occupy the bottom tail.
  - An up-hit needs pb(truth) ≥ 0.25, so the gene's predicted target level must reach about 0.25.
  - Negating the prediction moves true-up genes into the down slot: the score can go negative.
- **Expression-magnitude null.** The best biology-free ranking (by level, either orientation) sets the origin. DE genes are enriched among highly expressed genes, so chance is usually well above n_DE/G.
  - In the simulation, a prediction carrying only noise in its tails scored −0.37.
  - If the model's confident list is shorter than n_up + n_dn, fill the remaining tail positions with a weak level-based prior. Use the orientation favoured on the training transitions (E8.5→E9.5 now; E9.5→E10.5 after the P3 release; Mab21l2 and Gata4 for T3). Simulation: +0.08 to +0.13 in the right orientation; a change of −0.06 to −0.10 in the wrong one. This is a second-order lever **[Medium]**.
- **How to maximise the expected score.**
  - Rank the up-tail by the posterior P(dt ≥ +0.25) and the down-tail by P(dt ≤ −0.25), not by the expected change.
  - Equivalently, shrink uncertain genes' predicted changes toward 0 so they stay out of the tails.
  - Make sure at least as many genes are confidently placed in each tail as the truth is likely to have. Estimate n_up and n_dn from training transitions; it is a few hundred to about 1,000+ on T1 [Low on the exact numbers].
- **Cheap vs expensive.**
  - Cheap: shrinkage; the down-tail eligibility rule; a level prior in the right orientation.
  - Expensive: actually knowing which genes cross ±0.25. On T1 that means predicting emerging cell types; on T3, the knockout's targets.
- **Failure modes.**
  - Sign errors (the T3 references score −0.41 raw).
  - Global or level-proportional deltas, which simply are the null.
  - Too few strongly predicted genes, so ties and noise fill the tails.
  - A pseudobulk computed in a different space than the scorer's (it must be the mean of log1p values).
- **Noise.** Truth DE calls near the 0.25 threshold flip with sampling. That is why the half-vs-half ceiling is only 0.82 to 0.94. Across simulated worlds a fixed strategy varied by about ±0.015 to 0.02 raw (about ±0.3 T1 points). The DE set itself is fixed by the truth half, so the board is deterministic for a given file.

### 4.2 de_direction (T1 .25; T2 .125/board; T3 .25)
- **Construction [High].**
  - The rank partial correlation of rank(dp) and rank(dt), controlling for rank(pb(ref_sub)), computed over **all** genes (32,285 in T1) with average ranks for ties.
  - It returns 0 for no change (the guard) and 0 if |corr(rank dp, rank ref)| > 1 − 1e-6 (any monotone function of the reference).
- **Ties [High on samples].**
  - Genes undetected in both the truth half and ref_sub have dt = 0 exactly. That was 43.9% of genes in the 75-vs-150-cell sample; the real proportion at 1,249/1,706 cells is unknown but probably still tens of percent.
  - These genes are also tied at the bottom of rank(ref), so after partialling they form one point mass.
  - Predicting exactly 0 for genes never detected in the full reference keeps that block concordant. Leaking tiny expression into them cost only 0.004 in simulation, so this is minor.
- **Pseudobulk noise, the important part [Medium: algebra plus simulation].**
  - Write pb(ref_sub) = pb(ref_full) + ε_ref. Then dt = Δ + ε_true − ε_ref, and for any prediction built on the full released reference, dp = Δ̂ + pb(ref_full) − pb(ref_sub) = Δ̂ − ε_ref.
  - Both vectors contain −ε_ref. For noise-dominated genes, with dp = −ε_ref exactly, corr(dp, dt) = sqrt(v_ref/(v_ref + v_true)) ≈ sqrt((0.9/1706)/((0.9/1706) + (0.95/1249))) ≈ 0.64.
  - In a simulation calibrated to the T1 ceiling (half-vs-half 0.78):

| prediction (T1 sim) | de_direction | de_score |
|---|---|---|
| no change, pb = full-reference mean (zero biology) | 0.204 | −0.367 |
| no change, a different 5,118-cell sample of the reference | 0.182 | −0.370 |
| perfect mean model, calibrated | **0.839** (above the ceiling) | 0.948 |
| perfect mean model, sampled 5,118 cells | 0.827 | 0.938 |
| perfect mean model applied to the server subsample instead (no shared noise) | 0.766 | 0.934 |
| knows the top 1,500 genes (ρ=0.8), rest Δ=0, full reference | 0.294 | 0.524 |
| same, rest Δ=0 relative to the server subsample | 0.255 | 0.521 |
| same + small noise on the unknown genes | 0.217 | 0.524 |
| same + 5,118-cell sampling noise in pb | 0.254 | 0.521 |
| knows the top 6,000 genes (ρ=0.8) | 0.406 | 0.500 |

- **How to set deltas for genes you know nothing about.**
  - (i) Δ̂ = 0 relative to the **full** released reference, not to any subsample of it.
  - (ii) No added noise: noise on unknown genes cost about 0.08 raw, about 2 T1 points.
  - (iii) Calibrate the submission so its pseudobulk equals the intended vector; uncalibrated sampling noise cost about 0.04.
  - (iv) Exact zeros for genes never detected in the full reference, unless you predict they switch on.
  - (v) Any part of a level-dependent trend that is linear in rank is partialled out and does nothing. Only nonlinear trends learned from training transitions can add signal **[Medium]**.
  - The ceiling is beatable (clipped at 100), because a smooth calibrated prediction lacks the other half's ε_true noise.
- **Cheap vs expensive.** The hygiene steps above are cheap (days). Real gains above about 0.3 raw need a gene-wide predictive model: composition shift plus within-type change.
- **Failure modes.**
  - Generator noise across many genes.
  - Deltas collinear with the reference level (scored 0).
  - Sign-flipped programmes, which give negative values (T3).
  - Subsampling your own reference before computing deltas, which removes the shared −ε_ref.

### 4.3 severity_slope (T3 only, .25)
- **Construction [High].**
  - β = dt·dp/dt·dt over the truth DE genes.
  - R² = 1 − ‖dp − βdt‖²/‖dp − mean(dp)‖².
  - The value is log β if β > 0 and R² > 0.01; otherwise log(1e-3) = −6.908, which is the published floor. With no response at all (the std guard), it is also −6.908. It is NaN only if the board has fewer than 5 DE genes.
  - It is scored on |value|.
- **Consequences [High].**
  - Getting the sign wrong costs nothing relative to the floor (skill 0.5).
  - The map is extremely flat near the ceiling: |log β| = 0.25 gives 0.966; 0.5 gives 0.934; 1 gives 0.874; 2 gives 0.776.
  - The only below-floor outcomes are 0 < β < 1e-3 or β > 1000 with R² > 0.01.
  - Because de_score and de_direction are scale-invariant, you can rescale the predicted response freely to target β ≈ 1, for example by undoing the shrinkage of an MSE-trained model, without affecting them. It does affect mmd_u and the variogram (weight 0.20) and the cells themselves.
  - A constant offset across the DE genes lowers R² (the denominator is centred).
- **Lever.** Get the response direction right on the large-effect genes; that alone gives +5 to +12 points. Then calibrate the magnitude on Mab21l2 (train) and Gata4 (validation answers released at P3).

### 4.4 mmd_u (T1 .30; T2 .15/board; T3 .12)
- **Construction [High].**
  - PCA with 30 components, fitted on **all truth cells** (unscaled log-expression, centred on the truth mean, randomized SVD with random_state=0).
  - Up to 2,000 cells drawn from the prediction and from the truth (seed 0); T1 truth is 1,249, so all of it is used.
  - Bandwidth γ0 = 1/median of the squared pairwise distances among the truth cells. Kernels exp(−s·γ0·d²) for s ∈ {0.25, 0.5, 1, 2, 4}; the unbiased MMD² is averaged over the five scales.
- **Composition vs within-type variation.**
  - The PCs are dominated by high-variance genes, i.e. markers of distinct types such as blood, hepatocyte and cardiomyocyte programmes.
  - The broad kernel (s = 0.25) mostly sees mixture weights and big shifts. The narrow kernel (s = 4, bandwidth half the median distance) sees within-type structure.
  - The floor (previous stage) on T1 is dominated by composition: half of E10.5's types are absent or renamed at E9.5.
  - In my 75-cell probe, collapsing within-type spread to type means raised mmd_u by 30% (and the variogram 4×). A 2× weight on one type changed it by 0.001 to 0.011, inside the noise at that tiny size. At real sizes composition errors are large relative to the ceiling (0.004) **[Medium]**.
- **Cheap vs expensive.**
  - Cheap: generate from real cells, never smoothed or decoded means. Reweight composition toward the predicted target proportions, and keep within-type heterogeneity.
  - Expensive: realistic expression states for cell types that appear for the first time.
- **Failure modes.**
  - Mode collapse or means (ctrl_one_cell).
  - Global scaling (×2 doubled the T1-sample mmd).
  - Out-of-distribution synthetic cells (for example, additive shifts clipped at 0 pile mass on zero).
  - Missing new types.
- **Noise.** The ceiling value is itself one noise draw: ≈0.004 at T1 sizes, 0.003 at [redacted 2026-10-19: A] cells, 0.0001 at 2,000 cells. At the leaders' level (0.018 on T1), ±0.003 of noise is about ±0.8 T1 points. The truth's PCA basis is fixed per board, so the noise comes from the finite samples.

### 4.5 variogram (T1 .20; T2 .10/board; T3 .08)
- **Construction [High].**
  - Up to 1,500 cells per side (seed 0), then 20,000 random gene pairs i ≠ j.
  - The gene pairs are drawn from the same RNG after the two cell draws, so **the pair set depends on n_pred and n_truth**.
  - Score = mean over pairs of (E_pred|x_i − x_j|^0.5 − E_truth|x_i − x_j|^0.5)².
- **What dominates on 32K genes [High, samples].**
  - "One gene silent" pairs: about 50% of pairs, 53 to 56% of the score, exactly the per-gene E[√x] difference (correlation 1.000).
  - "Both detected, weaker gene detected in <20% of cells": a further 30 to 36%.
  - Pairs where both genes are detected in >20% of cells: only about 12%.
  - So in practice it measures **per-gene detection rates and levels** (marginal distributions), weighted uniformly over genes; that is why low-expressed genes matter.
  - Destroying all co-expression by independent per-gene shuffles raised it only 2.1% (T1) and 6.4% (MERFISH 500 genes). At real sizes that is a skill change of about 1%.
  - It is exactly invariant to adding the same constant to every gene of every cell. It is not invariant to per-gene shifts or to scaling: ×2 multiplies each term by (√2 − 1)²·v².
- **Implications.**
  - Smoothing or decoding is heavily penalised: kNN smoothing gave 7.4×, type means 4×.
  - "Add the mean delta to every cell" doubled it, because it fills zeros. That is the pseudobulk_shift failure (T1 variogram skill 0.39) and the shift_transfer failure (T3 0.155).
  - Change means through **composition and detection** (which cells and how many express), not through additive constants.
  - Real-cell resampling or count-sampling decoders (NB or ZINB counts, then log1p normalisation) preserve sparsity.
- **Noise and count.** E[score] = (true difference)² + Var(v_pred) + Var(v_truth).
  - With fewer than 1,500 submitted cells, Var(v_pred) rises, which is a small systematic penalty. Always submit ≥ 1,500 cells; T1 val allows 5,118.
  - The variability from gene-pair sampling across seeds was CV 1 to 2%.
  - Duplicated cells (bootstrap) raised it 40% at 75 cells.

### 4.6 Shape group (T2 only; .25 × gate; ⅓ each)
- **d2_shape [High].**
  - 200,000 random within-cloud pairs; distances divided by the cloud's own RMS radius; 1-D Wasserstein between the two distance distributions.
  - Blind to translation, rotation, reflection and scale; essentially independent of count.
  - It is sensitive to the elongation profile, lobes, holes and bends.
  - Probes: anisotropic squash (4, 1, 0.25) raised it 68× (0.0034 → 0.227); jitter of 0.3 RMS 5×; random cube 3×; mirror 0×; isotropic scaling 0×.
  - The heart-extrap gap is narrow (0.0093 → 0.0034), so small form errors swing its skill a lot. Leaders reach the ceiling by deforming real neighbouring-stage clouds.
- **occupancy_dice [High].**
  - Both clouds are subsampled to min(n_pred, n_true), then canonicalised (centre; PCA rotation forced right-handed; divide by RMS), placed on a shared 16³ grid over [−3, 3] RMS, and compared by Dice. The best of 4 proper sign flips is taken; **axis permutations are not searched**.
  - Implications:
    - (a) Submit ≥ truth_cells. Otherwise the truth is subsampled down to your count and the expected Dice falls.
    - (b) The predicted cloud's principal-axis **ordering** must match the truth's. If two axes have similar variance and are swapped, overlap collapses.
    - (c) At these sizes the voxel edge is only a few nearest-neighbour spacings: 4.0 in my synthetic [redacted 2026-10-19: A]-cell test; the real boards were not measured. The code's own resolution check says read Dice only when this ratio is above 10. Dice is largely a sampling-density statistic, which is why the ceilings are 0.766 ([redacted 2026-10-19: A] cells), 0.88 ([redacted 2026-10-19: A]) and 0.9465 ([redacted 2026-10-19: A]).
    - (d) Even point placement helps: a jittered lattice gave 0.924 against 0.904 for random sampling, with d2 unchanged. A pure lattice hurt d2 (0.0135 against 0.0066). Real cells have a minimum spacing, so even spacing is physically plausible.
    - (e) Mirror images are not penalised (laterality blind spot).
- **scale_log_ratio [High].**
  - log(RMS radius pred / RMS radius truth) in absolute coordinate units, scored on |value|.
  - A pure scalar: rescaling the coordinates fixes it with no effect on any other metric (veckit PR #6 says the same).
  - [redacted 2026-10-19: A]
  - For interpolation, interpolate log RMS radius between the bracketing stages. For heart E12.5, extrapolate the log RMS trend over E8.25…E10.5. For reference, the anchors imply E10.5/E9.5 heart RMS ≈ [redacted 2026-10-19: A] and E7.5/E7.25 embryo ≈ [redacted 2026-10-19: A]; [redacted 2026-10-19: A].
  - Do **not** derive the size from test-board anchors or score differences (section 8).
- **The structure gate.** The shape group is multiplied by min(1, 2·s_nb). Shape credit is lost only if neighborhood_mmd falls below the floor. ctrl_random_cube is the example: the right volume, but neighborhood_mmd about 18× the ceiling.

### 4.7 neighborhood_mmd (T2 .25/board)
- **Construction [High].**
  - For every submitted cell, the mean expression of its 15 nearest spatial neighbours (**including itself**), computed on the **full submission**.
  - Then mmd_u (n = 2,000, seed 0) in a 30-PC space fitted on the truth's neighbourhood profiles.
- **Joint dependence [High, synthetic probes].**
  - It is the only ranked term that sees whether expression and position belong together.
  - Expression shuffled over a correct cloud raised it 18× (to random-cube level) with every shape metric unchanged.
  - Jitter of 0.05 / 0.15 / 0.3 RMS multiplied it by 1.1 / 2.9 / 7.6. A squash raised it only 1.7× (local relations largely survive).
- **Density, an open question.**
  - The truth is a 5% sample, so its 15-NN neighbourhoods span roughly 20× the volume of a true 15-cell neighbourhood (radius about 2.7×).
  - In a synthetic tissue with coarse domains, more submitted cells **lowered** neighborhood_mmd: at 0.25 / 0.5 / 1 / 2 / 4× the truth count it was 0.016 / 0.009 / 0.005 to 0.007 / 0.0033 / 0.0029. This is mainly estimator variance.
  - In real tissue with thin layers (endocardium, epicardium, myocardium) a denser submission makes neighbourhoods more local than the truth's and could hurt. **Test on the released validation answers at P3** [Low/unresolved].
- **Levers.**
  - Start from real tissue: warp the nearest stage's cloud and move each cell's expression with it.
  - For interpolation, blend bracketing stages with an OT or deformation map.
  - Keep local composition gradients. Never generate coordinates and expression independently.
  - The floor (real previous-stage tissue) is strong: 16 of 200 teams fell below it on heart-extrap and lost shape credit through the gate.

### 4.8 How subsampling, the half-split and fixed seeds shape the noise (question 1, cross-cutting)
- **The truth is about 5% of the target.** On small boards that is very few cells: [redacted 2026-10-19: A] for the embryo board. This drives the low Dice ceiling (0.766), the higher MMD ceilings, and DE-threshold flips.
  - T1 test (E12.5): the file is 360 MB against 647 MB for E10.5 (24,980 cells), so the truth half could be only about 700 cells [Low: estimated from file size].
- **The reference is a 10% subsample.** Its sampling error ε_ref enters both dp and dt (4.2). Floor de-metrics are exactly 0.
- **Every random choice uses seed 0, so a given file always gets the same score.** But the draws depend on population sizes: mmd_u's 2,000-cell draw; the variogram's 1,500 cells and gene pairs; occupancy_dice's matched subsample; d2's 200,000 pairs.
  - Changing n_pred or the cell order re-rolls these draws, a noise of order the ceiling values (for example ±0.8 T1 points through mmd_u).
  - With only 2 shots per board in P3 nobody can exploit this. Local evaluation must average over replicates.
- **The anchors are single draws.** A submission can beat the ceiling (values clip at 100): smooth calibrated pseudobulks on DE metrics, even placement on Dice, and the absence of an extra half's noise on mmd_u and the variogram.

---

## 5. Empirical probes (question 7)

All under about 2 minutes locally. Scripts in `scratchpad/metrics_agent/work/`; the samples are veckit's 150-cell files.

| probe (script) | setup | key numbers |
|---|---|---|
| Guards and attacks (`probe_t1.py`) | E8.5 sample as reference, E9.5 halves (75/75), 32,285 genes | copy_last, ×2, ×0.99 and one-cell: de 0, dir 0. Random direction: dir −0.008. Negated delta: dir −0.979. Reference + dt added to every cell: dir 0.998 but variogram 2× floor (0.0048 vs 0.0024). kNN smoothing of the ceiling half: variogram 0.0167 vs 0.0023 (7.4×), mmd 0.139 vs 0.091. Mean-tiling: mmd 0.223, variogram 0.0208. Only 16 DE genes at this size (BH binds here, unlike on the real boards). |
| Variogram decomposition (`probe_vario.py`) | Same | One-silent pairs 49 to 50% of pairs, 53 to 56% of score, corr 1.000 with the marginal E√x proxy. Weak-detection pairs (<20%) another 30 to 36%. Seed CV 1.1 to 1.6%. |
| Co-expression sensitivity (`probe_shuffle.py`) | Per-gene independent shuffle | T1: variogram +2.1%, mmd_u ×2.4. MERFISH: variogram +6.4%, mmd_u ×2.9. |
| DE simulation (`sim_de.py`, `sim_de_lib.py`, `sim_strat.py`, `sim_fallback2.py`) | Pseudobulk-level: binomial/hypergeometric detection with the server's 10% and half-split, 32,285 genes, effect size calibrated to the T1 ceiling de_direction (0.78; simulated de_score ceiling 0.93 vs 0.846 published, so the DE noise is somewhat understated) | Table in 4.2; prior-fill results in 4.1 **[Medium; synthetic]**. |
| Spatial (`probe_spatial.py`) | 14-domain ellipsoid tissue populated with real MERFISH cells; truth [redacted 2026-10-19: A] | Count, jitter, squash, cube, mirror and shuffle results in 4.6 and 4.7. |
| Dice placement (`probe_dice.py`) | Ellipsoid, truth [redacted 2026-10-19: A] random points, 10 replicates | Random 0.904±0.015; 5,000 points matched down 0.906; jittered lattice 0.924±0.009; lattice 0.913 (d2 doubled); shrunk core plus halo 0.889. |
| mmd composition (`probe_mmd.py`) | E9.5 halves, 4 splits | Type means only: mmd +30%, variogram ×4. Spread ×0.5: mmd +21%, variogram ×3.4. +0.1 on all genes: variogram unchanged (exact invariance), mmd +6%. |
| Gate and weights (`gate_fit.py`, `weights_fit.py`) | Public API rows | T2 gate min(1, 2·s_nb): max error 1e-4 on 608 rows; alternatives fail (no gate: 0.76; multiply by s_nb: 10.3). T1/T3 weights: max error 6e-5. |

---

## 6. Server-faithful local evaluator (question 8)

Reference implementation: `research/vec_local_eval.py` (uses veckit's own functions). Smoke-tested on the sample files; run it on the DGX for real data, since T1 densifies n × 32,285.

**Protocol**
1. **Inputs.** The full released target stage (validation answers arrive at P3), the full reference (T1/T2: the preceding stage; T3: the matched WT at the same stage), and the submission.
2. **Validate like the server.** Exact panel and order; finite, non-negative X (never rely on renormalisation); cell count in [min, max]; spatial_3D finite with at least 3 columns.
3. **For r = 1…R (R ≥ 10):**
   - Draw ref_sub = 10% of the reference and tgt_sub = 10% of the target (`rng.choice` without replacement, k = round(0.1·N)).
   - Split tgt_sub with `split_half` (permutation; first n//2 is one half): truth = one half, ceiling prediction = the other half.
4. **Metrics** with veckit 0.1.1 and seed 0:
   - de_score and de_direction against ref_sub.
   - mmd_unbiased(n=2000) and variogram_score(20,000 pairs, 1,500 cells).
   - T3: severity_slope.
   - T2: d2_distance, occupancy_dice, scale_log_ratio, and neighborhood_mmd (15-NN including self, computed on the full set, then mmd_unbiased n=2000).
5. **Anchors per replicate.** Floor = the metrics of ref_sub itself as the prediction (with its coordinates). Ceiling = the metrics of the other half. Check that the ceiling beats the floor for every metric; the script warns if not.
6. **Skill.** Hyperbolic map, clipped at 1. |·| for scale_log_ratio and severity_slope (applied to value, floor and ceiling). NaN gives 0.
7. **Aggregate.** Weights and gate exactly as in section 0. Board score ×100; task = mean over boards; overall = sum.
8. **Report** mean ± sd over replicates, and compare candidate models by **paired** differences on identical replicates.
   - Sanity checks: the floor should score 50 on set-based metrics. Local anchors should be distributed around the published index.json anchors for the validation boards; for example, T1 local ceilings should fall near 0.846 / 0.790 / 0.0041 / 0.00016.

**Unknowns to keep visible**
- The server's 10% rule (RNG, seed, round vs ceil, whether stages share a generator); which half is the truth; the seeds bake.py passes. `score_h5ad` defaults to 0.
- Whether the submission itself is ever subsampled. The code says no.
- The test-phase reference for the embryo E7.75 board: probably E7.5 once its answers are released (copy_last = "last observed stage"), but it could stay E7.25. Not confirmed.
- Test-board anchors and cell bounds (published only at P3).
- Whether veckit PR #6 (severity_ratio, size_fidelity) or other changes reach the live board before 2026-12-04. The leaderboard notes reference rows "move whenever the panel does".
- Whether the T3 truth pools both knockout replicates before the 10% draw [Low: probably pooled].

---

## 7. Per-task metric playbooks

Effort scale: S = days; M = 1 to 2 weeks; L = 3+ weeks. Points are versus the floor, from section 3.2.

### Task 1 (test: E12.5 from E8.5, E9.5, plus E10.5 released at P3)
**Required properties**
- 32,285 genes in released order; log1p(1e4) values; non-negative.
- About max_cells cells, and at least 2,000 (fills the mmd draw and avoids the variogram's variance penalty).
- Every cell a realistic sparse single cell, not a mean.

| rank | lever | metrics | expected gain | effort | risk |
|---|---|---|---|---|---|
| 1 | Pseudobulk hygiene: deltas relative to the full released E10.5 reference; calibrate the submission's pseudobulk exactly; Δ = 0 and no noise on unknown genes; exact zeros for genes never detected | de_direction (+0.04 to 0.10 raw, about +1 to +2.5), protects de_score | +1 to +3 | S | low |
| 2 | Real-cell generation with predicted composition: resample real cells, reweight types toward predicted E12.5 proportions, keep within-type spread; no additive shifts or smoothing | mmd_u, variogram | +5 to +12 (mmd ½→¼ floor: +5 to +10; variogram ½: +3.5) | M | medium (new types) |
| 3 | Gene-level DE model ranking by P(cross ±0.25); down-tail only for genes with ref ≥ 0.25; enough confident genes; level prior in the training-favoured orientation | de_score | +2 to +5 (0.2 to 0.4 raw) | M | medium |
| 4 | New and maturing cell states from external references after E13.5 (allowed with disclosure; E10.5–E13.5 excluded) | mmd_u, de_score | uncertain, possibly large | L | rule and disclosure care |

**Pitfalls**
- Decoded or smoothed means (variogram 7×, mmd up).
- Adding a delta to every cell (variogram 2×; pseudobulk_shift is below the floor).
- Global scaling.
- A raw-count or wrong-normalisation matrix (not renormalised by the server).
- Fewer than 1,500 cells.
- Building deltas on a subsample.
- A short confident list that leaves the de_score tails filled with noise.

### Task 2 (test: embryo E7.75 interpolation and heart E12.5 extrapolation; each worth ½ of T2)
**Required properties**
- Board panel and order: embryo 498 (Casp4 and Pnliprp1 absent); heart 500.
- spatial_3D in the same physical units as the truth.
- At least truth_cells cells, preferably 1 to 3× truth, but check density effects on real data.
- Expression and positions generated together.

| rank | lever | metrics | expected gain (board pts; overall ÷2) | effort | risk |
|---|---|---|---|---|---|
| 1 | Set the RMS radius from the stage trend (interpolate or extrapolate log RMS) | \|scale\| | [redacted 2026-10-19: A] | S | E12.5 heart extrapolation error: [redacted 2026-10-19: A] |
| 2 | Keep the gate open and build neighbourhoods from real tissue: warp the nearest or bracketing stage clouds and carry expression with the cells | neighborhood_mmd, plus the gate | +4 to +8.5 | M | medium |
| 3 | Even (blue-noise) placement inside the predicted volume; match the principal-axis ordering; submit ≥ truth_cells | occupancy_dice | +1.4 to +4.2 | S | low |
| 4 | Nonrigid growth deformation consistent with stage-to-stage aspect-ratio trends | d2_shape, dice | +1.5 to +4.2 | M | medium |
| 5 | Real-cell expression with composition interpolation or extrapolation | mmd_u, variogram | +2.5 to +7 | M | medium |
| 6 | DE hygiene and ranking as in T1 (500 genes, ceilings 0.92 to 0.99, so harder to convert) | de_score, de_direction | +1 to +4 | S–M | low |

**Pitfalls**
- MSE regression of coordinates (the frame is ignored; the relational structure is what is scored).
- Independent generation of coordinates and expression (closes the gate).
- Uniform random fill (ctrl_random_cube).
- Anisotropic squash.
- Wrong units.
- Swapped principal axes.
- Too few cells (Dice falls).
- Mirroring is neither penalised nor credited, so no point optimising it.
- Morans_I_agreement and composition_JSD are not ranked.

### Task 3 (test: β-catenin KO at E8.75 against the matched WT)
**Required properties**
- 500 genes. Submit the mutant itself, not the delta.
- Coordinates are required but **not scored**; use the WT coordinates.
- About 4,000 to 7,000 cells (truth about 2,300).

| rank | lever | metrics | expected gain | effort | risk |
|---|---|---|---|---|---|
| 1 | Get the response direction right on the large-effect genes; then rescale to β ≈ 1 (calibrate on Mab21l2 and Gata4) | \|severity\| | +5 to +12.5; no loss if wrong | S once a direction model exists | the β-catenin direction is uncertain ("hardest") |
| 2 | Lineage-restricted perturbation: change only Mesp1-lineage cells, through state or composition changes; leave other WT cells intact | mmd_u, variogram | +1.5 to +4; avoids the −3 to −4 of additive global shifts (shift_transfer variogram skill 0.155) | S–M | low |
| 3 | Identity and sign of the responding genes (wrong sign costs about 1.5 to 3.5 points) | de_score (.30), de_direction (.25) | +3 to +10 | M–L | high: published references get the Gata4 direction wrong |

**Pitfalls**
- Copying the training knockout's response ([redacted 2026-10-19: B]).
- Adding Δ to all cells.
- Responses collinear with WT expression level (scored 0).
- A constant offset across the DE genes (lowers R²).
- Absurd magnitudes (β > 1000 or < 1e-3).

---

## 8. Rule risks and open questions (summary)

**Rule risks**
- (1) Using score returns, published test anchors or inverted skills to compute target properties (for example the E12.5 RMS radius) is prohibited. Line-search by trial submissions is allowed, but P3 allows only 2 shots per board. If test anchors appear in index.json, ask the organisers before using any of them.
- (2) Reproducing the server's reference subsample to align noise on purpose could be read as probing the evaluation infrastructure. Avoid it. The natural full-reference effect in 4.2 needs no such step.
- (3) External data windows (T1/T2 heart: nothing after E9.5 up to and including E13.5) and the disclosure duty.
- (4) Metric definitions may change (PR #6).

**Open questions** are listed in the structured summary: the subsample and split seeds; the E7.75 test reference stage; test anchors and sizes; neighborhood_mmd's density dependence on real tissue; whether PR #6 goes live; the real n_up/n_dn and null orientation on the test boards.

---

## 9. References
- veckit (official local scorer), github.com/aristoteleo/veckit: https://github.com/aristoteleo/veckit. Commit used: `cc3859a` (2026-10-02).
- veckit PR #6, "PSS/TSR scoring change" (open): https://github.com/aristoteleo/veckit/pull/6
- veckit PR #12, "scorer never renormalises" (open): https://github.com/aristoteleo/veckit/pull/12
- veckit issue #3 (T1 gene space): https://github.com/aristoteleo/veckit/issues/3
- veckit issue #7 (canonicalisation handedness): https://github.com/aristoteleo/veckit/issues/7
- Challenge pages:
  - https://virtualembryo.ai/challenge/baselines
  - https://virtualembryo.ai/challenge/data
  - https://virtualembryo.ai/challenge/evaluation (Scoring tab text is shipped in https://virtualembryo.ai/_nuxt/BPXbFRVP.js)
  - https://virtualembryo.ai/challenge/faq
  - https://virtualembryo.ai/challenge/leaderboard
  - https://virtualembryo.ai/challenge/rules
- Board index with anchors: https://virtualembryo.ai/challenge/panels/index.json
- Leaderboard API (public JSON):
  - https://api.virtualembryo.ai/challenge/leaderboard?task=T1&board=T1%3Aval
  - https://api.virtualembryo.ai/challenge/leaderboard?task=T2&board=T2%3Aheart%3Aval_extrap
  - https://api.virtualembryo.ai/challenge/overall
  - https://api.virtualembryo.ai/challenge/boards
- Competition proposal (PDF): https://virtualembryo.ai/challenge/papers/neurips-2026-competition-proposal.pdf
- Scheuerer M, Hamill TM (2015). Variogram-based proper scoring rules for probabilistic forecasts of multivariate quantities. Monthly Weather Review 143:1321–1334. https://doi.org/10.1175/MWR-D-14-00269.1 (DOI verified via Crossref)
- Gretton A et al. (2012). A kernel two-sample test. JMLR 13:723–773. https://jmlr.org/papers/v13/gretton12a.html (page resolves)
- Osada R, Funkhouser T, Chazelle B, Dobkin D (2002). Shape distributions. ACM Transactions on Graphics 21(4):807–832. https://doi.org/10.1145/571647.571648 (DOI verified via Crossref)

**Local artefacts**
- Report: `research/metrics.md`
- Evaluator: `research/vec_local_eval.py`
- Probe scripts: `metrics_agent/work/*.py`
- API snapshots: `metrics_agent/api/*.json`
- Downloaded site JS: `metrics_agent/js/`

All under `SCRATCH/`.
