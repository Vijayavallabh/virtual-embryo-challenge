# VEC 2026 (Human Team): shared architecture, cross-task transfer, Generality Award and portfolio strategy

Compiled 2026-10-18 (two days before P3 opens on 2026-10-20). Sources: the scratchpad copies of the challenge site, the accepted proposal, the `aristoteleo/veckit` scorer, the **live public leaderboard API** (fetched 2026-10-18 15:44 UTC), GitHub and the literature. Confidence labels: **[High]** means verified from code, data or primary text; **[Medium]** means reasoned from verified facts but not tested; **[Low]** means judgement or extrapolation.

---

## 0. Executive summary

**Recommendation: build one shared "transport-then-adjust" core with thin task heads (Candidate B below) and run every board through it.** This is unified in architecture but tuned per task. Designed this way, it should cost almost no per-task score: the core is a superset of the per-task pipelines, so `copy_last`, constant-shift, bracket interpolation and the WT-plus-response baseline are all special settings of it. It also keeps the team eligible for the $8K Generality Award. Do not build a neural end-to-end generator (Candidate A) as the scoring path. It is the riskiest way to spend six weeks, because both distributional metrics reward realistic, sparse, real-cell-like profiles.

Seven findings drive this recommendation:

1. **Board arithmetic [High].** In the final phase the overall score is T1 + mean(T2-embryo, T2-heart) + T3, which I verified exactly against the leaderboard JSON. **A point on T1 or T3 is worth twice a point on either T2 board.**
2. **The skill map is a logistic in log-distance [High].** skill = 1/(1+x), where x is the distance to the ceiling divided by the floor's distance. Halving the distance gives 0.667; doubling it gives 0.333.
3. **The T2 "structure gate" is shape × min(1, 2·skill(neighbourhood_mmd)) [High].** This fits all 600 T2 leaderboard rows exactly. The gate only bites if neighbourhood_mmd is worse than the floor.
4. **The validation board is inflated by search.** Most top teams sit at 99-100 skill on T3 `severity_slope` and T2 `scale_log_ratio`, which is consistent with line search using 8 scored submissions a day. The test phase allows **2 blind submissions per board**, so test scores will be lower and noisier, and decided by how well blind predictions generalise. That matters most on T1 (new cell populations) and T3 (response direction).
5. **Rule landmines [High].**
   - **(a)** The open veckit PR #6 docstring tabulates geometry statistics measured on the held-out heart stages, including the E12.5 test target. Using them would be held-out-data use. This report does not reproduce them. **Notify the organisers.**
   - **(b)** The anchors in `index.json` (floor and ceiling values per metric) are statistics of the held-out targets. Do not use them to set prediction parameters.
   - **(c)** In P2, T1's E8.5 scRNA is measured data from a held-out T2-heart stage. Under the plain wording ("by any route") it must not feed any T2-heart submission. In P3 the point is moot because E8.5 heart MERFISH is released.
   - **(d)** [redacted 2026-10-19: B] studied the **exact test genotype** (MesP1-cre; β-catenin lox/lox). Using its reported marker changes is plausibly "data from another allele of the same gene at a comparable stage". Ask before using it.
   - **(e)** The Data page says external single-cell data is allowed only for T1. The Rules allow external data generally. This conflict needs an answer.
6. **The metric definitions may change before the test.** PR #6 proposes swapping `severity_slope` for `severity_ratio = min(β, 1/β)` and `scale_log_ratio` for `size_fidelity` (tissue RMS radius plus median nearest-neighbour spacing). Build for both: aim for β≈1 on T3, and emit individual-cell coordinates at realistic spacing on T2.
7. **Portfolio.**
   - **Effort:** T1 ~35%, T3 ~25%, T2-heart ~15%, T2-embryo ~10-15%, shared core and validation harness ~10-15%.
   - **Submission policy:** with best-of-two scoring, order matters only through information and development time.
     - First submission on each board: the best honestly validated candidate, in by ~Nov 18-22. An empty board scores 0, so this is also insurance.
     - Second submission: chosen by a rule written down beforehand.
     - On T3, the first score shows whether the predicted direction was right. Flip if it is below 50, refine if it is above.
   - **Expected value of the conditional T3 policy:** about +5 to +9 points over a single submission on a full-weight board, under the assumptions in section 4.4.

---

## 1. Verified facts that drive strategy

### 1.1 Exact scoring (re-derived from 600+ live leaderboard rows) [High]

| Task | Headline formula (skills in [0,1], ×100) | Check |
|---|---|---|
| T1 | 0.25·de_score + 0.25·de_direction + 0.30·mmd_u + 0.20·variogram | max abs error 0.0001 over 200 rows |
| T2 (per board) | 0.25·(0.5·de_score + 0.5·de_direction) + 0.25·(0.6·mmd_u + 0.4·variogram) + 0.25·**[mean(d2_shape, occupancy_dice, scale_log_ratio) × min(1, 2·s_nbhd)]** + 0.25·s_nbhd | exact on all 3 val boards (600 rows); a plain product gate is rejected (errors up to 0.41) |
| T3 | 0.30·de_score + 0.25·de_direction + 0.25·severity_slope + 0.20·(0.6·mmd_u + 0.4·variogram) | max abs error 0.0001 |
| Overall | T1 + (sum of T2 boards)/n_boards + T3; a missing board counts 0 | exact on top 100 rows |

Here s_nbhd is the skill on neighbourhood_mmd.

**Final-phase consequence [High]:** the test phase has one T1 board (E12.5), two T2 boards (embryo E7.75 interpolation, heart E12.5 extrapolation) and one T3 board (β-catenin). Each T2 board therefore carries **0.5 weight**.

Overall points gained per +0.10 skill on one metric in the final phase, from `scratchpad/scripts/portfolio_sim.py`:

| Metric | Overall points per +0.10 skill |
|---|---|
| T1 mmd_u, T3 de_score | 3.00 |
| T1 de_score, T1 de_direction, T3 de_direction, T3 severity | 2.50 |
| T1 variogram | 2.00 |
| T2 neighbourhood_mmd, per board (plus its gate role) | 1.25 |
| T3 mmd_u | 1.20 |
| T3 variogram | 0.80 |
| T2 mmd_u, per board | 0.75 |
| T2 de_score and de_direction, per board, each | 0.62 |
| T2 variogram, per board | 0.50 |
| Each T2 shape term, per board (gate open) | 0.42 |

### 1.2 Skill map behaviour [High]

skill = d_floor/(d_floor + d), clipped at 1, where d is the signed distance from the ceiling. In normalised form x = d/d_floor, so skill = 1/(1+x) = sigmoid(−ln x):

| x | 0 | 0.25 | 0.5 | 1 (floor) | 2 | 3 |
|---|---|---|---|---|---|---|
| skill | 1.00 | 0.80 | 0.667 | 0.50 | 0.333 | 0.25 |

- For DE metrics the floor value is 0, so x = 1 − v/ceiling. A wrong-signed response costs about as much as a right-signed one of the same size gains (v = ±0.5·ceiling gives 0.40 or 0.667).
- DE metrics are **rank-based and scale-free** for positive rescaling. The only magnitude gate is std(dp) ≥ 0.01·std(dt).
- **Consequence [Medium]:** a "floor-plus" submission works as a cheap, low-variance candidate. It is the last released stage's real cells plus a small but correctly signed pseudobulk shift. It earns DE credit while leaving mmd_u and variogram at floor level. The arithmetic for T1: about 54 if the direction is right, about 48 if wrong (§4).

### 1.3 Live validation leaderboard (P2), 2026-10-18 15:44 UTC [High]

- **Field:** 388 teams (313 human, 75 agent). 150 of 313 human teams exceed 150 total.
- **Human totals:** #1 224.74, #3 222.63, #6 (last prize slot) 220.82, #10 216.61, #20 207.35. The agent leader is at 220.66.
- **Board tops (any track):** T1 76.0; T2-embryo 77.7; T2-heart-extrap 67.8; T2-heart-interp 75.4; T3-Gata4 80.5.
- **Published baselines:** `pseudobulk_shift` T1 48.6, embryo 55.2, heart-extrap 47.8. `shift_transfer` and `gene_ko` on T3 are 42.5 (both below the floor).
- **Per-metric skill of the top 10 (×100):**
  - T1: de_score ≈61, de_direction ≈73, mmd_u ≈80, variogram ≈80.
  - T3: severity ≈99 (saturated), de_score ≈70, de_direction ≈72.
  - Heart-extrap: neighbourhood ≈63, occupancy ≈56, scale ≈100 (saturated).
- **Among the 50 best human teams,** T1 shows the largest spread (sd 5.3), then T3 (4.2), then T2 (3.3, diluted by averaging). Task scores are correlated at about 0.65-0.69.
- **The saturated columns (T3 severity, T2 scale) are what one-dimensional leaderboard search produces.** The organisers say so in PR #6: "a prediction whose response is 2x too strong still scores 91/100 … that column reads 99.x for essentially every submission" (see §1.5). Expect test scores below the validation numbers.

### 1.4 What is cheap or expensive per board (scorer mechanics) [High unless marked]

- **T1 (32,285 genes).**
  - de_direction is a rank partial correlation over **all** genes, so the sign of thousands of small, low-expression changes matters.
  - The variogram uses 20,000 random gene pairs, mostly low-expressed genes, so per-gene sparsity and level dominate it.
  - **Real-cell resampling or transport beats parametric decoders on both distributional metrics [Medium].**
- **T2.**
  - neighbourhood_mmd carries 25% directly and controls the gate. It compares 15-NN-averaged expression within each cloud separately.
  - **[Medium] Spatial density probably matters.** With 3× the truth's cell count in the same volume, each 15-NN neighbourhood is physically smaller, so neighbourhood averages are less mixed. That likely shifts the distribution. Test this in week 1 (§5).
- **T3.** Shape is not scored, so WT coordinates are fine. Only the expression response and the mutant-state distribution count.
- **Gate [High].** If neighbourhood skill is at least 0.5, shape credit is unaffected. So keep coordinates and expression coherent (cells move with their neighbours), then work on shape and size.

### 1.5 Metric-change risk and a held-out-data leak in veckit PR #6 [High]

- `aristoteleo/veckit` PR #6 ("feat(metrics)!: PSS/TSR scoring change, sync with vec_baselines, bump to 0.2.0") was opened 2026-08-25 by the same organiser account that answers scorer issues, and is still **open**. It proposes:
  - **`severity_ratio = min(β, 1/β)`**, ranked in place of the log-slope. Undershoot or overshoot by a factor of 2 gives 0.5; inverted or no response gives 0.
  - **`size_fidelity`** = mean of folded tissue-RMS ratio and folded median-NN-distance ratio, ranked in place of `scale_log_ratio`. `median_nn_distance` subsamples to **4,000 cells**, so submitted cell count interacts with the cell-scale term.
- **Status [Medium-High]:** the live board still behaves like the old definitions (saturated severity and scale columns).
- The timeline page says "Datasets and exact metric formulations can still be corrected." **Design for both versions** (§2.4, §4).
- **LEAK.** The PR docstring includes a per-stage table of heart RMS radius and median NN distance that **covers the held-out E8.5 and E10.5 and the E12.5 test target**. That is measured data from held-out stages. Even its qualitative conclusions about how size and spacing vary across stages are partly derived from those stages. Actions:
  1. Whoever builds the T2-heart size and growth model should not open that PR.
  2. Derive any size or spacing prior from released stages and permitted external data only. Commit the derivation to git before the first test submission.
  3. Email the organisers privately so they can redact it and decide how to handle it.

### 1.6 Other operational facts [High]

- **Registration:** we have not registered or downloaded data, and sign-in is required. This is the most urgent action.
- **Phase API:** `api.virtualembryo.ai/challenge/phase` currently returns `nominations_open: false`. The timeline describes Dec 2 as "the last moment a submission can be made **or a final entry nominated**", and the client code exposes `/challenge/teams/final`. **Expect a nomination step; the Generality claim and method summary will probably travel with it.**
- **Disclosure:** every external source must be disclosed **with the submission**, so s1 must already carry the disclosure text.
- **Limits:** format checks are unlimited; a rejected file does not consume an attempt; 1,200 MB per file.
- **Panel overlap (verified locally):** all 500 MERFISH panel genes appear verbatim in the 32,285-gene T1 list. The embryo panel is the heart panel minus Casp4 and Pnliprp1. Not on the panel: Ctnnb1, [redacted 2026-10-19: B]. On the panel: Gata4/6, Mab21l2, [redacted 2026-10-19: B].

---

## 2. Q1: candidate unified architectures

### 2.1 Rules for the Generality Award [High, text]

From Rules s.14 and the Prizes page:

- One $8K award across both tracks, held in addition to placed prizes.
- "enters at least two tasks and answers them with one shared architecture. **The same model architecture must underlie every task it enters.** Training it separately per task is allowed, and so are task-specific input and prediction heads … What does not qualify is an unrelated model per task."
- Judged on the method summary, which must say what is common, what differs and where the boundary falls.
- "**Among eligible entries the award follows final standing.**" The floor condition applies.

Two consequences:

- **(i)** Under the literal reading, if we enter all three tasks, the shared architecture must underlie **all three**, not just two. Ask (Q-E in §3.4).
- **(ii)** The award goes to the eligible entry with the best final standing, so it is only reachable as a by-product of a top overall finish. EV is P(top eligible entry across both tracks) × $8K. **[Low]** I estimate 5-25% conditional on a top-10 human finish, because most top teams probably run per-task pipelines; this is unverifiable.

The award only pays if the shared architecture costs us roughly nothing per task.

### 2.2 Three candidates

#### Candidate A: conditional latent flow over cell sets (neural, generative)

- **Shared:**
  - A latent space z (d≈32-64) aligned across modalities with gimVI- or ENVI-style modality-specific encoders and decoders. The 500 panel genes are shared reconstruction targets, plus adversarial or MMD alignment.
  - A conditional, non-autonomous velocity field v(z, t, c) plus a growth field g(z, t, c), trained by unbalanced conditional flow matching on OT couplings between consecutive stages (VGFM/TIGON-style). c is a setting or genotype embedding.
- **Heads:**
  - scRNA decoder (32K genes); MERFISH decoder (500).
  - Coordinate head: positions inherited from the coupling plus a displacement field and a global growth scale, or a point-cloud generator such as Wasserstein Flow Matching.
  - Genotype head: latent shift Δz(c_KO).
- **Feasibility (6 weeks, 2-4 H100):** compute is fine (VAE on ~60K scRNA cells takes tens of minutes per run; latent CFM takes minutes). The real cost is **4-6 person-weeks** to reach parity with simple baselines.
- **Main risks:**
  1. Parametric decoders produce over-smooth, dense profiles that lose on variogram (sparsity) and mmd_u. The fix is retrieval decoding, which collapses A into B.
  2. Very few transitions to learn conditional dynamics from: T1 has 3 stages in P3, T2-heart 5, T2-embryo 4, T3 2 KOs.
  3. The organisers' own neural ODE does not beat a constant shift on T1.
- **Per-task cost: high risk** of being worse on every board than specialised pipelines.

#### Candidate B (recommended): "transport-then-adjust" population operator

Shared core, the same code for every board:

1. **Encode** each modality: PCA or an optional shared ENVI-style latent; modality-specific.
2. **Couple** adjacent populations (consecutive stages, bracketing stages, or WT vs KO) with **entropic unbalanced OT**. For spatial data use fused Gromov-Wasserstein (moscot `TemporalProblem`/`SpatioTemporalProblem`, or ott-jax/POT). This yields cell-level couplings and growth or death marginals.
3. **Move:**
   - Interpolation: McCann displacement interpolation along the coupling.
   - Extrapolation: a small learned non-autonomous velocity-plus-growth MLP fitted by flow matching on the couplings, falling back to a constant velocity if it does not help.
4. **Reweight:** resample real source cells by growth weights to change composition.
5. **Decode by retrieval:** each moved cell takes a gene-space displacement averaged over its k nearest real anchor cells, scaled by a shrinkage α and clipped at 0. This keeps real sparsity.
6. **Birth module:** for populations that appear, sample from the nearest real cells of the best permitted source (the later bracketing stage for interpolation; external >E13.5 references for T1 extrapolation), weighted by the growth model.
7. **Calibrate** α, composition shrinkage β and the response scale λ by **leave-one-transition-out** across all released transitions, scored locally with veckit.

Task heads:

- **T1:** births from permitted external data (§3.2); no coordinates.
- **T2:** coordinate head. Each output cell inherits its source cell's position, then gets a displacement from the spatiotemporal coupling (or a Spateo-style morphometric vector field) and a global size factor taken from a pre-registered growth model. The cell count is chosen so that spatial density mimics the target's working copy (§1.4, §3.4 Q-H).
- **T3:** genotype head. A lineage-restricted response, Δx = −λ·m(cell)·u, where u is a WT-derived Wnt-responsive panel signature (from scRNA, §3.1) and m(cell) is the Mesp1-lineage × Wnt-activity weight. An optional composition term depletes SHF-derived states through growth weights. Coordinates are copied from WT.

Assessment:

- **Feasibility:** core about 1.5-2 person-weeks, each head about 1 person-week, harness about 1 person-week. Every component exists under permissive licences: moscot 0.5.2 BSD-3, ott-jax 0.6.0 Apache-2.0, POT 0.9.7.post1 MIT, torchcfm 1.0.7 MIT, geomloss 0.3.1 MIT. One H100 is enough; a second runs calibration sweeps. Work at 10-20% subsamples, which matches the scorer's own 10% working copy.
- **Per-task cost ≈ 0 [Medium]:** α=0 gives `copy_last`; constant velocity per type gives `pseudobulk_shift`; bracket mode gives interpolation; the genotype head gives WT plus a signed response. The best per-task pipelines are all reachable settings, so specialisation lives in heads and calibration.
- **Risk:** organisers might read it as "a framework of heuristics" rather than "an architecture". Mitigations: the shared learned velocity-growth network is the same network class on every task; ask Q-E; describe it as one population-dynamics operator conditioned on (time, setting, genotype).

#### Candidate C: conditional mixture-of-real-cells generator ("retrieval set model")

- **Shared:** a cell bank of all released cells plus permitted external cells, with a conditional weighting network w(cell | target condition) and a residual shift network. Meta-trained to reproduce each released stage from the others (leave-one-stage-out) with Sinkhorn or MMD losses.
- **Heads:** modality encoders; coordinates re-embedded from canonicalised bank positions with a growth scale; genotype condition.
- **Feasibility:** about 3 person-weeks.
- **Risks:** it can only reweight and shift existing states (fine for interpolation and T3, weak for T1 extrapolation without external data); there are few meta-training examples, so overfitting is a risk; mixing positions from different embryos needs care.
- **Per-task cost:** low to medium. **Its leave-one-out calibration is worth borrowing for B.**

### 2.3 Comparison

| | A: latent CFM generator | **B: transport-then-adjust** | C: retrieval set model |
|---|---|---|---|
| Shared component | encoder interface, conditional v/g field, decoders | OT coupling and solver, conditional v/g field, retrieval decoder, calibration protocol | cell bank, conditional weighting network, calibration |
| Task-specific | 2 decoders, coordinate generator, genotype embedding | encoders (32K vs 500), coordinate head, genotype head, birth source | encoders, coordinate re-embedding |
| Person-weeks to parity | 4-6 | **2-3** | 3 |
| GPU need | 1-2 H100, hours | 1 H100, minutes per run | 1 H100 |
| Variogram/MMD realism | weak unless retrieval decoding | **strong (real cells)** | strong |
| Can create new states | yes (unreliable) | via birth module and external data | no |
| Expected per-task cost | high risk | **≈0 (superset of baselines)** | low-medium |
| Generality story | strongest narrative | good: one operator, task heads | good |
| Main failure mode | over-smooth outputs; too few transitions | the learned field adds nothing (fall back to constant velocity) | overfits few stages |

### 2.4 Make B robust to both metric versions (live and PR #6)

- **T3 magnitude:** target β≈1, not just the right sign. Calibrate λ so the predicted response norm on likely-DE genes matches the typical observed KO response norm (Mab21l2 in P2, plus Gata4 in P3), scaled by the expected affected-lineage fraction.
  - Under the live log-slope, β=0.4 scores 0.88 severity skill.
  - Under `severity_ratio` it scores about 0.63, and β=1 scores about 1.0. Severity ceiling assumed at 0.95; see `portfolio_sim.py`.
- **T2 size:** emit **individual-cell** coordinates at realistic spacing. Grow the tissue by adding cells, not by uniformly rescaling a dense cloud. Choose the cell count with density in mind (both clouds above 4,000 cells, or a matched count; ask Q-H).

---

## 3. Q2: cross-modal transfer and the rule check

### 3.1 What scRNA (T1) can do for T2 and T3, ranked by expected value

1. **WT-derived Wnt-responsive signature for T3: highest value [Medium].**
   - Ctnnb1 and the classic Wnt targets ([redacted 2026-10-19: B]) are **not** on the panel but are in scRNA.
   - Recipe:
     - Score Wnt activity per scRNA cell (E8.5/E9.5, plus E10.5 in P3) from a target set.
     - Within mesoderm and SHF-like populations (pSHF, aSHF, OFT/RV-CM, JCF…), regress the 500 panel genes on Wnt activity, controlling for cell type, to get a signed "Wnt-responsive panel signature" u.
     - Transfer to WT E8.75 MERFISH by kNN or label transfer on the shared panel genes (Tangram, or a gimVI/ENVI shared latent).
     - Predict the β-catenin cKO response as −λ·u in cells with high Wnt activity and Mesp1-lineage identity.
   - It uses only WT data, so it is genotype-clean. It is correlational, so the direction is a prior, not a certainty.
   - **Do not** use DE lists from Ctnnb1-loss datasets or Wnt-loss phenocopies (rule).
   - CellOracle-style GRN simulation is an option, but β-catenin is a co-activator rather than a TF, and CellOracle's mouse base GRN needs checking.
2. **Composition and lineage priors for T2-heart and T3: moderate value [Low-Medium].**
   - The scRNA ladder adds cardiac-lineage composition trends (CM subtypes, endocardium, proepicardium, SHF) at stages that bracket the MERFISH stages.
   - Two caveats: the dissection scope differs (heart-centred scRNA vs heart-scope MERFISH), and the annotation vocabulary is **not harmonised** (Endocardium vs Endocardium-1/-2).
   - Use it as a prior on lineage ratios, not as absolute proportions.
3. **Shared latent for joint dynamics: moderate value, higher risk [Low].**
   - Learn temporal dynamics on the denser scRNA ladder and apply them to MERFISH latents (ENVI and gimVI use a conditional VAE with modality-specific decoders).
   - It strengthens the Generality story. Its score value is unproven and modality mismatch (panel signal vs scRNA dropout) is real.
4. **Imputation or denoising of MERFISH predictions: likely negative for score [Medium].**
   - Every metric compares against the raw, noisy MERFISH truth (mmd_u, variogram, neighbourhood averages).
   - Denoised predictions look *less* like the truth.
   - Use imputation only as an internal representation (to compute Wnt activity, or for kNN mapping), never in the submitted `.X`.

MERFISH → T1:

- Heart E8.25, E8.75 and E9.5 (plus T3 WT E8.75/E9.5) add time points between and at T1's E8.5 and E9.5 on 500 genes. They are useful as extra supervision for the shared velocity field.
- Spatial information is irrelevant to T1 scoring.
- The embryo-setting stages (E6.75-E8.0) have little relevance to T1.

### 3.2 External data levers per board (only as allowed)

- **T1 [High text, Medium value]:**
  - Anything after E13.5 is allowed with disclosure; after E9.5 up to and including E13.5 is excluded.
  - Bracketing E12.5 between internal E10.5 and external ≥E14 data turns extrapolation into interpolation for the birth module.
  - The mouse time-lapse atlas (Qiu et al. 2024, 12.4M nuclei, E8 to P0) qualifies **after removing nuclei staged after E9.5 through E13.5** and documenting the removal. Subset heart and surrounding lineages to mimic the heart-centred dissection.
- **T2-heart:** the same window applies. External spatial data ≥E14 is permitted by the Rules, but the Data page conflict (§3.4 Q-C) needs resolving first.
- **T2-embryo:** gastrulation scRNA (Pijuan-Sala 2019, E6.5-E8.5) would need everything strictly between E7.25 and E8.0 removed. It is only usable if Q-C resolves in favour of external single-cell data for T2.
- **Pretrained models:** forbidden if pretraining touched held-out stages or genotypes.
  - Nicheformer: SpatialCorpus-110M, human and mouse, 73 tissues. The corpus list is not in the repo; **unaudited**.
  - Navigo: requires spliced/unspliced layers, which the challenge data does not provide. Its tutorial checkpoints come from a "mouse embryogenesis" atlas whose stages are unverified.
  - **Default: no pretrained single-cell models for T1/T2-heart unless the corpus is audited.**

### 3.3 RULE CHECK: cross-task use of competition data

Governing text (Rules s.10):

> "No measured data from a held-out stage or genotype may be used, **by any route** … The held-out stages and genotypes are: E10.5 and E12.5 for Task 1; E7.5 and E7.75 in the embryo setting and **E8.5**, E10.5 and E12.5 in the heart setting for Task 2…"

The Task 1 page adds that E7.75 scRNA is withheld because "**measured data at a held-out stage is not released for any task, whichever assay produced it**." Yet E8.5 scRNA *was* released for T1, even though E8.5 is a held-out T2-heart stage. The site never says whether competition data from one task counts as "external data" for another.

| Use | P2 (now to Oct 20) | P3 (Oct 20 to Dec 2) | Recommendation |
|---|---|---|---|
| **T1 E8.5 scRNA → any T2-heart submission** (val_interp E8.5 or val_extrap E10.5; the rule bans *use*, not only predicting that stage) | **Prohibited on the plain reading** (measured data at a held-out heart stage, "by any route"). The release for T1 does not authorise use for T2. | E8.5 heart becomes training input ("every heart stage except E12.5"), so it is no longer held out. **Permitted in spirit; low risk.** | P2: exclude T1-E8.5 from anything producing T2-heart files, including the shared model's training set. Validation boards carry no prize weight, but use of held-out data is a disqualification ground regardless of phase. P3: allowed; confirm (Q-A). |
| T1 E9.5 scRNA → T2-heart or T3 | allowed (E9.5 is a released stage for both) | allowed | OK. |
| **T1 E10.5 scRNA (released P3) → T2-heart E12.5 test** | not available | **Ambiguous.** As released competition data it looks fine (T2's own E10.5 is training input). As "external" data it falls inside the stated "after E9.5 up to and including E13.5" exclusion. | Keep it out of T2-heart, or build both variants, until Q-B is answered. |
| **T2-heart E10.5 MERFISH (released P3) → T1 E12.5 test** | not available | **Same ambiguity** | Same. |
| T2-heart E8.25/E8.75/E9.5, T3 WT E8.75/E9.5, Mab21l2 KO E9.5 → T1 | allowed (all ≤E9.5; exactly E9.5 is explicitly permitted) | allowed | OK. |
| T2-embryo E6.75-E8.0 (plus E7.5 in P3) → T1 | allowed | allowed | Low value anyway. |
| T1 scRNA (E8.5+) → T2-embryo E7.75 | outside the E7.25-E8.0 window, but the Data page says "Task 1 is the only task where training on external public single-cell data is allowed" | same | Probably fine (it is competition data, not external public data); confirm (Q-C). |
| T1 scRNA → T3 β-catenin | WT data, no genotype issue | same | OK; confirm alongside Q-C. |
| Gata4 KO (released P3) → T3 test | not available | intended as training data | OK. |
| Published Ctnnb1-cKO gene-level findings in cardiac mesoderm (e.g. [redacted 2026-10-19: B], MesP1-cre; β-catenin lox/lox) [gene-level content redacted 2026-10-18] | — | **Borderline**: "another allele of the same gene at a comparable stage" | Ask (Q-F). Default: WT-derived signature only. |
| index.json anchors (floor/ceiling metric values) | statistics of held-out targets | same for test boards if published | **Never use them to set prediction content** (Q-G). |
| veckit PR #6 docstring table | held-out-stage measurements | same | **Never use; report it.** |

In the final phase specifically, the T1→T2-heart E8.5 question disappears: the E8.5 MERFISH answers are released anyway and the test target is E12.5. What remains is the **E10.5 cross-modal question in both directions**, because E10.5 sits inside the stated exclusion window for both extrapolation targets.

### 3.4 Draft questions for the organisers (send by email and post the non-sensitive ones in Slack #virtual-embryo-challenge)

> **Q-A (cross-task, P2).** T1's released `E8.5_RNA.h5ad` is measured at E8.5, which is a held-out validation stage of the Task 2 heart setting. The Task 1 page says measured data at a held-out stage is not released for any task, whichever assay produced it. During P2, may T1's E8.5 scRNA be used in any way (including as training data for a shared model) in a pipeline that produces Task 2 heart submissions (val_interp E8.5 or val_extrap E10.5)? Until you reply we are assuming it may not.
>
> **Q-B (cross-task, P3).** In P3 the heart E8.5/E10.5 MERFISH and the T1 E10.5 scRNA become training input. (a) May T1 scRNA at E8.5, E9.5 and E10.5 be used for the T2-heart E12.5 test and for T3? (b) May T2/T3 MERFISH at E8.25-E10.5, including heart E10.5, be used for the T1 E12.5 test? (c) Does the "after E9.5 up to and including E13.5" exclusion for extrapolation targets apply to competition data released in P3 from another task, or only to external data?
>
> **Q-C (external single-cell data outside T1).** The Data page says "Task 1 is the only task where training on external public single-cell data is allowed", while Rules s.10 permits external public data for all tasks with disclosure. Which governs? Specifically: (a) external scRNA for T2-embryo, with all cells strictly between E7.25 and E8.0 removed; (b) external spatial or scRNA data after E13.5 for T2-heart; (c) T1 competition scRNA used for T2 or T3.
>
> **Q-D (test-board references).** For each P3 test board, which stage is the reference, i.e. the `copy_last` floor and the DE reference? T1-E12.5: E10.5 or E9.5? T2-heart-E12.5: E10.5 or E9.5? T2-embryo-E7.75: E7.5 or E7.25? Will the test boards' `index.json` publish truth_cells, ref_cells and anchors?
>
> **Q-E (Generality scope).** If a team enters all three tasks, must the shared architecture underlie all three, or is it enough that at least two entered tasks share it while the third uses a different model? Does a shared population-transport operator count as "one shared architecture"? We mean OT coupling, a conditional velocity/growth network, a retrieval decoder and a common calibration protocol, with task-specific encoders, a coordinate head (T2) and a genotype head (T3).
>
> **Q-F (prior knowledge, T3).** Is published, qualitative knowledge of the Ctnnb1 loss-of-function phenotype in cardiac mesoderm treated as "data from the held-out knockouts … or a perturbation that phenocopies them"? An example is marker changes reported for MesP1-cre; β-catenin lox/lox embryos ([redacted 2026-10-19: B]). Are general Wnt-target gene sets (e.g., [redacted 2026-10-19: B]) and signatures derived from wild-type data acceptable?
>
> **Q-G (published anchors).** The anchors in `index.json` (floor/ceiling values per metric) are statistics of the held-out target; for example, the floor's scale_log_ratio encodes the target's size relative to the reference. May they be used to set prediction parameters, or only to interpret scores? We are assuming the latter.
>
> **Q-H (cell count and density).** Cell count is described as unscored, but neighbourhood_mmd (15-NN averages) depends on spatial density. If `size_fidelity` from veckit PR #6 goes live, the median-NN term (subsampled to 4,000 cells) also depends on count. May we choose our cell count using truth_cells from `index.json`? And will PR #6's `severity_ratio`/`size_fidelity` be used in the final phase?
>
> **Q-I (second submission).** Is it acceptable to choose the second official submission from a candidate set fixed before the first submission, for example the same T3 response with the opposite sign, based on the first submission's score? We read this as selection, which s.10 permits ("trying a parameter value, seeing how it scores, and trying another"). We will record the candidate set and hashes beforehand.
>
> **Q-J (disclosure: possible leak).** Separately and privately: the docstring of open veckit PR #6 tabulates heart RMS radius and median NN distance for E8.5, E10.5 and E12.5. We have not used these values and wanted to flag them.

---

## 4. Q3: portfolio across T1, T2-embryo, T2-heart and T3

### 4.1 Board economics (final phase)

| Board | Weight in overall | Max gain over floor (overall pts) | Validation evidence | Cost to beat 50 | Ceiling headroom | Variance |
|---|---|---|---|---|---|---|
| **T1 E12.5** | 1.0 | 50 | top 76; widest spread among top teams (sd 5.3) | **cheap to reach ~52-56** via floor-plus (§1.2); expensive beyond, because ~55% of the weight is distributional and needs new populations | **largest** | medium-high |
| **T3 β-catenin** | 1.0 | 50 | val top 80.5, inflated by searched severity; every published model has the direction wrong | **cheap if the direction prior is right**; asymmetric payoff under the live metric: right ≈65, wrong ≈46 | large | **high** (direction) |
| **T2-embryo E7.75** | 0.5 | 25 | val top 77.7; median of the top 200 is 66 (without val answers) | **cheapest**: bracketed interpolation (E7.5 in P3 and E8.0) | small after ×0.5 | low |
| **T2-heart E12.5** | 0.5 | 25 | val-extrap top 67.8 (scale saturated by search) | hard: 2-day extrapolation, the imaged field of view varies by stage, size is unpredictable | moderate after ×0.5 | medium |

### 4.2 Effort allocation (about 6 weeks)

- **T1 ≈35%:** highest weight and headroom, widest spread between teams, and the external-data lever (bracketing with ≥E14 references).
- **T3 ≈25%:** full weight, asymmetric payoff, and the two-submission sign policy.
- **T2-heart ≈15%:** neighbourhood structure first (gate plus 25%), then composition, then the size prior.
- **T2-embryo ≈10-15%:** bracket interpolation gets most of the available gain quickly; diminishing returns after that.
- **Shared core, validation harness and submission tooling ≈10-15%:** this serves every board and is what makes the team Generality-eligible.

For a 3-person team (≈18 person-weeks): core 2.5, T1 6.5, T3 4.5, T2-heart 2.5, T2-embryo 2. For 2 people: core 2, T1 4, T3 3, T2-heart 1.5, T2-embryo 1.5.

### 4.3 How to use the two submissions per board

With best-of-two scoring, **the order matters only through (a) information and (b) development time.** E[max(s1, s2)] is symmetric for a fixed pair. So:

- **s1 = the best honestly validated candidate**, frozen and submitted around **Nov 18-22** on all four boards.
  - It doubles as insurance: an empty board scores 0, and any T2 board missing halves T2's contribution.
  - Never submit a pure copy. Floor-plus dominates it in expectation.
- **s2 = picked by a decision rule written before s1** (git tag plus file hashes), submitted around **Nov 28-30**, never in the last hour.
  - **T3:** if s1 < 50 (direction probably wrong), submit the alternative-direction candidate. That is the sign-flipped response, or a second biological hypothesis such as composition-driven (SHF depletion) instead of expression-driven. If s1 ≥ 50, submit the magnitude-refined variant (β closer to 1).
  - **T1:** if s1 < 52, fall back to floor-plus (DE-only shift, smaller α). If s1 is 52-60, submit the improved version of the same approach (two more weeks of development). If s1 > 60, submit the bolder variant (more external-data births, larger α).
  - **T2 boards:** the same pattern. The fallback is floor-plus with the size prior; the upgrade is the improved coordinate and composition head.
- **Use public information.** Competitors' test scores are published as they are made, which shows what is achievable per board. That counts as selection information only.

**Illustrative value** (assumptions in `scratchpad/scripts/portfolio_sim.py`: right ~N(63,5), wrong ~N(46,2), flipped ~N(58,6)). Expected best T3 score:

| P(direction right) | One submission | Two independent tries, same hypothesis | Conditional flip policy |
|---|---|---|---|
| 0.3 | 51.1 | 52.7 | 59.8 |
| 0.5 | 54.5 | 56.5 | 61.0 |
| 0.7 | 57.9 | 60.2 | 62.2 |

That is about **+5 to +9 points on a weight-1 board** **[Low: depends on assumptions]**.

### 4.4 Expected total and prize line [Low]

- On the validation board the human prize line (6th) is 220.8. Test totals will be lower: blind predictions, harder targets (E12.5, β-catenin) and no searched scale or severity.
- A team scoring roughly T1 60-65, T2 62-68 (mean of the two boards) and T3 58-63 lands at roughly **180-196**. That is likely competitive for a placed prize if the field drops as expected. This is a judgement, not a forecast.
- Because the task scores of strong teams are correlated (about 0.65-0.69), **the differentiators will be the high-variance blind calls (T1 births, T3 direction), not the T2-embryo interpolation everybody can do.**

---

## 5. Recommendation and staged plan

**Unified versus per-task:** go unified in architecture (Candidate B: one population-transport operator, task heads, one repository with `core/` and thin `tasks/`), specialised in heads and calibration. Commit to the Generality claim only if decision point D3 passes.

| Date | Stage | Deliverables | Decision point or kill criterion |
|---|---|---|---|
| **Oct 18-20** | 0: setup | Register the team (sign-in) and accept the Rules individually; download training data on the DGX; uv environment (moscot, ott-jax, POT, scanpy, veckit 0.1.2); send Q-A to Q-J; email Q-J privately; brief the team on PR #6 and the anchors. | **D0:** if registration or download fails, escalate immediately. |
| **Oct 20-24** | 1a: P3 intake | Pull P3 `index.json` (test boards: ref stage, truth_cells, max_cells, anchors present?); check whether metric definitions changed (PR #6 merged? severity_ratio or size_fidelity live?); add validation answers; build the **honest-validation harness**: T1 E9.5→E10.5, T2-heart E9.5→E10.5, T2-embryo (E7.25,E8.0)→E7.5, T3 Mab21l2→Gata4, plus leave-one-transition-out. | **D1:** freeze which reference stage each test board uses (from `index.json` or the organiser answer). |
| **Oct 24-31** | 1b: core v0 | B with floor-plus, bracket mode and the WT-plus-response head; density experiment (score a released stage's half A, subsampled 0.25× to 4×, against half B, for neighbourhood_mmd and size terms); first internal scores for all 4 boards. | **D2 (Oct 31):** every board beats its floor on internal validation; otherwise stop feature work and fix it. |
| **Nov 1-14** | 2: heads | T1 births from permitted ≥E14 external data (removal of E9.5-E13.5 documented); learned non-autonomous field (kept only if it beats constant velocity in leave-one-out); T2 coordinate head and pre-registered size model (released stages plus permitted external morphometrics only); T3 WT-derived Wnt signature, λ calibrated so β≈1. | **D3 (Nov 10-14), unified vs specialised:** if B is within about 1 point of the best specialised variant on every board's honest validation, commit to unified and the Generality claim. Otherwise use specialised code where it wins, and drop the claim unless Q-E confirms that 2-of-3 eligibility is enough. |
| **Nov 15-22** | 3: s1 | Freeze s1 per board; write the s2 decision rule and candidate set (tag plus hashes); disclosure text for every external source (version, licence, stages, removal procedure); format checks; **submit s1 on all 4 boards.** | **D4:** read the four scores; apply the pre-written rule. |
| **Nov 23-30** | 4: s2 | Improve the chosen s2 candidates; submit s2 by Nov 30 at the latest; write the method summary (Generality section: what is shared, what differs, where the boundary falls). | — |
| **Dec 1-2** | buffer | Nominate final entries (watch `nominations_open`); archive every submitted file and its code hash (reproducibility, 14-day report). | — |

---

## 6. Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Held-out leakage through PR #6 values, index.json anchors, T1-E8.5 in P2, or pretrained models | medium | disqualification | §1.5 and §3.3 rules; corpus audits; pre-registered derivations; email the organisers |
| Metric change before the test (severity_ratio, size_fidelity) | medium | ±3-8 points on T3, ±1-2 per T2 board | §2.4 dual-robust design; re-check the scorer on Oct 20 |
| T3 direction wrong | high (base rate: all published references are wrong) | −4 to −19 points | WT-derived signature; conditional flip (Q-I) |
| Generality ruled out (B "not an architecture", or a 3-task requirement) | medium | lose the $8K option only | ask Q-E early; the plan does not depend on it |
| External-data scope ruling (Q-C) | medium | removes a T2 lever | T2 plan does not rely on external data |
| Licence: CellFlow (`cellflow-tools` 0.0.9) is PolyForm-Noncommercial | — | rules require permitting licences | prefer BSD/MIT/Apache tools |
| Service outage near the deadline | low-medium | lost board | s1 by Nov 22; s2 by Nov 30 |
| Over-fitting the released val transitions (only a few) | high | test underperforms validation | leave-one-transition-out; few hyperparameters; shrinkage |

---

## 7. Open questions (unresolved here)

1. Which stage is the copy_last/DE reference on each test board (Q-D)? Will test anchors and truth_cells be published?
2. Will PR #6's metrics be live in P3?
3. Does Generality eligibility require all entered tasks to share the architecture? Does an OT-based operator qualify (Q-E)?
4. May competition data cross tasks at E10.5 in P3 (Q-B)? May external single-cell data be used for T2 and T3 (Q-C)?
5. Is literature knowledge of the MesP1-cre; Ctnnb1 cKO phenotype "held-out genotype data" (Q-F)?
6. How strongly does submitted cell count (density) move neighbourhood_mmd? Answer with the week-1 experiment.
7. Are the test "inputs released without labels" anything beyond stage and genotype identity? Unclear from the site; resolves on Oct 20.
8. Nicheformer and other foundation-model corpora: do they contain mouse E10-E13.5 or Ctnnb1-perturbation cells? Not audited here.

---

## References (links verified 2026-10-18 unless noted)

**Challenge sources (local copies under `scratchpad/`)**
- Rules s.10, s.11, s.13, s.14: https://virtualembryo.ai/challenge/rules (`site_text/challenge_rules.txt`); Prizes: https://virtualembryo.ai/challenge/prizes; Data: https://virtualembryo.ai/challenge/data; FAQ: https://virtualembryo.ai/challenge/faq; Task 1 page (E7.75 statement): https://virtualembryo.ai/challenge/tasks/temporal; Timeline: https://virtualembryo.ai/challenge/timeline; Reference rows (structure gate): https://virtualembryo.ai/challenge/baselines
- Leaderboard API (public, read-only): https://api.virtualembryo.ai/challenge/overall, `/challenge/leaderboard?task=T1`, `/challenge/leaderboard?task=T2&board=T2:heart:val_extrap`, `/challenge/phase` (snapshots in `scratchpad/work/api/`)
- NeurIPS competition proposal (Table 1 weights): https://virtualembryo.ai/challenge/paper/neurips-competition (`site_assets/proposal.txt`)
- Cao, Lu & Qiu, "Towards predictive virtual embryos with genomics and AI", Nature Methods (2026): https://doi.org/10.1038/s41592-026-03055-4
- veckit scorer (MIT): https://github.com/aristoteleo/veckit; PR #6 (metric change; contains held-out statistics, see §1.5): https://github.com/aristoteleo/veckit/pull/6; issue #3: https://github.com/aristoteleo/veckit/issues/3; PyPI: https://pypi.org/project/veckit/ (0.1.2)
- Navigo (BSD-2-Clause; probably the `dynode_flow` reference): https://github.com/aristoteleo/Navigo-release

**Methods** (journal articles retrieved via PubMed; DOIs linked)
- Klein et al., "Mapping cells through time and space with moscot", Nature 638:1065-1075 (2025): https://doi.org/10.1038/s41586-024-08453-2; PyPI moscot 0.5.2 (BSD-3): https://pypi.org/project/moscot/
- Schiebinger et al., Waddington-OT, Cell 176:928-943 (2019): https://doi.org/10.1016/j.cell.2019.01.006
- Sha, Qiu, Zhou & Nie, TIGON (unbalanced dynamic OT, growth), Nat Mach Intell 6:25-39 (2024): https://doi.org/10.1038/s42256-023-00763-w
- Wang et al., Joint Velocity-Growth Flow Matching (VGFM), NeurIPS 2025: https://arxiv.org/abs/2505.13413
- Atanackovic et al., Meta Flow Matching, ICLR 2025: https://arxiv.org/abs/2408.14608
- Haviv, Pooladian, Pe'er & Amos, Wasserstein Flow Matching, ICML 2025 (PMLR 267:22238-22258): https://proceedings.mlr.press/v267/haviv25a.html
- Zhang et al., scNODE, Bioinformatics 40(Suppl 2):ii146-ii154 (2024): https://doi.org/10.1093/bioinformatics/btae393
- He et al., Squidiff (one diffusion model for development and perturbation), Nat Methods 23:65-77 (2026): https://doi.org/10.1038/s41592-025-02877-y
- Haviv et al., ENVI/COVET (joint scRNA–spatial conditional VAE, two decoders), Nat Biotechnol 43:269-280 (2025): https://doi.org/10.1038/s41587-024-02193-4
- Lopez et al., gimVI (modality-specific heads, shared latent): https://arxiv.org/abs/1905.02269; docs: https://docs.scvi-tools.org/en/stable/user_guide/models/gimvi.html
- Biancalani et al., Tangram, Nat Methods 18:1352-1362 (2021): https://doi.org/10.1038/s41592-021-01264-7
- Kamimoto et al., CellOracle, Nature 614:742-751 (2023): https://doi.org/10.1038/s41586-022-05688-9
- Tejada-Lapuerta, Schaar et al., Nicheformer, Nat Methods 22:2525-2538 (2025): https://doi.org/10.1038/s41592-025-02814-z; code (BSD-3): https://github.com/theislab/nicheformer
- Ahlmann-Eltze, Huber & Anders, "Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines", Nat Methods 22:1657-1661 (2025): https://doi.org/10.1038/s41592-025-02772-6
- Qiu X. et al., Spateo (3D E9.5/E11.5 atlas; morphometric vector fields), Cell 187:7351-7373 (2024): https://doi.org/10.1016/j.cell.2024.10.011

**Data and biology**
- Qiu C. et al., mouse time-lapse E8 to P0 (12.4M nuclei), Nature 626:1084-1093 (2024): https://doi.org/10.1038/s41586-024-07069-w
- Pijuan-Sala et al., gastrulation atlas E6.5-E8.5 (116,312 cells), Nature 566:490-495 (2019): https://doi.org/10.1038/s41586-019-0933-9
- [redacted 2026-10-19: B]

**Packages (PyPI JSON, 2026-10-18):** ott-jax 0.6.0 (Apache-2.0); POT 0.9.7.post1 (MIT); scvi-tools 1.5.1 (BSD-3, Python ≥3.12); torchcfm 1.0.7 (MIT); geomloss 0.3.1 (MIT); scanpy 1.12.4; anndata 0.13.4; cellflow-tools 0.0.9 (**PolyForm-Noncommercial-1.0.0**); spateo-release 1.1.1 (BSD).

**Not verified:** Nicheformer's corpus stage coverage; Navigo checkpoint training stages; whether PR #6 metrics will be live in P3; the reference stage on each test board; the exact neighbourhood_mmd sensitivity to cell count.
