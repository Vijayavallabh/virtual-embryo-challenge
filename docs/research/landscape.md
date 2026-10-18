# Virtual Embryo Challenge (VEC, NeurIPS 2026): competitors, organiser reference models and announcements

Research date: **2026-10-18** (leaderboard pulled 15:41 to 15:44 UTC; P3 starts 2026-10-20). Scope: Human Team track.
Confidence labels: **[High]** means checked directly against primary data (API JSON, repo code, the paper's own text). **[Medium]** means a well-supported inference, or a single secondary source. **[Low]** means a search-engine snippet or a speculative reading.

---

## 0. Executive summary

1. **The live board is public and needs no sign-in [High].** The leaderboard page is a Nuxt app. It calls `GET https://api.virtualembryo.ai/challenge/boards`, `/challenge/overall?track=&task=` and `/challenge/leaderboard?task=&board=&track=`. As of 2026-10-18, **388 teams** have at least one scored board: **313 Human, 75 Agent**. Of the 313 Human teams, 171 have entered all three tasks and 150 clear the 150-point prize-eligibility line.
2. **Where the top currently sits (validation, best of each team):**
   - T1: 76.0
   - T2 boards: heart-interp 75.4, heart-extrap 67.8, embryo-interp 77.7
   - T3: 80.5 (both tracks)
   - Best Human totals: **224.74** (EmbryoForge), 224.53 (Ivan & Arnav), 222.63, 221.93, 221.59, **220.82** (6th, the last prize slot), 216.61 (10th).
   - Per-task Human top-6 lines: T1 ≥ 73.7, T2 ≥ 69.3, T3 ≥ 78.2.
   - Between the 10-15 and 10-18 snapshots, the top moved ≤ 0.9 points per board, while the 10th place rose 0.05 to 2.3 points. The validation board is saturating at the top.
3. **The organisers' published models sit below or near the floor [High].**
   - On the live board, `pseudobulk_shift` scores 48.6 (T1), 47.8 (T2 heart-extrap) and 55.2 (T2 embryo, the only reference above 50).
   - `shift_transfer` and `gene_ko` both score about 42.5 on T3.
   - The ODE and flow references (`shift_ode`, `neural_ode`, `spatiotemporal_ode`, `dynode_flow`, `perturb_ode`) are described on the task pages but **are not on the public board, and no numbers are published for them**. The 8 adversarial `ctrl_*` rows are not exposed by the public API either.
   - The organisers' text says: the neural ODE "does not beat the one-line constant shift"; the spatiotemporal ODE's score "comes almost entirely from an analytic growth-rate extrapolation"; and on T3 "every [reference] model gets the direction of the response wrong".
4. **How the headline is built (re-derived from 1,000 public rows, [High]).**
   - T1 and T3 are exactly the published weighted sums.
   - On T2, the "structure gate" is exactly `shape_group × min(1, neighborhood_mmd_skill / 0.5)`. Shape credit is only cut when `neighborhood_mmd` is **below the floor**; this happens for 16 of 200 heart-extrap rows in the 15:41 UTC pull. The formula also fits archived snapshots from 10-15 and 10-18 exactly.
   - Saturated metrics: `scale_log_ratio` = 100 for 85, 109 and 41 of 200 rows on the three T2 boards; embryo `occupancy_dice` = 100 for 25 rows; T3 `severity_slope` ≥ 95 for 71 of 200 rows.
5. **Scorer versions are a trap [High].**
   - `pip install veckit` installs **0.1.2**, which is byte-identical to the unmerged PR #6. It reports `severity_ratio` and `size_fidelity`, which the board does **not** rank.
   - GitHub `main` @ `cc3859a` (0.1.1 plus the issue #7 right-handed-PCA fix) matches the board's documented metric set.
   - The scorer the server actually runs is not published.
   - PR #6 (open since 2026-08-25, authored from the organisers' own branch) would replace `severity_slope` and `scale_log_ratio` in the ranking. It shows the organisers consider both exploitable. A late switch is a real risk.
6. **dynode_flow is probably Navigo-like, but this is unverified [Medium/Low].**
   - Navigo (bioRxiv 10.64898/2026.06.18.733286; Qiu and Li labs, including veckit maintainer Zehua Zeng) is flow matching with reflow over spliced/unspliced RNA, plus RNA-kinetics velocity.
   - It does zero-shot knockouts by clamping genes during ODE rollout. This matches "reuses the Task 2 checkpoint through its perturbation path".
   - Against that: it needs `Ms`/`Mu` (spliced/unspliced) layers, which the challenge data (especially MERFISH) does not carry.
7. **What works for competitors (public journals and repos, [Medium]).**
   - Real-cell resampling and placement beats generated cells: OT-CFM, VAE and diffusion runs "never stably exceed copy-last" and blow up the variogram.
   - Decoupling helps: DE metrics depend only on the pseudobulk; distribution metrics depend on cells ("mean graft").
   - Geometry: set the tissue scale (free points), keep the shape intact, and get neighbourhood structure right.
   - T3: mix real KO-like cells into wild-type carriers and calibrate the response magnitude.
   - Some teams already use **late-stage external data (E14.5 Feng 2022 GSE193346, MOSTA E14.5)**. The rules allow it after E13.5, and it can partly bracket the E12.5 targets, subject to batch, dissection and platform differences. Teams using it are mid-table on validation (T1 61 to 66), so it is not yet a proven win.
8. **Rule and scoring announcements since launch [High].**
   - Rules §10 data-window restatement effective 26 Aug 2026.
   - The structure gate was added after `ctrl_random_cube` beat the floor.
   - veckit PR #9 merged on 2026-10-02.
   - FAQ and Data pages updated on 2026-10-16: the outright E9.5 to E13.5 exclusion window covers **Task 1 and Task 2 heart** extrapolation, and the cell minimum is "the board's own minimum".
   - No X/Twitter posts by organisers could be retrieved.

---

## 1. Live validation leaderboard (public API)

### 1.1 How the data was obtained [High]
- `https://virtualembryo.ai/challenge/leaderboard` is a Nuxt 3 app. It is not Next.js, so the bundles live under `/_nuxt/`, not `/_next/static/`. Runtime config: `apiBase: "https://api.virtualembryo.ai"`.
- In the API client (bundle `DxvW8jU0.js`), the helper `q()` adds `Authorization: Bearer` **only if the visitor is signed in**. Anonymous visitors call the same endpoints without auth. The leaderboard component (`DivB5mT6.js`) calls:
  - `GET /challenge/boards` → `{"split":"val","boards":[5 keys]}`
  - `GET /challenge/overall?track={human|agent}&task={T1|T2|T3}` (both optional) → per-team per-board and per-task scores, `total`, `rank`, `reference`, `baseline`, `counts`
  - `GET /challenge/leaderboard?task=T1&board=T1:val&track=human` → per-metric skills, `model` label, `scored_at`
  - Also used: `GET /challenge/phase` → `{"phase":"p2","daily_quota":8,"split":"val","accepts_submissions":true,...}`. `GET /challenge/community/register` is called by the public Community-resources page.
- Raw JSON (fetched 2026-10-18 ~15:41 to 15:44 UTC) is saved under `SCRATCH/landscape_work/api/`.
- **The per-board endpoint is capped at 200 rows.** `counts.teams` is 200 on every board, while `overall` lists more teams. A community study reports that `limit/page/offset` have no effect. Top-50 statistics are unaffected.
- The API host was `kg.virtualembryo.ai` until early October; it had a 502 outage on 2026-10-03 (veckit issue #10). Data manifest and download links (`/challenge/data/manifest`, `/challenge/data/link?key=`) require sign-in and were not touched.

### 1.2 Field size [High]

| | Human | Agent | Total |
|---|---|---|---|
| Teams with ≥1 scored board (overall) | 313 | 75 | 388 |
| T1:val | 269 | 66 | 335 |
| T2:heart:val_interp | 207 | 61 | 268 |
| T2:heart:val_extrap | 207 | 61 | 268 |
| T2:embryo:val_interp | 211 | 62 | 273 |
| T3:gata4 | 200 | 60 | 260 |

Human teams by task combination: all three, 171; T1 only, 72; T2 only, 23; T1+T2, 18; T3 only, 11; T2+T3, 10; T1+T3, 8.

Human teams above the floor, per board:

| Board | Above 50 |
|---|---|
| T1 | 169 / 269 |
| Heart-interp | 197 / 207 |
| Heart-extrap | 151 / 207 |
| Embryo | 207 / 211 |
| T3 | 130 / 200 |

An unverified search-engine snippet attributes to a social post the figure of "401 registered researchers and 295 teams within 4 days" of opening **[Low]**.

### 1.3 Top scores per board and task (both tracks) [High]

| Board | #1 overall | Best Human | Best Agent | 10th (both tracks) | Human top-10 / top-50 |
|---|---|---|---|---|---|
| T1:val | 76.00 EmbryoForge (H) | 76.00 | 71.81 workbuddy-science | 71.81 | 71.65 / 61.73 |
| T2 heart val_interp (E8.5) | 75.39 Omazili Manbo (H) | 75.39 | 72.32 Ablation | 71.15 | 70.75 / 67.57 |
| T2 heart val_extrap (E10.5) | 67.75 Ablation (A) | 67.31 Omazili Manbo | 67.75 | 64.95 | 63.35 / 56.82 |
| T2 embryo val_interp (E7.5) | 77.68 Omazili Manbo (H) | 77.68 | 77.26 Ablation | 75.87 | 73.97 / 68.03 |
| T3 gata4 | 80.51 YB-Salvectors (A) | 80.27 Bat | 80.51 | 77.51 | 76.96 / 68.11 |

Both-track 50th place at 15:41 UTC: T1 64.29, heart-interp 68.45, heart-extrap 58.76, embryo 69.52, T3 70.31. The 200th listed row: 50.89 / 58.76 / 50.03 / 59.03 / 47.82.

Human-track task scores (team best, task = mean over its boards):

| Task | Top-1 | Top-3 | Top-6 | Top-10 | Top-20 | Median |
|---|---|---|---|---|---|---|
| T1 | 76.00 | 74.29 | 73.70 | 71.65 | 70.55 | 51.92 |
| T2 | 73.46 | 71.37 | 69.27 | 68.87 | 67.61 | 58.95 |
| T3 | 80.27 | 79.90 | 78.17 | 76.96 | 73.00 | 57.43 |

The sum of the three Human per-task maxima is 229.73. No single team holds all three.

### 1.4 Human-track overall standings (prize slots are 1st, 2× 2nd and 3× 3rd, so top 6) [High]

| Rank | Team (affiliation as listed) | T1 | T2 | T3 | Total | heart-i / heart-e / embryo |
|---|---|---|---|---|---|---|
| 1 | EmbryoForge (Guangzhou Medical Univ.) | 76.00 | 71.52 | 77.22 | **224.74** | 73.7 / 66.4 / 74.5 |
| 2 | Ivan & Arnav (Independent) | 74.29 | 71.37 | 78.87 | 224.53 | 73.7 / 66.5 / 73.8 |
| 3 | Independent | 73.94 | 68.79 | 79.90 | 222.63 | 71.4 / 63.0 / 71.9 |
| 4 | Partytime (CUHK) | 73.70 | 68.82 | 79.41 | 221.93 | 71.3 / 63.3 / 71.9 |
| 5 | Omazili Manbo (Heidelberg) | 70.63 | 73.46 | 77.51 | 221.59 | 75.4 / 67.3 / 77.7 |
| 6 | Boom (Univ. Macau) | 74.05 | 66.87 | 79.90 | **220.82** | 70.7 / 60.2 / 69.6 |
| 7 | Bat (Peking Univ.) | 71.63 | 67.28 | 80.27 | 219.18 | 70.7 / 57.7 / 73.4 |
| 8 | Rapunzel (ANU) | 70.41 | 69.27 | 78.17 | 217.85 | 69.2 / 62.5 / 76.0 |
| 9 | Alpha (independent) | 72.60 | 69.47 | 75.16 | 217.23 | 71.2 / 66.1 / 71.2 |
| 10 | stBIOT (UCSF) | 70.59 | 69.88 | 76.13 | 216.61 | 72.0 / 61.8 / 75.9 |

- Lower thresholds: rank 20 = 207.35; rank 50 = 185.68; rank 100 = 165.85. The Human median among teams with all three tasks is 170.2.
- Agent track top: syl 220.66, Ablation 220.13, YB-Salvectors 219.09. Agents are roughly level with Human ranks 6 to 8.
- The `baseline` object in `/overall` (the "beat the best published reference" line, ceiling excluded) is 50 / 50 / 50 / 55.22 / 50, for a total of 151.74.

### 1.5 Per-metric anatomy, Human top 10 per board [High]

Per-metric values are skills on the 0 to 100 scale. Abbreviations: DE = de_score, dir = de_direction, vario = variogram, nbr = neighborhood_mmd, sev = |severity_slope|.

#### T1:val (human track; fetched 2026-10-18T15:43:48+00:00)

| # | team | headline | DE | dir | mmd | vario |
|---|---|---|---|---|---|---|
| 1 | EmbryoForge | 76.00 | 63.8 | 74.1 | 84.7 | 80.5 |
| 2 | DPath | 74.83 | 62.1 | 73.0 | 82.1 | 82.1 |
| 3 | Ivan & Arnav | 74.29 | 59.2 | 78.5 | 78.2 | 82.1 |
| 4 | Boom | 74.05 | 60.9 | 70.9 | 83.1 | 80.8 |
| 5 | Independent | 73.94 | 60.8 | 71.6 | 82.8 | 80.1 |
| 6 | Partytime | 73.70 | 61.6 | 71.5 | 81.9 | 79.3 |
| 7 | Alpha | 72.60 | 58.1 | 76.7 | 74.9 | 82.0 |
| 8 | luka | 72.15 | 59.9 | 71.9 | 78.3 | 78.4 |
| 9 | CellVR | 71.82 | 60.6 | 70.6 | 77.6 | 78.7 |
| 10 | Druggable | 71.65 | 59.2 | 72.0 | 76.2 | 80.0 |
| ref | ceiling | 100.00 | 100.0 | 100.0 | 100.0 | 100.0 |
| ref | copy_last | 50.00 | 50.0 | 50.0 | 50.0 | 50.0 |
| ref | pseudobulk_shift | 48.58 | 46.2 | 53.9 | 52.5 | 38.9 |
| stat | human top-10 median | - | 60.7 | 72.0 | 80.1 | 80.3 |
| stat | best of listed human rows | - | 63.8 | 78.5 | 84.7 | 83.3 |

#### T2:heart:val_interp (human; fetched 15:43:50 UTC)

| # | team | headline | DE | dir | mmd | vario | d2 | dice | scale | nbr |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Omazili Manbo | 75.39 | 64.9 | 71.8 | 79.9 | 62.8 | 94.7 | 47.8 | 100.0 | 79.3 |
| 2 | Ivan & Arnav | 73.73 | 63.0 | 71.3 | 72.2 | 61.2 | 93.3 | 62.7 | 100.0 | 74.6 |
| 3 | EmbryoForge | 73.67 | 67.7 | 69.2 | 74.4 | 59.3 | 97.8 | 50.9 | 100.0 | 75.0 |
| 4 | Druggable | 72.81 | 63.6 | 70.1 | 68.4 | 59.0 | 92.8 | 71.0 | 100.0 | 71.8 |
| 5 | stBIOT | 72.02 | 61.2 | 70.2 | 72.1 | 58.9 | 85.8 | 60.0 | 100.0 | 73.7 |
| 6 | Independent | 71.43 | 63.6 | 68.2 | 70.1 | 56.6 | 92.4 | 57.0 | 100.0 | 72.0 |
| 7 | SLYZ | 71.40 | 68.5 | 71.8 | 75.6 | 60.2 | 55.9 | 62.4 | 83.2 | 78.9 |
| 8 | Partytime | 71.31 | 63.0 | 68.7 | 69.6 | 57.2 | 92.4 | 57.0 | 100.0 | 71.6 |
| 9 | Alpha | 71.15 | 61.8 | 68.9 | 66.5 | 58.2 | 94.8 | 61.1 | 100.0 | 70.8 |
| 10 | VECOmics | 70.75 | 63.0 | 68.5 | 68.5 | 62.5 | 91.1 | 46.8 | 100.0 | 71.8 |
| ref | copy_last | 50.00 | 50 | 50 | 50 | 50 | 50 | 50 | 50 | 50 |
| stat | top-10 median | - | 63.3 | 69.7 | 71.1 | 59.2 | 92.6 | 58.5 | 100.0 | 72.8 |
| stat | best listed | - | 68.5 | 71.8 | 79.9 | 62.8 | 100.0 | 71.0 | 100.0 | 79.3 |

#### T2:heart:val_extrap (human; fetched 15:43:53 UTC)

| # | team | headline | DE | dir | mmd | vario | d2 | dice | scale | nbr |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Omazili Manbo | 67.31 | 59.1 | 60.3 | 63.6 | 58.9 | 100.0 | 58.5 | 100.0 | 61.6 |
| 2 | John Roos | 67.01 | 59.1 | 61.3 | 68.1 | 59.2 | 93.1 | 54.3 | 100.0 | 60.8 |
| 3 | Ivan & Arnav | 66.54 | 59.1 | 60.6 | 66.4 | 58.6 | 91.2 | 54.3 | 100.0 | 61.2 |
| 4 | EmbryoForge | 66.43 | 59.6 | 58.0 | 67.7 | 57.8 | 73.4 | 54.6 | 100.0 | 67.2 |
| 5 | Alpha | 66.05 | 58.6 | 61.0 | 64.7 | 57.9 | 87.9 | 51.2 | 100.0 | 62.8 |
| 6 | MorphoFlux111 | 65.02 | 57.5 | 56.8 | 58.9 | 59.0 | 100.0 | 59.3 | 100.0 | 57.5 |
| 7 | Eternal Existence | 63.36 | 57.5 | 61.5 | 67.1 | 62.1 | 50.8 | 51.6 | 91.2 | 64.3 |
| 8 to 10 | 哈吉米拯救世界 / 11111 / hti (identical rows) | 63.35 | 61.9 | 61.2 | 63.1 | 62.4 | 45.8 | 49.1 | 95.0 | 65.7 |
| ref | copy_last | 50.00 | 50 | 50 | 50 | 50 | 50 | 50 | 50 | 50 |
| ref | pseudobulk_shift | 47.76 | 40.9 | 49.5 | 50.1 | 44.7 | 50.0 | 50.0 | 50.0 | 48.9 |
| stat | top-10 median | - | 59.1 | 61.1 | 64.2 | 59.1 | 80.7 | 53.0 | 100.0 | 63.6 |
| stat | best listed | - | 61.9 | 65.6 | 68.1 | 62.4 | 100.0 | 60.2 | 100.0 | 67.2 |

#### T2:embryo:val_interp (human; fetched 15:43:55 UTC)

| # | team | headline | DE | dir | mmd | vario | d2 | dice | scale | nbr |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Omazili Manbo | 77.68 | 58.9 | 71.8 | 80.2 | 67.6 | 84.5 | 100.0 | 100.0 | 75.4 |
| 2 | GISL | 76.10 | 59.5 | 68.6 | 74.5 | 60.1 | 81.6 | 100.0 | 99.7 | 77.8 |
| 3 | Rapunzel | 76.02 | 60.9 | 71.3 | 77.8 | 64.8 | 72.9 | 100.0 | 99.3 | 74.7 |
| 4 | stBIOT | 75.87 | 57.0 | 70.8 | 79.8 | 63.0 | 80.8 | 100.0 | 91.9 | 75.6 |
| 5 | 关注嘉然顿顿解馋 | 75.38 | 55.8 | 71.5 | 72.0 | 68.7 | 81.7 | 100.0 | 100.0 | 73.3 |
| 6 | khj1222 | 74.82 | 53.0 | 69.5 | 74.9 | 64.7 | 81.3 | 98.9 | 96.8 | 74.9 |
| 7 | Druggable | 74.62 | 57.6 | 71.7 | 71.0 | 65.2 | 80.2 | 100.0 | 98.7 | 72.2 |
| 8 | EmbryoForge | 74.45 | 60.9 | 69.8 | 81.1 | 65.9 | 62.1 | 79.2 | 100.0 | 77.0 |
| 9 | SP | 74.19 | 57.0 | 71.0 | 67.7 | 65.5 | 90.9 | 100.0 | 100.0 | 68.9 |
| 10 | 哈吉米拯救世界 | 73.97 | 57.6 | 70.0 | 77.0 | 66.4 | 74.8 | 89.6 | 99.4 | 71.4 |
| ref | pseudobulk_shift | 55.22 | 48.6 | 50.0 | 61.5 | 57.6 | 50.0 | 50.0 | 50.0 | 61.7 |
| ref | copy_last | 50.00 | 50 | 50 | 50 | 50 | 50 | 50 | 50 | 50 |
| stat | top-10 median | - | 57.6 | 70.9 | 76.0 | 65.4 | 81.1 | 100.0 | 99.6 | 74.8 |
| stat | best listed | - | 62.4 | 73.3 | 81.1 | 69.8 | 90.9 | 100.0 | 100.0 | 77.8 |

#### T3:gata4 (human; fetched 15:43:58 UTC)

| # | team | headline | DE | dir | sev | mmd | vario |
|---|---|---|---|---|---|---|---|
| 1 | Bat | 80.27 | 74.1 | 72.7 | 99.5 | 71.0 | 80.8 |
| 2 | Boom | 79.90 | 71.4 | 73.7 | 98.3 | 74.6 | 81.7 |
| 3 | Independent | 79.90 | 71.4 | 73.7 | 98.3 | 74.6 | 81.7 |
| 4 | Partytime | 79.41 | 70.2 | 74.0 | 98.5 | 73.0 | 81.1 |
| 5 | Ivan & Arnav | 78.87 | 67.8 | 73.9 | 98.5 | 75.8 | 79.0 |
| 6 | Rapunzel | 78.17 | 69.0 | 68.7 | 99.9 | 75.6 | 78.2 |
| 7 | Omazili Manbo | 77.51 | 67.8 | 70.5 | 97.1 | 77.1 | 75.1 |
| 8 | VCDesign | 77.40 | 67.8 | 72.1 | 98.5 | 69.0 | 76.5 |
| 9 | EmbryoForge | 77.22 | 69.0 | 68.0 | 97.7 | 77.1 | 73.0 |
| 10 | SP | 76.96 | 62.5 | 70.3 | 99.4 | 75.3 | 84.5 |
| ref | wt_identity | 50.00 | 50 | 50 | 50 | 50 | 50 |
| ref | gene_ko | 42.50 | 40.4 | 48.7 | 50.0 | 37.2 | 15.5 |
| ref | shift_transfer | 42.48 | 40.4 | 48.7 | 50.0 | 37.1 | 15.5 |
| stat | top-10 median | - | 69.0 | 72.4 | 98.5 | 74.9 | 79.9 |
| stat | best listed | - | 74.1 | 74.0 | 100.0 | 80.4 | 84.5 |

**Implied raw metric values at the top [High arithmetic, Medium interpretation].** These invert the skill map `skill = d_floor/(d_floor + d)` with the public `index.json` anchors, and give local scoring targets to aim for:

| Board | Metric | Implied raw value at the top | Anchors / note |
|---|---|---|---|
| T1 | de_score | ≈ **0.37** | ceiling 0.846 |
| T1 | de_direction | ≈ 0.51 to 0.57 | ceiling 0.790 |
| T1 | mmd_u | ≈ **0.018** | floor 0.084 / ceiling 0.004 |
| T1 | variogram | ≈ 0.0012 to 0.0014 | floor 0.0052 / ceiling 0.00016 |
| Heart-interp | de_score | ≈ 0.38 to 0.44 | |
| Heart-interp | mmd_u | ≈ 0.0059 | |
| Heart-interp | nbr | ≈ 0.018 | |
| Heart-interp | occupancy_dice | 0.66 to 0.80 | ceiling 0.88 |
| Heart-extrap | de_score | ≈ 0.29 to 0.36 | |
| Heart-extrap | mmd_u | ≈ 0.012 to 0.014 | |
| Heart-extrap | nbr | ≈ 0.035 to 0.045 | |
| Embryo | de_score | ≈ 0.25 to 0.33 | |
| Embryo | mmd_u | ≈ 0.022 | |
| Embryo | nbr | ≈ 0.068 to 0.076 | |
| Embryo | occupancy_dice | at or above ceiling | anchors span only 0.705 to 0.766, so the metric is easy |
| T3 | de_score | ≈ **0.57** | ceiling 0.870 |
| T3 | de_direction | ≈ 0.58 to 0.60 | |
| T3 | \|severity_slope\| | ≈ 0.04 | fitted response slope within about 4% of 1 |
| T3 | mmd_u | ≈ 0.008 to 0.013 | |

**Where the points are.**
- **T1:** `de_score` is the weakest group everywhere. 94 of the top-200 rows are below the floor on it, and the best is 63.8. Distribution metrics are already 80 to 85.
- **T2:** `neighborhood_mmd` (25% weight, top about 70 to 79) is the largest open block. `scale_log_ratio` is saturated. `occupancy_dice` on the heart boards (top about 60 to 71) is the hardest shape term.
- **T3:** `severity_slope` is essentially solved by the leaders (≥ 97). The spread is in `de_score` (top 74) and `de_direction` (top 74). 189 of 200 listed T3 rows have `de_direction` > 50, so the field gets the Gata4 direction right even though the organisers' references do not.

### 1.6 Reference rows on the board vs what is described [High]
- **Present in `reference` (API) and in the page's tooltip dictionary:** `ceiling`, `copy_last`, `wt_identity`, `pseudobulk_shift` (T1, T2 heart-extrap, T2 embryo; absent on heart-interp, which has only one earlier stage), `shift_transfer`, `gene_ko`. Numbers are in the tables above.
- **Absent from the public API and the page code:**
  - `shift_ode`, `neural_ode`, `spatiotemporal_ode`, `dynode_flow`, `perturb_ode`. They appear only as one-line descriptions on the task pages.
  - All 8 `ctrl_*` controls: `ctrl_scale_ref`, `ctrl_shrink_ref`, `ctrl_one_cell`, `ctrl_random_cube`, `ctrl_squashed_ref`, `ctrl_scale_wt`, `ctrl_shrink_wt`, `ctrl_random_dir`. The Reference-rows page says these are "scored on every run and left on the board", but the public endpoints do not return them. They may be admin-only. This could not be verified.
- The page also supports a "Partner reference" block (EigenLab teams, not prize-eligible). It is empty now.
- "Reference baseline scores are regenerated on every scoring run, so they move whenever the panel does" (leaderboard page).

### 1.7 Headline composition and the structure gate (re-derived) [High]
- **T1:** headline = .25·DE + .25·dir + .30·mmd + .20·vario (max error 1e-4 over 200 rows).
- **T3:** headline = .30·DE + .25·dir + .25·sev + .20·(.6·mmd + .4·vario) (max error 1e-4).
- **T2:** headline = .25·(DE+dir)/2 + .25·(.6·mmd + .4·vario) + .25·**shape·min(1, nbr/50)** + .25·nbr, where shape = (d2 + dice + scale)/3 on the 0 to 100 scale. This matches all 600 T2 rows to ≤ 1e-4. The ungated formula is off by up to 0.76 points on 16 heart-extrap rows, all with nbr < 50.
- It also fits both archived community snapshots (2026-10-15 and 2026-10-18) to ≤ 1e-4 on all 1,200 T2 rows. Ungated, the 10-15 snapshot is off by up to 1.01 points on 20 heart-extrap rows and 0.12 on 1 heart-interp row.
- Two community dashboards could not explain these residuals: i-habib's "Leaderboard Anatomy" says "no cause … found", and WaelYTz guessed `count_log_ratio`. This formula closes them.
- **Implication:** keep `neighborhood_mmd` at or above its floor, and the full shape credit is available.

### 1.8 Trends [High for the numbers, Medium for the interpretation]
- 2026-10-15 10:28 UTC to 2026-10-18 11:07 UTC (WaelYTz snapshots):

  | Board | Best | 10th |
  |---|---|---|
  | T1 | 75.95 → 76.00 | 70.72 → 71.81 |
  | Heart-interp | 74.43 → 75.30 | 70.18 → 70.75 |
  | Heart-extrap | 67.67 → 67.75 | 62.69 → 64.95 |
  | Embryo | 77.66 → 77.68 | 75.82 → 75.87 |
  | T3 | 80.51 → 80.51 | 76.13 → 77.34 |

- The tubager repo records the board #1 on 2026-09-20 evening as heart-interp 73.5, heart-extrap 63.0, embryo 72.1. Since then the top moved +1.9, +4.8 and +5.6 respectively.
- Several of today's top rows were scored in September: EmbryoForge heart-interp 73.67 on 09-22; stBIOT embryo 75.87 on 09-09; GISL embryo 76.10 on 09-22.
- The board reflects each team's best of many submissions, at 8 per day per task. That is a lot of selection on validation noise. One methods journal measured noise across 6 scorer seeds on the same file: `mmd_u` ±5.75 skill points, `variogram` ±2.2, while the DE, direction and severity metrics are deterministic. Expect some regression to the mean on the hidden test.

### 1.9 Signals about competitors' methods [Medium]
- **Model labels** in the per-board rows are free text, often informative:
  - External data: "t1mean_Feng2022GSE193346_E95E145…" (APC_TDC, T1 65.9); "LongWz1018 [Feng2022 WT E9.5+E14.5; MOSTA E14.5]" (T1 61.1); "t2hx_extsame_…Feng2022-GSE193346" (heart-extrap 56.9). This shows teams use **E14.5 references**, which the rules allow after E13.5.
  - GRN or foundation-model T3: "CellVR R4 STRING11.5-GATv2/PyG-MIT IGVFtf ESM…UniProtGOcc" (T3 73.0).
  - OT / flow: "moscot" (EmbryoFlow: 62.8 / 55.2 / 64.3 on the T2 boards); "PopTransport-v4" (T1 69.2); "Source-bracket spatial state field"; "scNODE-NicheFlow-POT"; "full-native-scvi-scnode" (59.5 heart-interp).
  - LLM agent names on Agent rows: Grok 4.7, Codex, claude-opus-5, kimi-k3.
  - The top-5 Human teams mostly leave the label blank.
- **Identical metric vectors across differently named teams** suggest shared pipelines:
  - Boom = Independent on T3 (both "Carina"). Partytime = VCDesign on heart-extrap (both "Nunki"). These four teams all use star-name labels (Pictor, Tucana, Wezen, Menkar, Bellatrix, Canopus, Kochab, Rigel, Diphda, Shaula, Elnath, Dorado).
  - "local_algorithm_treatment" is shared by 4 teams (11111 / hti / in / 哈吉米拯救世界).
  - "agent_sota_v5_*" is shared by IAI Lab (HIT) and linlab.
  - One row says "authorized cross-team reuse: jiaoyifeng/i112" (team REVISE).
  - Rules prohibit code sharing across teams. Whether organisers will act is unknown. If shared clusters are merged or disqualified, the prize line moves.

---

## 2. Baseline Colab notebook (Drive file `11kHZjd4IawoC6uwCqiHNuDGAYpoPlUQO`) [High]
The notebook was downloaded via `drive.usercontent.google.com` (public, about 28 KB) and is saved at `SCRATCH/landscape_work/colab/nb_download.bin`. Title: *"The shift baseline — a real submission for T1, T2 and T3 in ~10 lines each."* Contents:
- **Setup:** `pip install -q veckit`. Note that this currently installs **0.1.2**, the PR #6 metric variant (see §6). The setup then downloads the six 150-cell sample h5ads from `aristoteleo/veckit/data`.
- **`pseudobulk_shift(prev, last, damp=1.0)`** (T1/T2):
  - Δ_k = mean_k(last) − mean_k(prev) per cell type; x̂ = max(0, x + damp·Δ_{type}).
  - Types absent in prev get Δ = 0.
  - Builds from `last.copy()`, so `obsm['spatial_3D']` is carried over unchanged for T2.
- **`shift_transfer(train_wt, train_ko, target_wt, ko_gene=None, damp)`** (T3):
  - One global Δ = pb(KO_train) − pb(WT_train), added to the target WT and clipped at 0.
  - With `ko_gene`, that gene is zeroed, which is `gene_ko`.
  - **For the test KO, Ctnnb1 is not on the 500-gene panel, so `gene_ko` reduces to `shift_transfer`** (checked against the panel file).
- `copy_last` is shown as the floor. No ODE or flow models are in the notebook.
- **Published numbers in the notebook are not challenge scores.** They are raw veckit metrics on the 150-cell samples with `--target` pointed at the input stage itself, as the notebook says explicitly:

  | Run | Values |
  |---|---|
  | T1 shift | de_score 0.5455, de_direction 0.8618, mmd_u −0.00443, variogram 0.00201 |
  | T1 copy_last | de_score 0.5152, de_direction 1.0, variogram 0 |
  | T2 shift | de_score 0.510, de_direction 0.993, d2_shape 0.00126, occupancy_dice 1.0 |
  | T3 shift_transfer (same pair) | de_score 0.641, de_direction 0.998, severity_slope −0.094 |

- **The organisers' suggested "knobs", in order:**
  1. `damp = Δt_target/Δt_observed`
  2. finer cell types
  3. predict geometry (growth and deformation)
  4. per-gene, GRN-conditioned KO responses
  5. per-cell velocity, which is "where the ODE baselines (`neural_ode`, `dynode_flow`) begin"
- The real-data scores of these baselines are the live-board reference rows in §1.6.

---

## 3. Organisers' reference models: what they are and what they score

| Reference | Task | Definition (site wording) | Live-board score |
|---|---|---|---|
| copy_last | T1/T2 | last observed stage verbatim, with its coordinates | 50 (floor) |
| wt_identity | T3 | matched WT verbatim | 50 (floor) |
| pseudobulk_shift | T1/T2 | per-type constant velocity, real cells, clip 0, coordinates unchanged | T1 **48.58**; heart-extrap **47.76**; embryo **55.22**; not on heart-interp |
| shift_ode | T1 | "kept because its minimiser is that same constant shift" | not shown |
| neural_ode | T1 | "autonomous field in gene space, flow-matched on a global entropic-OT coupling"; "does not beat the one-line constant shift" | not shown (implied < 48.6) |
| spatiotemporal_ode | T2 | "the reference dynamics model … score comes almost entirely from an analytic growth-rate extrapolation rather than from the ODE" | not shown |
| dynode_flow | T1, T2, T3 | "the external flow-matching reference model, unmodified"; on T3 "reuses the Task 2 checkpoint through its perturbation path, with no separate training" | not shown |
| shift_transfer | T3 | training KO's pseudobulk Δ added to target WT | **42.48** (dir 48.7 means wrong direction; vario 15.5) |
| gene_ko | T3 | shift_transfer plus zeroing the KO gene (task page: "adds GRN conditioning on the perturbed gene") | **42.50** |
| perturb_ode | T3 | "an integrated perturbation field rather than a fixed transferred shift" | not shown |

All [High], from the site text, the JS bundle and the API.

Reading: on these validation boards the organisers' **published** references are at or below the floor, except `pseudobulk_shift` on the embryo board. The proposal also promised "pseudotime and OT-coupled flow-matching extrapolation" (T1), "optimal-transport alignment and integration" (T2) and "conditional generative baselines" (T3). None of these is on the public board [High, proposal §1.6].

**dynode_flow and Navigo.**
- For the Navigo reading:
  - the description fits a population flow-matching model with a zero-shot perturbation mode;
  - Navigo is from the organiser's lab, with co-authors Xiaojie Qiu and Zehua Zeng (the veckit maintainer);
  - the Colab names `dynode_flow` together with `neural_ode` as per-cell-velocity models.
- Against it:
  - "dynode" is not "Navigo";
  - Navigo's code requires `layers['Ms']`, `layers['Mu']` and `obs['time']`, and the challenge's MERFISH T2/T3 data has no spliced/unspliced split. Running it "unmodified" on T2 is therefore not possible as published, unless the organisers set u = s or use raw-read layers that participants do not get.
- Verdict: **plausibly Navigo or a close relative; unverified [Medium/Low]. No public scores exist for it.**

---

## 4. Navigo, Spateo and dynamo: relevant capabilities

### 4.1 Navigo [High unless noted]
- **What it is:** "Generative Modeling of Mouse Embryogenesis for Fate and Disease Prediction", Fan Y., Liu X., Wang Y., Zeng Z., Li L., **Qiu X.**, **Li Y.** (corresponding author, CUHK CSE). bioRxiv, posted 2026-06-24, DOI 10.64898/2026.06.18.733286, CC BY-NC-ND. Code: `aristoteleo/Navigo-release` (BSD-2-Clause, created 2026-03-19, last push 2026-06-10, 19 stars). `aristoteleo/Navigo-docs` is the Sphinx/ReadTheDocs mirror of `docs/`, synced by a GitHub Action (last push 2026-04-07). It contains no extra results. The `navigo` package is **not on PyPI**; install from git.
- **Model:**
  - State x = (s, u) ∈ R^{2M} (spliced and unspliced moments); dx/dt = v(x, t).
  - Network `MLPTimeGRN`: fc1(s) + fc1_t(s)·t → ReLU → fc2 → ReLU → fc3 → `alpha_predictor` (ReLU). β and γ are per-gene sigmoid parameters.
  - v_u = α(s,t) − β∘u; v_s = β∘u − γ∘s. Transcription α is cell- and gene-specific; splicing and degradation are gene-specific.
- **Training:**
  - Population flow matching with **iterative reflow**. Adjacent time points start with random pseudo-pairs.
  - Linear-path flow-matching loss (scaled by Δt).
  - Forward-integrate, then re-pair cells by nearest neighbour, then refit.
  - Metacells are aggregated within embryo × cell type × time before training.
  - `submission/main_navigo.py` defaults: hidden 5012×5012, 10 rounds × 200 steps, batch 25, lr 4e-5, 10 ODE steps, per-feature min-max scaling. Input `.h5ad` must have `obs['time']`, `layers['Ms']`, `layers['Mu']`.
- **Inference:**
  - Forward simulation (Euler integration, clipping at 0).
  - **Aligned interpolation:** forward-simulate from t_i to t_{i+k}, match to observed cells by nearest neighbour, then blend α·x^(i) + (1−α)·x_match with α = (k−m)/k.
  - Denoising.
- **Perturbation:** set gene g ∈ G to m × its cell-type mean (m ≤ 0 means a knockout), clamp it during rollout, and compare with the wild-type rollout. GRN analysis uses the Jacobian ∂α_i/∂s_j.
- **Reported results:**
  - Atlas: Qiu C. et al. 2024 *Nature* (12.4 M nuclei, 43 time points).
  - "Zero-shot" KO predictions are validated on **six knockout genotypes**. The public tutorial compares them qualitatively (Gli1/2/3 and Scn-family lineage effects versus experimental DEG burden), with no DE-overlap numbers.
  - The interpolation tutorial uses an 84,696-cell, 21-time-point subset (E8.75 to E18.75 in half-day steps), 3,902 HVGs (input dimension 7,804), every other stage held out.
  - EMD of the Navigo prediction vs copying the source stage: better in **8/10** held-out stages. Vs copying the next stage: better in 8/10. Better than both: 6/10.
  - Mean EMD 6.03 vs 6.68 (source) and 6.57 (next). Per-stage DEG Jaccard overlaps are only 0.07 to 0.30.
  - **Takeaway:** even in dense, half-day **interpolation**, the organisers' own flow model beats "copy the nearest stage" by only about 10% in EMD. That matches the observation that copy_last is a hard floor.
- **Relevance to us:**
  - The aligned-interpolation recipe (simulate, match to real cells, blend) is directly usable for the T2 embryo E7.75 interpolation.
  - The clamp-and-rollout perturbation idea is relevant for T3, but **Ctnnb1 is not on the 500-gene panel**. Any clamp-the-gene approach has to run in a gene space that contains Ctnnb1 (for example the T1 32,285-gene scRNA, where Ctnnb1 and [redacted 2026-10-19: B] are all present) and then be mapped onto panel genes.
  - The spliced/unspliced requirement is a blocker for running it unmodified on challenge files.

### 4.2 Spateo (spateo-release 1.1.1 on PyPI, 2025-03-25; BSD-2) [High]
- Qiu X. et al., "Spatiotemporal modeling of molecular holograms", *Cell* 187:7351–7373 (2024), doi:10.1016/j.cell.2024.10.011 (erratum *Cell* 188:1744, 2025). Built on an 8-million-cell Stereo-seq E9.5/E11.5 embryo atlas. Co-author Qingquan Zhang (UCSD) also generated the challenge data.
- Capabilities relevant to T2:
  - `st.align.morpho_align` / `paste_align`: partial, non-rigid alignment of slices and stages, with `BA_transform`, TPS and grid deformation, mesh correction.
  - `st.tdr` 3D models: point cloud, mesh, voxel; `model_morphology` (volume, area, density); `pairwise_shape_similarity`.
  - **Morphometric vector fields**: `morphofield_sparsevfc`, `morphofield_gp`, `cell_directions`, `morphopath`.
  - Differential geometry: velocity, acceleration, curvature, curl, divergence, torsion, Jacobian.
- These give a principled way to estimate per-cell displacement and growth between released stages (alignment, then a vector field). The extrapolated field can then predict E12.5 or E7.75 geometry.
- **Caveat:** community teams report that every attempt to deform geometry for **extrapolation** lowered `d2_shape` and `occupancy_dice` on the board; their best heart-extrap entries only rescale. See §8.

### 4.3 dynamo (dynamo-release 1.5.3, 2026-03-17; BSD-3) [High]
- Qiu X. et al., *Cell* 185:690–711 (2022), doi:10.1016/j.cell.2021.12.045.
- Capabilities: vector-field reconstruction (SparseVFC) from velocities; analytic Jacobian; `dyn.pd.perturbation` (δY = J·δX propagated through the learned field); `dyn.pd.KO`; least-action paths; fate prediction.
- **Requirements:** spliced/unspliced or metabolic-labelling velocities, which the challenge files lack. A community team flagged this as a "data-contract" blocker. However, the vector-field and Jacobian machinery accepts **any** velocity estimate, for example OT-displacement vectors between stages in PCA space.
- This is a plausible route to a Jacobian-based T3 response model trained on WT development. It is untested here.

---

## 5. Nature Methods comment: what the organisers value [High for content, Medium for implications]
Cao N., Lu Y., Qiu X. "Towards predictive virtual embryos with genomics and AI." *Nat Methods* (2026), doi:10.1038/s41592-026-03055-4. Text read from `site_assets/natmeth.txt`. Key positions:
- **Predictive, not descriptive.** The aim is a "data-backed digital twin" simulating embryogenesis under normal and genetically or environmentally perturbed conditions. Atlases (TOME, MOSTA, the 3D reconstructions) are "largely statistical", unlike predictive biophysical simulation.
- **Named method families:**
  - transformer foundation models (Geneformer, scGPT), whose embeddings are "transferred or fine-tuned within Virtual Embryo frameworks";
  - **flow matching inspired by optimal transport**, which "learn vector fields that link cell states over developmental time, enabling interpolation of intermediate states and estimation of continuous migration paths";
  - **Spateo morphometric vector fields** linking tissue deformation (curvature, divergence, torsion) to morphogen gradients and transcriptional programs;
  - multiscale embeddings (gene → cell → niche → tissue → organ) with cross-attention;
  - coarse-graining ("modeling cells at the tissue or coarse level") to manage complexity.
- **Challenges they emphasise:**
  - spatial registration to a standard embryo reference, and temporal alignment;
  - scarce time-resolved 3D data;
  - scaling;
  - **interpretability and agreement with known biology**, with the example "perturbing Tbx5 in a particular time frame causes heart malformation";
  - scarcity of systematic in vivo perturbation data (Mouse Perturbation Atlas).
  - They propose **annual Virtual Embryo Challenges**, a CASP analogue, "to promote community standards and rigorous benchmarking".
- **Implications:**
  - Ranking is purely metric-based. The method report (due within 14 days for winners), the Generality Award (one architecture across at least 2 tasks, $8K) and any post-competition paper will be read against this vision.
  - A unified spatiotemporal generative model would align well: flow or OT dynamics on expression, a vector field on geometry, and a perturbation path for T3.
  - Expect the organisers to value **mechanism-consistent T3 responses** (Wnt/β-catenin biology) and explicit handling of cell populations appearing.

---

## 6. veckit (official local scorer): all issues and PRs, versions, behaviour changes [High]
Repository `aristoteleo/veckit` (MIT): created 2026-08-10, 16 stars, 7 forks. Commits on `main`:
- `16d905e`, `27ae299` (0.1.1), `46d41e6` on 2026-08-10 (Zehua Zeng)
- `a905603` / `9b27c34` / merge `cc3859a` on 2026-10-02 (PR #9)

No tags or releases on GitHub. PyPI has 0.1.0, 0.1.1 (both 2026-08-10) and **0.1.2 (2026-08-25)**. Branches: `main`, `feat/metric-hardening-0.2.0` (PR #6), `docs/library-size-ratio-no-renormalisation` (PR #12). The `aristoteleo/vec_baselines` repo referenced in PRs #6 and #12 is **not publicly accessible** (HTTP 404, presumably private). It is likely the internal tree that the Reference-rows page cites (`bake.py`, `T1/baselines.py`, `T2/model.py`, `T3/model.py`); this is unverified.

| # | Type / state | Date | Gist | Takeaway for us |
|---|---|---|---|---|
| 1 | Issue, open | 08-18 | Asks for the canonical scorer target. **Maintainer (Zehua Zeng):** "use `veckit==0.1.1`"; for T2 "the eight board-anchor metrics are the leaderboard metrics"; `energy_distance`, `sliced_wasserstein`, `count_log_ratio` are diagnostics (`count_log_ratio` audit-only); `--setting` is metadata-only; local scores do not predict leaderboard scores. | Rank on the 8 / 4 / 5 anchored metrics only |
| 2 | PR, open | 08-19 | Memory-safe sparse T1 scorer (peak RSS 14.5 → 6.5 GiB, identical metrics) | Local T1 scoring needs about 15 GB RAM with the stock code |
| 3 | Issue, open | 08-19 | T1 docstring says "500-gene"; maintainer: stale; T1 is the ~32k-gene 10x Multiome release; no score affected | — |
| 4 / 5 | PRs, open | 08-22 | Synthetic T2 conformance tests; executable synthetic versions of the 8 controls | Useful test fixtures |
| 6 | PR, open (organisers' branch) | 08-25 | **BREAKING:** `severity_slope` → `severity_ratio = min(β, 1/β)` (old form "compressed every realistic prediction into the top few points — a response 2× too strong still scored 91/100"). `scale_log_ratio` → `size_fidelity` = tissue scale **and** cell scale (median NN distance) ("rescaling the floor's coordinates took the old score from 50.0 to 100.0 with every other shape metric unchanged"). Version bump to 0.2.0. | **Not on the board as of today** (anchors still `severity_slope` / `scale_log_ratio`), but PyPI 0.1.2 ships it. If adopted later, coordinates need a realistic cell spacing as well as the right radius, and T3 magnitude scoring tightens |
| 7 | Issue, closed | 09-17 | PCA canonicalisation sign ambiguity: about 50% of proper rotations scored like mirrors in `occupancy_dice` and `sliced_wasserstein` (Dice 1.0 → 0.627 on the sample) | Fixed in #9 |
| 8 | PR, open | 09-24 | MMD memory/speed (identical scores) | — |
| 9 | PR, **merged 10-02** | 10-02 | Force a right-handed PCA frame (`det(Vt)<0` → flip), guarded for n<3. Reviewer note: "any precomputed floor/ceiling baselines need regenerating after this merges" | Board anchors may or may not include this; `index.json` anchors did not change after 10-02 (rules watch) |
| 10 | Issue, open | 10-03 | Backend (`kg.virtualembryo.ai`) returned 502; sign-in down | Outage precedent; rules §12 allows remedies |
| 11 | PR, open | 10-06 | Moran's I self-neighbour bug with coincident coordinates (`morans_I_agreement` is reported, not ranked) | Avoid duplicate coordinates anyway |
| 12 | PR, open (organisers' branch) | 10-15 | Docs: "the scorer never renormalises a submission"; raw counts are scored as if log-normalised; released files have `expm1(X).sum(1) == 1e4` exactly | Always submit log1p CP10k |

**Scorer variants (community audit, i-habib `scorer-versions`, 2026-10-14) [Medium-High]:**
- PyPI 0.1.2 = PR #6 scoring modules, without the #7 fix.
- `main` @ `cc3859a` = 0.1.1 plus the #7 fix only. Install with `pip install "git+https://github.com/aristoteleo/veckit@cc3859ab2a2305efc4544815d3f44f660efa77d4"`.
- T1 metrics are identical in every variant.
- Without the #7 fix, 5 of 12 random rotations of a perfect E9.5 prediction dropped `occupancy_dice` from 1.0 to as low as 0.83. That costs up to 1.7 points on heart-interp and 3.4 on heart-extrap.
- "The version the leaderboard itself runs is not published."

---

## 7. Announcements and rule or scoring changes (timeline) [High unless noted]

| Date | Change | Source |
|---|---|---|
| 2026-07-28 | NeurIPS 2026 competitions announced. Listed as "The Virtual Embryo Challenge: Generative Modeling of Embryogenesis Across Space, Scale and Time" under Healthcare, Biology and Humanitarian Applications (no description given) | NeurIPS blog |
| 2026-08-10 | Terms effective; veckit 0.1.0/0.1.1 released | terms page, PyPI |
| 2026-08-15 | P2 development phase (validation boards) | timeline |
| 2026-08-25 | veckit 0.1.2 on PyPI (PR #6 metrics; not the board's) | PyPI, PR #6 |
| from 2026-08-26 | Rules §10 external-data restatement "corrects the previous wording rather than tightening it"; the embryo interval is wider than before | rules page |
| undated (P2) | **Structure gate** added to T2 shape credit after `ctrl_random_cube` beat the floor on both interpolation boards (formula re-derived in §1.7) | Reference-rows page |
| undated | DE metrics changed to measure a null "built from exactly that attack [scaled copy of reference], in both directions" | Task 1 page |
| 2026-10-02 | veckit PR #9 merged (right-handed PCA frame) | GitHub |
| 2026-10-03 | Backend outage (kg host, 502) | veckit #10 |
| 2026-10-16 | FAQ and Data: "For the extrapolation targets (**Task 1, and Task 2 in the heart setting**) the window is stated outright" (E9.5 < t ≤ E13.5 excluded; E9.5 allowed; after E13.5 allowed with disclosure). Evaluation pages: cell minimum is "the board's own minimum" (embryo: [redacted 2026-10-19: A]), not 1,000 | community rules-watch changelog |
| 2026-10-20 | P3: validation answers released for all tasks; test inputs released; **2 official scored submissions per board** for the whole phase | timeline, rules §11 |
| 2026-12-02 / 12-04 / 12-11 | Final submissions / official evaluation / winners at NeurIPS | timeline |

Not found: X/Twitter posts by the organisers.
- @Xiaojie_Qiu is linked from the site, but no tweets were retrievable.
- A nitter mirror for @Nanguage (a veckit reviewer) returned 502.
- Announcements are said to go "to Slack and Discord first" (FAQ). Those channels need joining and were not accessed.
- The proposal's original schedule (P3 on 10-25, final on 11-02) was superseded by the current timeline.

---

## 8. Public participant resources and what they reveal [Medium]
The organisers' Community Contribution Award ($20K pool, up to $200 per contribution) has produced **15 reviewed and 357 filed** resources (public list at `GET /challenge/community/register`). Most useful for us:
- **Scorer and leaderboard tooling:**
  - i-habib/virtual-embryo-community: scorer-versions, metric-lens, leaderboard-anatomy, val-to-test, history-rescore, count-recovery (exact integer UMIs from T1 log1p files), panel-normalise, data-safety, submission-blender, ensemble-lab, morph-to-truth.
  - xxx12e/vec-community-kit (validator, local scoring wrapper pinned to `46d41e6`, **daily rules watch** at `rules-watch/CHANGES.md`).
  - AJK0921/vecbench (noise-floor calibration).
  - WaelYTz/vec-leaderboard-insights (daily snapshots).
  - gh-dv-openclaw/vec-submit-check and wqty123/vec-scalecheck (`expm1(X).sum(1)==1e4` check).
- **Note on authorship:** the #2 Human team "Ivan & Arnav" appears to be **Ivan Habib + Arnav Mana**. Their public tools emphasise ensembling ("mean graft": X_B + (μ_A − μ_B); "quantile graft"; "spatial transplant": expression from model A, coordinates from model B by Hungarian matching), validation-to-test transfer, and metric attribution. This is an inference from names and tooling, not a method disclosure.
- **wqty123/virtual-embryo-challenge-methods** (full journal of one mid-table team, with 2026-09-23 scores T1 53.7, embryo 65.3, heart-interp 65.0, heart-extrap 54.5, T3 66.6). Lessons:
  - Global shifts cap `de_score`.
  - A Markov/flux transfer of the E8.5→E9.5 rule falls below the floor (44.2).
  - Noise or anonymous program mixing destroys the variogram.
  - Missing a cell type collapses `mmd_u`.
  - Geometry changes on heart-extrap always hurt (best is "scale_only").
  - T3 global Δ transfer gives 41.0 ([redacted 2026-10-19: B]).
  - Mixing real KO cells into WT carriers gives 64.6. Resampling by cell type and diluting with real WT cells gives 66.6.
  - "VGFM / CVAE / diffusion … never stably exceed copy-last on the T2 ladder."
  - Submitting a full WT copy scored 45.8 (< 50), with DE 38.8 and dir 46.4. My interpretation, not verified: the scorer's reference is a 10% subsample, so a full-WT copy has a small noisy dp, and the DE nulls score that noise below zero. The official floor row, by contrast, has dp = 0 exactly.
  - Their "structure theorem": de/dir/sev depend only on the pseudobulk; the variogram depends on within-cell gene-pair structure; the best entries keep real cells and leave expression untouched.
- **tubager/vec-t2** (claims T1 #1 around 09-20; T2 lock 69.7 / 69.7 / 57.6 on 09-20). Lessons:
  - A panel OT-CFM "variogram explodes 4–6×", so flow is disabled.
  - "Don't lerp expression; real cells + one geometry cloud."
  - Cluster-wise Hungarian placement of real bracket cells onto an anchor cloud ("OT-place"): embryo 66.4, heart-interp 68.3.
  - Cluster partial-mean expression edits: 69.1 to 69.7.
  - Per-submission row-order noise is about ±0.6 headline from subsampling; DE and SLR are noise-free.
- **ShaoliZhao/virtual-embryo-EmbryoMatch** (stage-aware external dataset registry); i-habib external-data catalog (GSE282547 E14.5 heart Visium covers 500/500 panel genes with coordinates; GSE247450 E9.5 MERFISH; sc3D GSE197353 E9.0 Slide-seq).

---

## 9. What score levels are competitive, and implications
- **Validation prize line (Human):** about **221** total (6th = 220.82); about **224.5 to 225** for 1st. Per task: T1 about 74 to 76, T2 about 69 to 73, T3 about 78 to 80. Getting into the top 10 needs about 216.6.
- **The test will be different [Medium, inference]:**
  - **T2 test = embryo E7.75 (interpolation) plus heart E12.5 (extrapolation).** Heart-interp, today's easiest heart board, drops out. Heart-extrap is the lowest-scoring board (top 67.8), and the E10.5 → E12.5 step is 2 days rather than 1. Expect T2 test task scores to fall.
  - T1 test is E12.5 from E8.5/E9.5/E10.5. The 2-day step is longer, but **external data after E13.5 is legal** (E14.5 Feng 2022 GSE193346, MOSTA E14.5 Stereo-seq, Qiu 2024 atlas cells after E13.5) and allows bracketing.
  - T3 test (β-catenin) is described as hardest: Ctnnb1 and canonical Wnt targets are off-panel, but there are now **two training KOs** (Mab21l2 @ E9.5 and Gata4 @ E8.75 × 2 replicates; Gata4 is at the test stage).
  - Each team gets only **2 scored test submissions per board**, so validation-board hill-climbing is worth less. Community pseudo-pair analysis found validation-to-test rank agreement of Spearman about 0.88 for simple candidate families.
- **The organisers' own models do not set the bar.** Their best published row is 55.2 (embryo); everything else is ≤ 50. The bar is set by the field: about 20 to 30 points above the floor at the top.
- **Robust levers seen in public:**
  1. Keep real cells and edit the pseudobulk separately (graft) to win DE without destroying the variogram or MMD.
  2. Model cell-type composition change and the appearance of new populations; missing populations are heavily penalised.
  3. Geometry: right RMS scale, intact shape, neighbourhoods at or above floor. Avoid large deformations on extrapolation unless validated.
  4. T3: calibrate magnitude (sev ≥ 97 is table stakes) and win on DE and direction.
  5. Hedge against a PR #6 switch: keep a realistic cell-to-cell spacing (median NN distance), not only the right RMS radius, and do not overshoot the response magnitude.

---

## 10. Open questions and unverified items
- What `dynode_flow`, `neural_ode`, `spatiotemporal_ode`, `shift_ode` and `perturb_ode` score on the real boards, and whether `dynode_flow` is Navigo. Not public.
- Which veckit build the server runs: whether it includes the #7 fix, and whether PR #6 will be adopted for P3 or the final evaluation. "Datasets and exact metric formulations can still be corrected" (timeline page).
- Whether the adversarial-control rows are still scored (the page says yes; the API does not show them).
- Test-board keys, cell bounds and floors (in `index.json` on 10-20). How the two official files per board are designated. Whether P3 returns per-metric breakdowns.
- Whether the rules §10 held-out list stops applying to validation stages once their answers are released.
- Discord/Slack announcements and any X posts. Not accessed.
- Methods of the top-5 Human teams (no public write-ups found).

---

## References
- Leaderboard page: https://virtualembryo.ai/challenge/leaderboard ; API (public GET): https://api.virtualembryo.ai/challenge/boards , https://api.virtualembryo.ai/challenge/overall , https://api.virtualembryo.ai/challenge/leaderboard?task=T1&board=T1:val , https://api.virtualembryo.ai/challenge/phase , https://api.virtualembryo.ai/challenge/community/register
- Board contract: https://virtualembryo.ai/challenge/panels/index.json ; Reference rows: https://virtualembryo.ai/challenge/baselines ; Tasks: https://virtualembryo.ai/challenge/tasks/temporal , /spatial , /perturbation ; Rules: https://virtualembryo.ai/challenge/rules
- Baseline Colab: https://colab.research.google.com/drive/11kHZjd4IawoC6uwCqiHNuDGAYpoPlUQO (raw: https://drive.usercontent.google.com/download?id=11kHZjd4IawoC6uwCqiHNuDGAYpoPlUQO&export=download)
- Competition proposal PDF: https://virtualembryo.ai/challenge/papers/neurips-2026-competition-proposal.pdf
- NeurIPS 2026 competitions: https://blog.neurips.cc/2026/07/28/neurips-2026-competitions-announced/
- Cao, Lu & Qiu, Nat Methods (2026): https://doi.org/10.1038/s41592-026-03055-4 (PDF: https://virtualembryo.ai/challenge/papers/towards-predictive-virtual-embryos.pdf)
- veckit: https://github.com/aristoteleo/veckit ; issues: https://github.com/aristoteleo/veckit/issues?q=is%3Aissue ; PR #6: https://github.com/aristoteleo/veckit/pull/6 ; PR #9: https://github.com/aristoteleo/veckit/pull/9 ; PyPI: https://pypi.org/project/veckit/
- Navigo: https://github.com/aristoteleo/Navigo-release ; preprint https://doi.org/10.64898/2026.06.18.733286 (https://www.biorxiv.org/content/10.64898/2026.06.18.733286v1)
- Spateo: https://github.com/aristoteleo/spateo-release ; Qiu et al., Cell 2024 https://doi.org/10.1016/j.cell.2024.10.011 (erratum https://doi.org/10.1016/j.cell.2025.02.026)
- dynamo: https://github.com/aristoteleo/dynamo-release ; Qiu et al., Cell 2022 https://doi.org/10.1016/j.cell.2021.12.045
- Qiu C. et al., Nature 2024 (12.4 M nuclei atlas): https://doi.org/10.1038/s41586-024-07069-w
- Feng W. et al., Nat Commun 2022 (GSE193346): https://doi.org/10.1038/s41467-022-35691-7
- Community: https://github.com/i-habib/virtual-embryo-community (scorer-versions, val-to-test, leaderboard-anatomy, submission-blender, ensemble-lab, count-recovery) ; https://github.com/xxx12e/vec-community-kit (rules-watch/CHANGES.md) ; https://github.com/wqty123/virtual-embryo-challenge-methods ; https://github.com/tubager/vec-t2 ; https://github.com/WaelYTz/vec-leaderboard-insights ; https://github.com/AJK0921/vecbench ; https://github.com/ShaoliZhao/virtual-embryo-EmbryoMatch
- Local copies (this run): `SCRATCH/landscape_work/api/*.json` (leaderboard), `SCRATCH/landscape_work/colab/nb_download.bin`, `SCRATCH/landscape_work/gh/veckit_*.json`, `SCRATCH/landscape_work/navigo/*.ipynb`, scripts in `SCRATCH/landscape_scripts/`. SCRATCH = SCRATCH
