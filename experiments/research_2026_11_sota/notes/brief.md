# Brief: state-of-the-art research for team Mithril, 2026-11-06

AS_OF 2026-11-06. Repository: /mnt/md0/IITM/BackUp/Home/vijayavallabh/virtual-embryo-challenge (read-only for you).
Everything in this brief comes from the team's own clean records (`docs/CHALLENGE.md`, `METHOD_SUMMARY.md`,
`experiments/harness/TRIES.md`, `docs/PLAN.md`, `docs/DATA_POLICY.md`, `results/upload_gate.json`).

## 1. The ask

The captain: research state-of-the-art papers and reports, then synthesise solutions that maximise our chance of a
top-3 place on the Human Team leaderboard; his target is a total above 209. Earlier searches (2026-10-19 sweep,
2026-11-04 rebuild design, 2026-11-05 DGX campaign design, 2026-11-05 cross-disciplinary research) found only small
levers (about +1 to +3 Overall), mostly because changes need legal, counted evidence. This round goes after what is
NEW (2025 to Nov 2026) or deeper than before, and above all after anything that could create legal evidence or a
large gain.

## 2. The challenge (Virtual Embryo Challenge, NeurIPS 2026; Track 1, Human Team)

Data: ~1M mouse cells, 11 stages E6.75-E12.5, unpublished (Neil Chi lab, UCSD); organisers include the Xiaojie Qiu lab
(GitHub `aristoteleo`: veckit scorer, dynamo, spateo, Navigo).

| Test board | What is predicted | Reference / floor (scores 50) | Counts in Overall |
|---|---|---|---|
| T1:test | a SET of E12.5 cells, scRNA (Multiome-chemistry nuclei, whole-cell-like RNA), all 32,285 genes, log1p CP10k; heart-centred dissection ("only the heart chamber", organiser note) | copy_last = E10.5 (validation stage, unlocked for us) | 1 |
| T2:heart:test_extrap | E12.5 3D MERFISH, 500-gene panel + 3D coordinates | copy_last = E10.5 with its coordinates | 0.5 |
| T2:embryo:test_interp | E7.75 3D MERFISH (whole embryo, gastrulation), 498 genes + coords; brackets E7.5 (val, unlocked) and E8.0 | copy_last = E7.5 | 0.5 |
| T3:bcatenin | beta-catenin conditional KO (Mesp1-Cre; Ctnnb1 F/F) at E8.75, 3D MERFISH 500 genes; Ctnnb1 is NOT on the panel | wt_identity = matched WT E8.75 | 1 |

Training/validation material we may use: T1 E8.5, E9.5, E10.5; heart MERFISH E8.25, E8.5, E8.75, E9.5, E10.5; embryo
MERFISH E6.75, E7.25, E7.5, E8.0; T3 Mab21l2 KO at E9.5 (training), Gata4 KO at E8.75 (Mesp1-Cre; Gata4 F/F;
Gata6 F/+, two replicates; validation, unlocked), WT E8.75 and E9.5.

**Metrics** (veckit; pb = pseudobulk mean of log-normalised X; ref = reference stage or matched WT; dp = pb(pred) -
pb(ref); dt = pb(true) - pb(ref)):
- de_score: truth DE genes = Mann-Whitney + BH + |lfc| >= 0.25, split n_up/n_dn; the prediction contributes its top
  n_up and bottom n_dn genes by dp; signed overlap, chance-corrected by an expression-magnitude null (ranking by
  +pb(ref) or -pb(ref)). Rank-based: rescaling dp is a no-op.
- de_direction: rank partial correlation of dp vs dt over all genes, controlling for the rank of pb(ref).
- severity_slope (T3 only): log of the zero-intercept OLS slope of dp on dt over truth DE genes; floored at log(1e-3)
  if slope <= 0 or R^2 <= 0.01; scored on |value|.
- mmd_u: unbiased multi-kernel RBF MMD^2 on 2,000 subsampled cells in a 30-PC space fitted on the TRUTH; bandwidths
  0.25-4x median heuristic.
- variogram: 20,000 random gene pairs, 1,500 cells, E|x_i - x_j|^0.5 compared between prediction and truth (gene-gene
  co-variation and per-gene sparsity).
- T2 shape: d2_shape (1-D Wasserstein of within-cloud pairwise distances / RMS radius), occupancy_dice (canonicalised
  16^3 voxels), |scale_log_ratio| (RMS radius); shape credit multiplied by a neighbourhood-structure gate;
  neighborhood_mmd (mmd_u on 15-NN-averaged expression).
- Skill per metric = d_floor / (d_floor + d), floor -> 0.5, ceiling (split-half of the target) -> 1.0; NaN -> 0.
- Weights. T1: de_score .25, de_direction .25, mmd_u .30, variogram .20. T2: four groups .25 each (expression change =
  de_score + de_direction; cell state = mmd_u .6 + variogram .4; shape & growth = d2, Dice, scale 1/3 each; local
  organisation = neighborhood_mmd). T3: de_score .30, de_direction .25, |severity| .25, mutant-state .20 (mmd_u .6,
  variogram .4).
- Overall = T1 + (heart + embryo)/2 + T3. Two scored uploads per board for the whole final phase; the better counts.

## 3. Where we stand (our own committed results)

- Estimate (each board at the better of its two uploads): committed outlook 184 (170-198); the coded reading 191.5
  (T1 69.6 [an upper reading], heart 53.9, embryo 60.9, T3 64.6). First round alone about 170.
- A screened read of overall totals puts a Human top-3 place near 202-210 (overall totals, calibration only: never a threshold or evidence for a try, adoption or method choice). Our gate: 209 (it gates uploads only). [Label added 2026-11-06 after the counter-review; the brief went to the agents without it.]
- Points left to the split-half ceiling, Overall units (headroom tables): T3 response identity (de_score) 15.2;
  T1 distribution (mmd_u) 15.1; T1 DE recovery 13.0; T1 covariance (variogram) 11.5; T3 direction 11.3; T3 mutant
  state 10.1; T1 direction 9.5; heart groups about 6 each; embryo cell state 5.8. Per board: T1 49, T3 39 (26.7 with
  s1's own factor stand-in), heart 24.5, embryo 20.

Our entries (beyond the baseline):
- **T1** s1: E10.5 real cells restricted to chamber categories at E10.5's own chamber mix (cardiomyocytes 26.5% of
  heart cells, blood about 24% of chamber cells); count-space shifts along Feng 2022's E9.5 -> E14.5 per-type
  direction at 2/5 for cardiomyocytes, endocardium, endothelium, epicardium; blood, macrophages, Early VIC keep E10.5
  expression; 23 genes held flat (sex-linked etc.). s2: same with cardiomyocytes at 60% of heart cells. On the
  backtest (E10.5 from E9.5) the stand-in sits about +1 over the floor; all of its gain is on de_direction.
- **Heart** s1: E10.5 cells and coordinates, half the E9.5 -> E10.5 per-lineage change in count space, core grown by a
  pre-registered size rule. Every checked metric within 1.3 points of copy_last. s2: the half trim of surrounding cells.
- **Embryo** s1: real cells half from E7.5 and half from E8.0, positions moved half-way along fused Gromov-Wasserstein
  partners, size by PCHIP of log RMS radius; about +11 over the floor on backtests. s2: E7.5's geometry.
- **T3** s1: WT E8.75 cells with the Gata4 per-gene factor (measured in non-Mesp1 compartments, applied in count space),
  Gata4's per-type depth step, and five Mesp1-lineage states depleted (OFT/RV-CM, aSHF, pSHF, aPHM, pPHM to 0.38x).
  s2: WT cells with the same states at 0.61x, no factor (or, if the organisers say a floored severity scores zero, the
  Mab21l2 factor without depth). The response's direction is the open problem: "on the published references every
  model gets the direction of the T3 response wrong" (organisers).

## 4. Hard constraints on any solution

- **Data windows** (by any route: directly, via a pretrained model trained on it, via an atlas containing it):
  T1 and heart: nothing measured after E9.5 up to and including E13.5 (exactly E9.5 allowed; after E13.5 allowed with
  disclosure). Embryo: nothing strictly between E7.25 and E8.0. Heart interpolation: nothing strictly between E8.25 and
  E8.75. Somite counts / Theiler stages are mapped to E-days by the midpoint rule (TS22 straddles E13.5; excluded).
  Genotype: nothing from the beta-catenin or Gata4 knockouts, another allele of those genes at a comparable stage, or
  a perturbation that phenocopies them. Public atlases spanning a window (MOCA E9.5-E13.5, MOSTA, the JAX time-lapse
  between E9.5 and E13.5) only after removing in-window cells by our own code.
- **Organiser ruling 2026-10-15:** generic pathway and gene-set databases are allowed except gene sets built from
  perturbing the held-out gene; gene-level results of losing beta-catenin or Gata4 in embryos are not allowed, even
  from a paper; phenotype-level knowledge is allowed. T3 gene directions come only from released data.
- **Post-exposure block** (team rule after an exposure): no change to cell numbers, sizes or dissection scope on any
  board; T3's composition is closed; T1's composition may change only by "the legal route": quotas at the frozen
  count, from a prior the captain writes from legal data, never importance weights.
- **A15 (working answer):** no new cross-modal input (e.g. scRNA-derived direction applied to MERFISH) enters an upload
  until the organisers answer.
- **Frozen rules:** the s2 rules for T1/heart/embryo and the pre-registered T3 rules move only by the captain's ruling.
- **Adoption rule:** a change is kept only if it gains (>= 0.3 points, paired over 10 seeds) on counting backtests it
  was not tuned or selected on: T1 = E10.5 predicted from E9.5 (and from E8.5), with a known composition confound (the
  inputs carry E9.5's mix, ~64% cardiomyocytes vs ~19% in the chamber-cut truth); heart = E10.5 from E9.5 and from
  <= E8.75; embryo = E7.5 (t = 1/3) and E7.25 (t = 2/3); T3 = the two cross-knockout folds (Gata4 -> Mab21l2,
  Mab21l2 -> Gata4); the Gata4 replicate fold never counts. A test score is never evidence.
- **Scores choose, never compute:** no inverting the skill map, no differencing scores, no leaderboard figure as a
  threshold.
- **Timeline:** s1 methods freeze 2026-11-27 00:00 UTC; uploads from then, each on the captain's go-ahead; final
  submissions through 2026-12-02.

## 5. Already covered: do not redo, go newer or deeper

- 2026-10-19 sweep (`docs/research/digest_2026-10-19.md`, screened; read it if useful): scNODE, scIMF, scTimeBench
  2026, Sunbear (Genome Res 2025), WOT, moscot, TrajectoryNet, MIOFlow, VGFM, DeepRUOT, OT-CFM, MMFM, Light SB,
  PRESCIENT, TIGON, CFGen, CellFlow; foundation models (UCE, TranscriptFormer, scPRINT, Nicheformer, GeneCompass,
  Mouse-Geneformer, C2S-Scale, Navigo/Chreode weights: plausibly contaminated; human-only scGPT, Geneformer,
  scFoundation, State SE, TF-Sapiens: allowed, near-useless); Ahlmann-Eltze 2025, Csendes 2025, Systema, VCC 2025 (only
  "won by hybrid statistical/deep systems tuned to the metrics"), GEARS/CPA/scGen/CellOT/biolord/PerturbNet,
  CellOracle, scTenifoldKnk; moscot spatiotemporal, POT FGW, Spateo morpho_align/morphofield, McCann interpolation,
  NicheFlow (NeurIPS 2025), Wasserstein Flow Matching (ICML 2025), STORIES, LION, SpaDiT. Legal data already used:
  Feng 2022 (GSE193346; E9.5 CD1 and >= E14.5), Feulner 2026 whole hearts (E14.5-P7, its TS22 sample out), JAX/Qiu
  2024 (embryos outside the window), ENCODE heart bulk from E14.5, GSE208162 (T3 only).
- 2026-11-05 cross-disciplinary research (`experiments/research_2026_11/REPORT.md`, `notes/registry.md`, 128 works):
  proper kernel scores, ranking under uncertainty, contest economics, fluctuation-response, TCF/LEF motif priors,
  herding/low-discrepancy draws, NB resampling, Gaussian-OT decoding, time warps, Potts sorting, assay depth.
- Pruned ideas with reasons: `experiments/harness/TRIES.md`, section "Pruned before running". A proposal that repeats
  one must say what is new.

## 6. Barred classes (never search for them, never record them)

- **A, hidden-target properties:** any size, radius, cell count, file size, composition or other measured property of
  a test target (E12.5 T1 or heart, E7.75 embryo, the beta-catenin KO samples); test-board limits or anchors.
- **B, gene-level knockout findings:** which genes change (or by how much) when beta-catenin/Ctnnb1 or Gata4 (incl.
  Gata4/Gata6 compound) is lost in embryos or cardiac lineages; gene sets or regulons built from perturbing them. Do not
  run searches like "Ctnnb1 knockout differentially expressed genes". Phenotype-level descriptions are allowed.
- **W, in-window findings:** measured data or results inside a window (e.g. anything measured at E10.5-E13.5 in mouse
  heart or embryo, or between E7.25 and E8.0), including qualitative stage-by-stage accounts that would set a prior.
  A paper that used in-window data may be cited for its METHOD only; record no in-window result.
- **L, per-board leaderboard figures:** any team's per-board or per-task score, rank, gap or range. Do not open the
  challenge leaderboard.
- **Y, partner material:** the Yukon/Eigen Labs team (yukon.org, its forms, posts, code, models).
- **Q, quarantined artefacts:** veckit PyPI 0.1.2 and veckit PR #6; the challenge's board files and endpoints
  (`panels/index.json`, `/challenge/boards`), the Data page; the module docstring and `de_genes` docstring of veckit's
  `common/core_metrics.py`; session logs, transcripts, harness-saved tool outputs and memory files under
  `~/.claude-account` or `~/.claude`.

If you meet a barred item, do not copy it; note only its class and the source in your compliance log.

## 7. Tools

Load WebSearch and WebFetch with ToolSearch ("select:WebSearch,WebFetch"). WebSearch mode "standard" by default,
"extended" for very recent (2026) or niche items. PubMed, bioRxiv and Consensus MCP tools are available through
ToolSearch ("pubmed", "biorxiv", "consensus"). Prefer HTML renderings (arxiv.org/abs, ar5iv, PMC, openreview.net/forum)
over PDFs; if a fetch saves a binary copy under the account directory, do not open it. No computation, no toy runs, no
repository edits; write only your own output file.
