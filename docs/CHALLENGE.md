# Virtual Embryo Challenge (VEC), NeurIPS 2026: challenge reference

Compiled 2026-10-18 from virtualembryo.ai/challenge (all subpages), the accepted NeurIPS competition proposal (PDF),
and the official local scorer `aristoteleo/veckit` (MIT). Our team competes in **Track 1 (Human Team) only**.

Sources (keep in sync when the site changes):
- Challenge site: https://virtualembryo.ai/challenge (tasks, data, evaluation, baselines, rules, FAQ, timeline).
- Accepted proposal with the metric definitions and Table 1 weights: `docs/reference/neurips-2026-competition-proposal.pdf`.
- Nature Methods comment (Cao, Lu, Qiu 2026): `docs/reference/towards-predictive-virtual-embryos.pdf`.
- Board specs, gene panels and T1 composition: `boards/`. The validation boards keep their floor and ceiling anchors.
  Since 2026-11-04 (634b5fd) it also holds the P3 test boards, taken in with `vec boards update --from-file` from a
  copy the captain saved (nobody opened the saved file), their anchors stripped in memory. The test boards' cell
  limits are read by code only: nobody opens `boards/index.json` or quotes them (`docs/DATA_POLICY.md`). The online
  refresh still answers HTTP 403 (since 2026-10-27, not worked around; asked in the fourth organiser email).
- Official local scorer: https://github.com/aristoteleo/veckit (MIT), pinned in `pyproject.toml`.
- Partner rules (Yukon), https://virtualembryo.ai/challenge/partner-rules, and the Privacy Policy, both published
  2026-10-20. Checked 2026-10-20 17:10 UTC: the Rules are unchanged since 2026-10-18 apart from a Yukon banner
  ("Ordinary-team rules are unchanged"); the partner exceptions "apply to the partner category only". On 2026-10-28
  (a screened read of the site) the final-phase Rules had changed: s.10 and s.11 now cover the opt-in validation
  answers and the P3 start time (Timeline).

## Timeline (compiled 2026-10-18, updated to 2026-10-28; schedule change of 2026-10-20 included)
- P2 (until the final phase opened on 2026-10-27): validation leaderboard, 8 scored submissions/day/task. Starter kit +
  reference baselines are public (Colab: https://colab.research.google.com/drive/11kHZjd4IawoC6uwCqiHNuDGAYpoPlUQO).
  P2 standings were frozen at the P3 cutoff (the "Frozen P2" board); the live Validation board stays, but our team can
  no longer submit to it (the unlock, below).
- **P3 started 2026-10-27**, postponed by a week from 2026-10-20 (organisers' email to registered teams, 2026-10-20,
  citing a new partnership with Eigen Labs; scoring, submission limits and the other decisions for later dates are
  unchanged). Files committed before 2026-10-20 (pre-registered procedures, configs, reports) say 2026-10-20: they mean
  this release. P3 opened on 2026-10-27 at 09:00 America/Los_Angeles (16:00 UTC; Rules s.11); the test leaderboard
  opened then. As it turned out (2026-10-28), no test inputs are distributed, and the validation answers come only
  after a team member unlocks them on the Data page. Unlocking permanently ends P2 submissions for the whole team,
  even without a download, and uses no P3 attempt. Our team unlocked them on 2026-10-28 at 09:32 UTC. Rules s.10 makes
  those official files an exception to the held-out-data restrictions "solely for preparing P3 predictions"; the user
  read this as allowing them as model inputs (A20, passed on 2026-10-28; organiser confirmation asked in the fourth
  email). **Only 2 official scored submissions per board for the whole phase**; Rules s.11 says the board ranks a team
  by its best. The account page adds: "In the final phase the captain picks which submission counts" (a "Final" mark
  per submission, until 2026-12-02). The user is the captain; what counts on a board with no mark is not stated (asked
  in the fourth email). Each official submission "is scored and published when it is made" (Rules s.11), so the test
  leaderboard is not hidden until the end. Format checks unlimited and free; a separate daily upload limit bounds
  storage (`vec submit status` shows it; our format checks of 2026-10-20 did not count toward it).
- Final submissions due **through 2026-12-02 (UTC)** (the time zone comes from the organisers' email of 2026-10-20
  22:03 IST, which gives no hour); official evaluation 2026-12-04; winners at NeurIPS 2026-12-11 (unchanged by the
  postponement). Per the same email and the partner rules (s.5), P3 methods, prediction files and detailed evaluation
  feedback stay private from other participants until the 2026-12-11 announcement (the organisers may process
  submissions for scoring and verification). The sentence sits in Yukon's final-phase section and is absent from the
  Rules, but its wording is not limited to Yukon, so we follow it.
- **Yukon partner team** (announced 2026-10-20; the organisers' email of 2026-10-20 22:03 IST turned it into a
  collaborative track): Eigen Labs runs a Yukon team (https://www.yukon.org/) whose authorised members share
  executable models, fitted weights, methods and evaluation feedback with each other. It appears on the leaderboard as
  an unranked reference, takes no official prize position and displaces no eligible team. Its separate prize pool
  ($10K-$30K, administered by Eigen Labs) depends on its final reference position among VE teams; half goes to the
  top three Yukon results (50/30/20), half by contribution points. Entry is invite-only for now (an interest form asks
  for the applicant's current VE score; admitted members first submit the methodology behind that score). The
  collaborative phase was to end before P3 opened on 2026-10-27; members then build their P3 submissions independently.
  Moving from an official team is one-way for whoever moves (no return to an ordinary team this edition, and no joining
  or creating one, even after leaving Yukon; earlier submissions and scores stay with the original team, and a team
  left empty is archived and excluded from later rankings and awards: partner rules s.3) and makes that person
  ineligible for the Challenge's official prizes. Our entry is an official Human Team. Yukon team material is
  quarantined (`docs/DATA_POLICY.md`).

## Data (~1M cells, 11 stages E6.75-E12.5, mouse; generated by Qingquan Zhang & Neil Chi, UCSD; unpublished)
### Task 1: scRNA (dissociated, heart-centred dissection), whole transcriptome
- Files: E8.5_RNA.h5ad (571 MB, 16,787 cells, 18 types), E9.5_RNA.h5ad (590 MB, 17,057 cells, 21 types). E10.5 = val (647 MB), E12.5 = test ([redacted 2026-11-04: A]). No E9.25. E7.75 scRNA exists but is withheld (it is T2's embryo test stage).
- .X float32 log1p-normalised (log(1+1e4*c/sum c)), 32,285 genes, obs["celltype"]; no coordinates.
- Annotation vocabulary is NOT harmonised across stages (Endocardium at E9.5 -> Endocardium-1/-2 later). E8.5 and E9.5 share 11 of 18/21 types; E9.5->E10.5 shares 13 of 23 types. Dissection scope varies by stage (Neural Tube 4.6% at E8.5, 0% at E9.5; paraxial mesoderm and EXEM likewise). E9.5 includes Hepatocyte 7%, ST 5%, pPHM 7%, Proepicardium 5%, Endocardium 9%.
### Task 2: 3D MERFISH (serial sections reassembled into 3D), 500-gene panel, log-normalised, obsm["spatial_3D"] per-embryo local frame (not registered across stages), obs["celltype"]
- Heart setting: E8.25 train (E8.25_late.h5ad 225 MB), E8.5 val-interp ([redacted 2026-10-19: A]), E8.75 train (84 MB; also T3 WT reference), E9.5 train (164 MB), E10.5 val-extrap ([redacted 2026-10-19: A]), **E12.5 test-extrap** ([redacted 2026-10-19: A]). In the final phase every heart stage except E12.5 is training input (once a team unlocks the validation answers; ours did on 2026-10-28).
- Embryo setting (whole embryo, gastrulation): E6.75 train, E7.25 train, E7.5 val-interp, **E7.75 test-interp**, E8.0 train. Embryo panel = 498 genes (Casp4, Pnliprp1 missing at E6.75/E7.25).
### Task 3: conditional KOs, 3D MERFISH, same schema; Mesp1-Cre driven (mesoderm lineage)
- Train: Mab21l2 KO @ E9.5 (458 MB). Val: Gata4 KO @ E8.75 = Mesp1-Cre; Gata4 F/F; Gata6 F/+ (2 replicates, 374+409 MB; answers opt-in in P3, unlocked by our team on 2026-10-28). **Test: beta-catenin KO @ E8.75 = Mesp1-Cre; Ctnnb1 F/F** (2 replicates; broadly expressed gene, diffuse effect, "hardest").
- Matched WT @ E8.75 (84 MB) and E9.5 (164 MB). Every DE metric is computed against the matched WT at the same stage.
- **Ctnnb1 is NOT on the 500-gene panel** (Gata4, Gata6, Mab21l2 are). Classic Wnt targets [redacted 2026-10-19: B] are absent; [redacted 2026-10-19: B] are present.

## Tasks, boards and what is submitted
- T1: predict a SET of cells (AnnData, .X over all 32,285 genes in released order, 1,000-5,118 cells for val; no labels). Train E8.5, E9.5 -> val E10.5, test E12.5 (nothing observed between E10.5 and E12.5).
- T2: AnnData with .X on the board's panel + obsm["spatial_3D"] (n,3), any frame. Boards: heart val_interp (E8.5), heart val_extrap (E10.5), embryo val_interp (E7.5); test: embryo E7.75 (interp), heart E12.5 (extrap). Each board scored separately; task score = mean over its boards.
- T3: AnnData for the mutant (not the delta), .X 500 genes + obsm["spatial_3D"]. Board per held-out KO.
- Gene list must match the board panel exactly in order (use the published panel files). Cell count: any number within [min_cells, max_cells]; count is not ranked. A Task-1 cap is 3x max(truth half, ref) with floor 5,000.
- Scorer subsamples every stage to 10% first; target is then split in half (one half = truth, other half = ceiling prediction). E.g. T1 val: truth_cells 1,249, ref_cells 1,706 (10% of E9.5 = reference stage). The heaviest metrics subsample anyway (mmd_u 2,000 cells; energy distance and variogram 1,500).

## Metric panel (exact code in veckit; weights from proposal Table 1)
Notation: pb = pseudobulk (mean over cells) of log-normalised X; ref = preceding stage (T1/T2) or matched WT (T3); dp = pb(pred)-pb(ref); dt = pb(true)-pb(ref).
- **de_score**: truth DE genes = two-sided Mann-Whitney U (true vs ref cells) + BH at 0.05 + |lfc| >= 0.25 (in practice just |lfc|>=0.25), split into n_up/n_dn. Prediction contributes its top-n_up and bottom-n_dn genes by its own dp; signed overlap ov. de_score = (ov - chance)/(1 - chance), chance = max(ov using +pb(ref), ov using -pb(ref)) (expression-magnitude null). Rank-based: rescaling dp is a no-op; negation is punished. If std(dp) < 0.01*std(dt): score 0.
- **de_direction**: rank partial correlation of dp vs dt controlling for rank of pb(ref), over ALL genes (32,285 for T1). 0 for no change.
- **severity_slope** (T3 only): log of zero-intercept OLS slope of dp on dt over truth DE genes; worst = log(1e-3) if beta<=0 or R^2<=0.01 or no response. Scored on |value| (0 is ideal).
- **mmd_u**: unbiased multi-kernel RBF MMD^2 on 2,000 subsampled cells in a 30-PC space fitted on the TRUTH only; bandwidths 0.25-4x truth median heuristic.
- **variogram**: 20,000 random gene pairs, 1,500 cells; compare E|x_i - x_j|^0.5 between pred and truth; mean squared difference. Constrains gene-gene co-variation and, on 32K genes, per-gene sparsity/level.
- **d2_shape**: 1-D Wasserstein between distributions of within-cloud pairwise distances / RMS radius (alignment-free, reflection-blind).
- **occupancy_dice**: canonicalise (centre, PCA-rotate, scale by RMS radius, right-handed), voxelise on shared 16^3 grid over [-3,3], Dice maximised over 4 proper sign flips.
- **scale_log_ratio**: log(RMS radius pred / RMS radius true); scored on |value|. Coordinates' absolute units matter here.
- **neighborhood_mmd**: mmd_u between distributions of 15-NN spatial-neighbourhood-averaged expression profiles.
- Reported but NOT ranked: energy_distance, sliced_wasserstein, count_log_ratio, pb_rel_err, library_size_ratio, variance_ratio, composition_JSD (frozen logistic-regression probe trained on the truth), pseudobulk_pearson, morans_I_agreement.
- **Skill map** per metric: d_floor = |ceiling - floor|; skill = d_floor/(d_floor + d(value)) where d = signed distance from ceiling (worse = positive); clipped at 1. Floor -> 0.5, ceiling -> 1.0, never negative. NaN/undefined -> skill 0. Floor = copy_last (T1/T2: last observed stage verbatim, incl. its coordinates) or wt_identity (T3). Ceiling = half-vs-half of the target.
- **Weights** (task score S_t = 100 * sum_g w_g sum_m w_m|g skill(m)):
  - T1: de_score .25; de_direction .25; mmd_u .30; variogram .20.
  - T2: expression change .25 (de_score .5, de_direction .5); cell-state distribution .25 (mmd_u .6, variogram .4); morphology & growth .25 (d2_shape, occupancy_dice, |scale_log_ratio| 1/3 each) **with a structure gate: shape credit multiplied by neighbourhood-structure skill** (added after ctrl_random_cube beat the floor); local spatial organisation .25 (neighborhood_mmd).
  - T3: de_score .30; de_direction .25; |severity_slope| .25; mutant-state distribution .20 (mmd_u .6, variogram .4). Shape terms NOT scored for T3.
- Overall standing = sum of the three task scores (unattempted task = 0). Prize eligibility requires beating the floor overall (>150 if all three tasks entered). Ties: hidden test, then DE group, then distributional group.
- Validation anchors (`boards/index.json`), e.g. T1 val: de_score ceiling 0.846, de_direction ceiling 0.790, mmd_u floor 0.0836 / ceiling 0.0041, variogram floor 0.00522 / ceiling 0.000158. T3 gata4: de_score ceiling 0.870, de_direction 0.925, severity floor -6.91 / ceiling -0.009, mmd_u floor 0.0314 / ceiling 0.0004.

## Reference rows (published baselines; numbers are on the live leaderboard at https://virtualembryo.ai/challenge/leaderboard, fetched client-side from https://api.virtualembryo.ai)
- copy_last (T1/T2 floor = 50); wt_identity (T3 floor = 50). ceiling = 100 (clipped).
- The published baseline of each final-phase test board, which every entry is measured from (the captain's guiding
  principle, 2026-11-04; `AGENTS.md`):

  | Test board | Counts in Overall | Published baseline (scores 50) | What it submits |
  |---|---|---|---|
  | T1:test (E12.5 scRNA) | 1 | copy_last | the last released stage verbatim (E10.5, expected) |
  | T2:heart:test_extrap (E12.5 MERFISH) | 0.5 (one of the two boards in the T2 mean) | copy_last | E10.5 (expected), with its coordinates |
  | T2:embryo:test_interp (E7.75 MERFISH) | 0.5 (one of the two boards in the T2 mean) | copy_last | E7.5, the earlier bracket (expected), with its coordinates |
  | T3:bcatenin (beta-catenin knockout, E8.75) | 1 | wt_identity | the matched wild type at E8.75, unchanged |

  The reference stages are assumption A3's expectation: the P3 specs publish no ref_cells (D1, ae6d6b3). The floor
  scores exactly 50 because the scorer measures its anchors on its own working copy of the reference stage; a copy we
  submit is another draw and can score below 50 (the copy rows in `experiments/harness/TRIES.md`, "Where the entries
  stand against the published baselines").
- pseudobulk_shift (T1/T2): per-cell-type constant velocity between the last two stages added to the last stage's real cells, clipped at 0; types absent earlier stay put; coordinates unchanged.
- shift_ode (T1), neural_ode (T1: autonomous field in gene space, flow-matched on a global entropic-OT coupling; does NOT beat the constant shift), spatiotemporal_ode (T2: its score comes almost entirely from an analytic growth-rate extrapolation), dynode_flow (external flow-matching reference model, unmodified; likely the Qiu lab's Navigo, github.com/aristoteleo/Navigo-release), perturb_ode (T3), shift_transfer (T3: add the training KO's pseudobulk delta to the target WT), gene_ko (shift_transfer + zero the KO gene).
- Organisers' notes: "On the published references every model gets the direction of the [T3] response wrong." "Clearing the [T1] floor plausibly needs supervision across the whole stage ladder, a non-autonomous field, and a treatment of populations appearing rather than only existing ones moving." Fitting constant-velocity extrapolation lands below copy_last on every setting. On interpolation boards copy_last copies a stage [redacted 2026-10-19: A] too small by volume.
- Adversarial controls (must stay below floor): ctrl_scale_ref (x2), ctrl_shrink_ref (x0.99), ctrl_one_cell, ctrl_random_cube, ctrl_squashed_ref, ctrl_scale_wt, ctrl_shrink_wt (x0.75), ctrl_random_dir.

## Rules that constrain methods (Official Rules s.10; Human Team)
- External public data, pretrained models and published code allowed if licences permit and every source is disclosed. Expressly encouraged for T1.
- Forbidden: measured data from a held-out stage/genotype by any route (directly, via a pretrained model trained on it, via a public dataset containing it).
- Held out: T1 E10.5, E12.5; T2 embryo E7.5, E7.75; T2 heart E8.5, E10.5, E12.5; T3 Gata4 and beta-catenin KOs at E8.75.
- Final-phase exception (Rules s.10, as read 2026-10-28): the official validation files, once a team unlocks them on the Data page, may be used "solely for preparing P3 predictions"; external data from the validation stages and genotypes stays prohibited, and the test windows are unchanged. The user read this as allowing them as model inputs too (A20, passed 2026-10-28; `docs/DATA_POLICY.md`).
- Extrapolation targets (T1 and T2 heart): external data **excluded from after E9.5 up to and including E13.5**; exactly E9.5 allowed; **anything after E13.5 allowed** with disclosure of source and stages.
- Interpolation targets: external data strictly between the bracketing released stages is excluded: T2 heart strictly between E8.25 and E8.75; T2 embryo strictly between E7.25 and E8.0. Bracketing stages themselves allowed. Somite counts / Theiler stages are mapped onto this scale.
- Genotype: data from the held-out KOs, from another allele of the same gene at a comparable stage, or from a perturbation that phenocopies them, counts as held out.
- General resources spanning held-out windows may be used only after removing the in-window cells (and describing the removal).
- Score returns may be used to choose among predictions, not to compute them (no inverting the skill map / differencing scores to solve for target properties).
- Code sharing across teams prohibited (FAQ; Rules s.4 bars attempts to access another team's work). One team per person. Rate limits per team/task/day. File <= 1200 MB.
- Generality Award ($8K, one across both tracks, in addition to any placed prize; Rules s.14): a team entering >= 2 tasks answers them with one shared architecture, the same under every task it enters (separate training per task and task-specific input and prediction heads allowed; an unrelated model per task does not qualify). Judged on the method summary, which must ask to be considered and describe what is common, what differs and where the boundary falls; among eligible entries it follows final standing; the floor condition applies; the organisers may ask to see code privately. Unawarded, it goes to the Community Contribution Award pool.
- Winners must provide a written method report within 14 days; results must be reproducible; organisers may ask for code privately.
- Prizes ($112K in all; Prizes page, 2026-10-20): per track 1st $8K, 2x 2nd $5K, 3x 3rd $3K, filled only by eligible entries (beating the floor), with fewer than six eligible the rest not awarded; the Generality Award above; travel awards $30K (15-20 grants for early-career researchers to attend the NeurIPS workshop; no call posted yet); Community Contribution Award $20K (up to $200 per public contribution such as tools, tutorials, notebooks, documentation or fixes; anyone may file, competing teams included; rolling until the fund is committed or 2026-12-11; unawarded prize money is added to it). Leading teams of both tracks are invited to co-author the organisers' Nature Methods analysis (the publish page).

## Our compute
- PrakashDGX_H2: 8x H100 80GB (several GPUs busy with the user's other jobs; typically 2+ free), 1.6 TB free on the shared /home (2026-10-28), uv available, Python 3.10 system (uv-managed 3.11/3.12 fine), CUDA driver 13.0. Local machine is overloaded: editing/git only.
- BLAS threads on the DGX (2026-10-28): scikit-learn/numpy jobs need OPENBLAS_NUM_THREADS set (<= 64); without it a run crashed ("OpenBLAS warning: precompiled NUM_THREADS exceeded", then "malloc(): corrupted top size"). The repo's shell runners set OMP/OPENBLAS/MKL_NUM_THREADS themselves (at most 16 per process; `docs/COMPUTE.md`).
- Environment policy: uv only (no conda/pip venvs).
