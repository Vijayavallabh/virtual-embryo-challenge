# Assumption ledger

The bets the entry rested on, with their status at the end of the season (ledger started 2026-10-18). Each one's
test, evidence and log are in git history (`git show c10e079:ASSUMPTIONS.md`); code docstrings cite these IDs.

| ID | Assumption | Status |
|---|---|---|
| A1 | All released files (training stages, P3 validation answers, test inputs) can be downloaded to the DGX | confirmed |
| A2 | The local board emulator reproduces the live board | revised |
| A3 | Each test board uses the last released stage as copy_last and DE reference (T1 and T2 heart: E10.5; T2 embryo: E7.5) and publishes truth_cells | open |
| A4 | Real-cell resampling plus count-space edits beat additive or decoded generation on mmd_u and variogram | open |
| A5 | Legal post-E13.5 external heart data gives T1 a better E10.5->E12.5 direction than internal momentum | revised |
| A6 | The E12.5 T1 composition (heart chamber only, per the organisers' dissection note) can be predicted well enough to halve mmd_u versus copy_last | revised |
| A7 | FGW displacement interpolation of real E7.5 and E8.0 tissue predicts the embryo E7.75 board well | revised |
| A8 | The E12.5 heart size and scope can be predicted from released stages well enough to beat copy_last on scale_log_ratio | revised |
| A9 | Legal external datasets (Feng 2022 by per-cell stage, Feulner 2026 (not its sample coded Theiler stage 22, X006), JAX embryos 1-22 and 54-74 by embryo ID, extended gastrulation atlas; label-only E9.5 sources such as Spateo E9.5 only after the organisers' somite convention) can be acquired and window-filtered on the DGX within 150 GB | confirmed |
| A10 | The released knockouts' per-gene factor (non-Mesp1 compartments vs WT) transfers to the beta-catenin knockout, i.e. the knockout series shares processing or the factor is biology | killed |
| A11 | Depleting second-heart-field / OFT-RV-lineage cells (broad phenotype background, allowed per the 2026-10-15 ruling) gives the correct sign of the beta-catenin pseudobulk response | killed |
| A12 | A Wnt-responsive programme learned from wild-type T1 scRNA adds signal to within-type KO responses | open |
| A13 | Submitted cell count moves T2 board scores, and matching the truth-half density is near-optimal | revised |
| A14 | Final evaluation keeps the current metric definitions (severity_slope, scale_log_ratio), not veckit PR #6 | open |
| A15 | Released data may cross tasks in P3 (T1 scRNA for T2/T3; MERFISH for T1) | open |
| A16 | One shared transport-then-adjust design (`vec.core` on T1 and T2, feat-020; on T3 the pre-registered script's own primitives, five pinned to copies in the core) builds every board within 3 points of overall standing of the best specialised pipelines (Generality Award claim) | revised |
| A17 | Backtest gains transfer to the hidden test boards | revised |
| A18 | The blind entry places in the Human top 6 on the test boards | open |
| A19 | Harmonised composition shifts between a released knockout and its matched wild type measure genotype effects well enough to set T3 scenario magnitudes | confirmed |
| A20 | Released validation answers (T1 E10.5, T2 heart E8.5 and E10.5, T2 embryo E7.5, T3 Gata4) may be used as training data and model inputs for the test predictions once the final phase opens | confirmed (user's reading) |
| A21 | The organisers count our shared core with task heads (a 32,285- or 500-gene encoder, a coordinate head, a genotype head) as one shared architecture for the Generality Award | open |
| A22 | A board's weakest weighted component on its counting backtests is also its weakest at test (the captain's guiding principle, 2026-11-04) | open |
| A23 | The beta-catenin knockout file shares Gata4's shallow depth against WT E8.75 (U2), so T3 s1's depth step matches the truth | open |
| A24 | In the final evaluation a floored severity (wrong sign, beta <= 0 or r^2 <= 0.01) scores skill 0.5, as the scorer gives it, not 0 (the "undefined" reading our questions attributed to the Evaluation page; no record quotes it, `docs/PLAN.md` decision log 2026-11-11 (e)) | open |
| A25 | The hidden E7.75 embryo sample was measured like E8.0: it shares E8.0's per-gene detection profile and depth (the bet behind request R-A, embryo s2 = fgw_aligned_series) | open |
| A26 | The hidden E12.5 T1 library is a dissected-heart sample no shallower than about 0.7 of E10.5's library 3 (the bet behind request R-B, T1 s1 on library 3's cells, and since row R7, at about half strength, T1 s2 on the base half-way to library 3) | open |
| A27 | Part of the beta-catenin response scored against WT E8.75 is WT E8.75's own per-gene deviation, which the other released wild-type samples estimate well enough to clear the scorer's severity gate (the bet behind adding a wild-type-only factor to T3's second upload; not requested) | open |
| A28 | The DGX is reachable on the upload days (27 to 30 Nov) and its numerical libraries and processor are those of 2026-11-13, so the freeze-morning rebuild reproduces the ten reference upload files | open |
| A29 | The heart's practice gain of the trend over the last step carries to the test, whose last step (E9.5 to E10.5) crosses the change from sectioned embryo portions to dissected hearts | open |
