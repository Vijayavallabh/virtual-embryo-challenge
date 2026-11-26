# Method summary

Team Mithril, Virtual Embryo Challenge 2026, Human Team track. The full working draft (data provenance, compliance
notes, every departure) is in git history: `git show c10e079:METHOD_SUMMARY.md`.

## The shared design: transport-then-adjust on real cells

No model is trained to generate cells: every predicted cell is a real released cell. Each entry runs one recipe and
switches on the steps it needs:

1. allocate a target composition to whole cell counts;
2. draw real cells from released stages;
3. change their expression in count space only (per-type, per-metacell or per-gene factors applied by binomial thinning
   or Poisson top-up of recovered counts), or leave it unedited;
4. write log1p CP10k;
5. where positions must move, couple two stages with positions (fused Gromov-Wasserstein coupling and its rotation).

The steps live in a common core (`src/vec/core.py`); each task has its own head. Task 3 runs the same steps in its own
pre-registered script, with five functions pinned by tests to the core's copies. Factors are fitted on released stages;
nothing is trained jointly across tasks.

## The submitted files

**T1, heart E12.5** (`T1_s2rcMSBm10.h5ad`). Real E10.5 cells of the released dissected-heart library, chamber
categories only, at that library's own composition. For cardiomyocytes, endocardium, endothelium and epicardium, three
count-space steps, each keeping a cell's expected depth:
1. each category's mean change from E9.5 to E10.5, carried forward (damping 1.0);
2. the E9.5-to-E14.5 log ratio of mean expression from the Feng et al. 2022 heart atlas, on cell groups recomputed from
   published marker genes, at 1.5 times the E10.5-to-E12.5 share of that interval, with fixed lists of technical genes
   held flat;
3. for cardiomyocytes and endocardium, within-type changes fitted by entropic optimal transport from E9.5 to E10.5,
   damped by half.

**T2 heart, E12.5** (`T2heart_organ_ts125_tk_half.h5ad`). Real E10.5 heart cells within the heart-organ scope.
Expression moves along each lineage's least-squares trend of log mean expression across the released heart stages
(endocardium and vascular endothelium with separate trends), scaled by 1.25. Coordinates are the E10.5 tissue's, grown
by a pre-registered size procedure. A uniform random half of the cells is kept.

**T2 embryo, E7.75** (`T2embryo_s2x_ovd.h5ad`). Real cells, half from E7.5 and half from E8.0. E7.5 cells keep their
positions; E8.0 cells keep theirs, turned into E7.5's frame by the rotation of a fused Gromov-Wasserstein coupling. The
cloud is scaled to a pre-registered size, and E8.0 cells are thinned to the depth midway (geometric mean) between the two
stages.

**T3, beta-catenin knockout, E8.75** (`T3_c2_s2_g15.h5ad`). Real wild-type E8.75 cells, drawn to a composition in which
second-heart-field, outflow-tract and right-ventricle-lineage states are mildly depleted (pre-registered rules). Each
gene's counts are scaled by one per-gene factor that combines the released Mab21l2 and Gata4 knockouts' factors against
wild type: Gata4's log factor is weighted 1.5, and the sum is held at the Mab21l2 factor's size. Nothing about the
held-out gene was used.

**Choosing files.** Every candidate was fully built before upload. Among the uploaded candidates, files were chosen by
leaderboard headline score only; no prediction was computed from scores.

## Data

- **T1:** released stages E8.5, E9.5 and E10.5. External:
  - Feng W, Bais A, He H, et al., Nat Commun 13:7960 (2022), GEO GSE193346, E9.5 and E14.5 cells only;
  - the genes held flat, from van den Brink SC, et al., Nat Methods 14:935-936 (2017), and Young MD, Behjati S,
    GigaScience 9:giaa151 (2020); Bakken TE, et al., PLoS One 13:e0209648 (2018) informed that step's premise;
  - marker-gene symbols for the cell-group labeller, from eight papers;
  - Qiu C, et al., Nature 626:1084-1093 (2024), for diagnostics only.
  
  Full citations are in `EXTERNAL_DATA.md`.
- **T2 heart:** the released heart stages and the unlocked E10.5 validation file. No external data.
- **T2 embryo:** the released embryo stages, including the unlocked E7.5 validation file and E8.0. No external data.
- **T3:** released wild-type E8.75 cells, the released Mab21l2 knockout and the released Gata4 validation knockout files.
  No external data.

No pretrained model is used. Code: open-source libraries called as libraries, namely POT (optimal transport), NumPy,
SciPy, scikit-learn, pandas, AnnData, scanpy and h5py. The organisers' veckit is used for local scoring only.
