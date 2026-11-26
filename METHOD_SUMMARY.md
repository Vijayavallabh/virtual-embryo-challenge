# Method summary

Team Mithril, Virtual Embryo Challenge 2026, Human Team track. The full working draft is in git history
(`git show c10e079:METHOD_SUMMARY.md`).

## Design

Every predicted cell is a real released cell. Each entry draws real cells to a target composition, edits their
expression only in count space (binomial thinning or Poisson top-up by per-type, per-metacell or per-gene factors
fitted on released stages), and places them in space. Nothing is trained jointly; the shared steps live in
`src/vec/core.py` and each task has its own head.

## Best solution (overall 205.62)

Each file was built twice with identical bytes, recorded with its sha256 before any score
(`experiments/builds/results/`), and chosen among uploaded files by headline only. Each command runs from the
repository root at the commit shown, as `uv run --frozen --no-sync python -I -B <script> <options> --out <file>`.

**T1:test, 74.75: `T1_s2rcMSBm10.h5ad`** (commit 0cdcc6f). Real E10.5 cells of the released dissected-heart library,
chamber categories only. Cardiomyocytes, endocardium, endothelium and epicardium take three count-space steps:
- each category's E9.5-to-E10.5 change, carried forward;
- the Feng et al. 2022 atlas's E9.5-to-E14.5 change, on cell groups recomputed from published markers, at 1.5 times
  the E10.5-to-E12.5 share, with technical genes held flat;
- for cardiomyocytes and endocardium only, within-type changes from entropic optimal transport.

`experiments/builds/build_test.py t1 --g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.0 --captain-ruling R51`

**T3:bcatenin, 64.22: `T3_c2_s2_g15.h5ad`** (726aaac). Real wild-type E8.75 cells, with second-heart-field,
outflow-tract and right-ventricle-lineage states mildly depleted. Each gene is scaled by one factor combining the
released Mab21l2 and Gata4 knockouts: Gata4's log factor is weighted 1.5, and the sum is held at the Mab21l2 factor's
size. Nothing about the held-out gene was used.

`experiments/t3_decomp/build_two_knockout.py build experiments/t3_decomp/results/gata4/gata4_decisions.json --wt data/raw/T2/E8.75.h5ad --run-dir <decomposition run of 2026-10-28> --board T3:bcatenin --seed 1000001 --form c2-s2 --gata4-exponent 1.5 --captain-ruling R76`

**T2:heart:test_extrap, 62.45: `T2heart_organ_ts125_tk_half.h5ad`** (8a8bc3a). Real E10.5 heart cells within the
heart-organ scope. Each lineage moves along its trend across the released heart stages, scaled by 1.25, with
endocardium and vascular endothelium handled separately. Coordinates grow by a pre-registered size procedure, and a
random half of the cells is kept.

`experiments/builds/build_test.py heart --entry organ_scope --trend-strength 1.25 --trend-form category --z-mode raw --captain-ruling R55`, then `experiments/builds/halve.py <file> <half> --seed 1000001`

**T2:embryo:test_interp, 70.85: `T2embryo_s2x_ovd.h5ad`** (ee22b75). Real cells, half from E7.5 and half from E8.0.
E8.0 positions are turned into E7.5's frame by a fused Gromov-Wasserstein coupling's rotation, and the cloud is scaled
to a pre-registered size. E8.0 cells are thinned to the depth midway between the two stages.

`experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact --shape-interp overlay --shape-interp-test-ok --expression depth --captain-ruling R75`

## Data and code

- **T1:** released E8.5, E9.5 and E10.5, plus the Feng et al. 2022 heart atlas (GEO GSE193346, E9.5 and E14.5 cells
  only). Published gene lists and marker statements are also used, and one atlas for diagnostics only. Full citations
  are in `EXTERNAL_DATA.md`.
- **T2 heart, T2 embryo and T3:** released challenge data and unlocked validation files only.
- **Models and code:** no pretrained models. Libraries only: POT, NumPy, SciPy, scikit-learn, pandas, AnnData,
  scanpy and h5py, plus the organisers' veckit for local scoring.
