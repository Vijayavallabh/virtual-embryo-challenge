# Data notes (structure only)

From `vec data inventory` on the DGX, 2026-10-18, re-run on 2026-10-28 with the six validation answers
(`runs/inventory.json` there, the first run kept as `runs/inventory.2026-10-18.json`; QC in `runs/qc.json` and
`runs/qc.csv`; log `runs/logs/intake_2026-10-28.log`). Checksums of all 15 files in `data/raw/SHA256SUMS` (the nine of 2026-10-18 are
unchanged; the first file is kept as `SHA256SUMS.2026-10-18`); the site publishes no checksums to compare with.
Structural facts only; no expression values. No test inputs are distributed.

## Files

| File | Cells | Genes | X | Layers | Notes |
|---|---|---|---|---|---|
| T1/E8.5_RNA.h5ad | 16,787 | 32,285 | CSC | - | obs: celltype (18 types); obsm: harmony UMAP |
| T1/E9.5_RNA.h5ad | 17,057 | 32,285 | CSC | - | 21 types |
| T2/E6.75.h5ad | 7,093 | 498 | CSC | counts, log1p | embryo setting; 18 types incl. Unknown |
| T2/E7.25.h5ad | 13,295 | 498 | CSC | counts, log1p | embryo setting |
| T2/E8.0.h5ad | 31,671 | 500 | CSC | counts, log1p | embryo setting; 26 types; has cm_celltype |
| T2/E8.25_late.h5ad | 58,716 | 500 | CSC | counts, log1p | heart setting; 33 types, many non-cardiac (extra-embryonic endoderm, neural tube, somites) |
| T2/E8.75.h5ad | 24,826 | 500 | CSC | counts, log1p | heart setting and T3 matched WT; 33 types, still much neural tube and NCC |
| T2/E9.5.h5ad | 53,742 | 500 | CSC | counts, log1p | heart setting and Mab21l2 matched WT; 22 types, heart-focused |
| T3/E9.5_mab21l2_ko.h5ad | 50,294 | 500 | dense | counts, log1p, normalized_counts | genotype all `mab21l2` (no littermate WT); 37 sample_ids; one slide |
| T1/E10.5_RNA.h5ad | 24,979 | 32,285 | CSC | - | validation answer; 23 types; obsm: harmony UMAP |
| T2/E7.5.h5ad | 11,665 | 500 | dense | counts, log1p | validation answer, embryo setting; 18 types; no cm_celltype |
| T2/E8.5.h5ad | 23,710 | 500 | CSC | counts, log1p | validation answer, heart setting; 33 types; 13 cm_celltype levels |
| T2/E10.5.h5ad | 167,864 | 500 | CSC | counts, log1p | validation answer, heart setting; 23 types; 9 cm_celltype levels |
| T3/E8.75_gata4_ko_1.h5ad | 42,257 | 500 | dense | counts, log1p, normalized_counts | validation answer; genotype all `gata4`, condition all `KO` (no littermate WT); 46 sample_ids; slide M164; 20 types; also `celltype_clean` |
| T3/E8.75_gata4_ko_2.h5ad | 46,396 | 500 | dense | counts, log1p, normalized_counts | validation answer; as replicate 1 (same slide M164), 44 sample_ids, no `celltype_clean` |

All obs names in the six validation files are unique. The Gata4 files carry the Mab21l2 KO's obs columns (slide,
section, embryo, sample_id, z_coord, time, genotype, perturb_batch, leiden ...) and obsm `spatial` (2D) beside
`spatial_3D`; the WT MERFISH files carry `spatial_2D` and `spatial_3D`.

Every file satisfies `expm1(X).sum(1) == 1e4` (log1p CP10k; the validation answers too, median 10,000). Stage cell counts are exactly 10x the board `ref_cells`
in `boards/index.json` (rounded), which confirms the server's 10% working copies.

## Quirks to handle in loaders

Settled by `experiments/qc/` (2026-10-19; `vec data qc`) and `experiments/coords/`:

- **z is section position x a per-file step**: 10 (E6.75, E7.25, Mab21l2 KO), 14 (E8.75), 15 (E8.0, E8.25_late, E9.5).
  E9.5's 98 raw z values are 40 levels exactly 15 apart; the extra values are floating-point copies (spread at most
  1.1e-13), one per separately placed tissue piece (E6.75 likewise: 65 raw values on 37 levels). Snap with
  `vec.spatial.snap_z`, never by rounding. xy are slide micrometres (each WT section is a pure translation of the slide
  coordinates; 30 of 36 KO sections are rigid copies, 6 were warped). The stored step is most likely the section pitch;
  use native coordinates and never rescale z across stages. Every file has 30-40 sections, so the z extent is set by
  acquisition, and the stack cuts the heart in E9.5 (start) and the KO (end).
- **Rows are stored section by section** in descending z (Spearman of row vs z -0.9994 to -0.9996): always subsample at
  random. T2/E7.5 runs the other way (ascending z, +0.9995).
- **Validation answers, intake QC (2026-10-28, `vec data qc`):**
  - z steps: 15 for T2/E8.5, T2/E10.5 and both Gata4 replicates; 14 for T2/E7.5.
  - Level counts (raw z values / levels): E8.5 81/25 and E10.5 157/44, floating-point copies again (snap with
    `vec.spatial.snap_z`); E7.5 34/34; Gata4 replicates 46/46 and 44/44.
  - Gaps: E10.5 has 9 single-level gaps, so from -433.272 to -163.272 its levels are 30 apart, not 15. The other new
    files have none.
  - Duplicates: no duplicated or interpolated records, no non-integer count rows, and no identical count rows in any of
    the six; T1 E10.5 has no identical rows.
  - Dense levels (cells against the mean of the two neighbouring levels; the QC reads about 2x as a hint of two
    sections on one level, while ordinary variation reaches about 1.8x): E10.5 -433.272 (1.73x); Gata4 replicate 1
    -585 (1.88x), -390 (3.26x) and -360 (1.71x); replicate 2 -615 (1.63x), -555 (2.35x) and -420 (1.64x).
  - Section and sample_id labels: one label per level and one level per label in both Gata4 files (the Mab21l2 KO
    fails this), so the 3.26x and 2.35x levels are each a single labelled section.
  - Shared rows between the two replicates are a pre-registered check of its own (runbook step 4.2).
- **Synthetic records** (none in E6.75, E7.25, E9.5, T1 or the validation answers):
  - copies, identical integer count rows in the same or an adjacent level, at the same xy (E8.0 443, E8.25_late 1,164,
    E8.75 317) or moved (2,460 / 563 / 518 more; the E8.75 ones are an exact affine image of the level above);
  - interpolated records with half-integer counts, each the average of two cells in neighbouring levels (E8.0 200,
    E8.25_late 2,558, E8.75 614);
  - Mab21l2 KO: sample_id 35 stored twice within level -320 (1,404 records) and sample_ids 18 and 19 one section at
    two levels, -210 and -220 (1,407 records); level -320 also holds sample_id 30, and the 30-unit gap to -290 means two
    sections are missing.
  `load_stage(drop_duplicates=True)` drops all of these for prediction inputs (E8.0 2,903 + 198, E8.25_late
  1,727 + 2,558, E8.75 835 + 614, KO 2,811); the emulator loads targets as released.
- **Counts:** `load_stage` keeps int32 counts when all values are integers and float32 otherwise (the interpolated
  rows), never truncating.
- **Mab21l2 KO:** non-unique obs names (all on the copies; made unique on load); coordinates not centred (mean
  -5439, -1464, -182; WT files are centred); different cell-type vocabulary (vCM, LPM, DPM, Paraxial Mesoderm,
  Branchial Arch...): map into the WT vocabulary with `vec.harmonise`; `configs/lineage_map.yaml` remaps the KO's LPM
  and DPM. Dense X.
- **Vocabularies are not harmonised across stages** (E8.25_late and E8.75 share label definitions; E9.5 does not;
  `cm_celltype` exists only in E8.0, E8.25_late, E8.75 and E9.5, where "Unknown" means "no subtype given"). Only CM
  transfers cleanly between E8.75 and E9.5; use the merged lineage level for cross-stage work.
- All MERFISH files carry a `counts` layer; T1 counts must be recovered from `.X`.
- A layer stored under the literal key `None` appears in every file (an artefact; ignore it).

## Loader checks (2026-10-18, `scripts/check_prep.py`)

- T1 counts recover exactly from `.X` (max distance from an integer < 6e-4); median library about 10.9K (E8.5) and
  11.1K (E9.5) UMIs per cell, p10-p90 about 6.2K-21K. T1 looks like whole-cell 10x data, not nuclei.
- MERFISH `counts` layers reproduce `.X` exactly (max difference 1e-6). Median counts per cell: E6.75 30, E7.25 35,
  E8.0 212, E8.25_late 231, E8.75 130, E9.5 75, Mab21l2 KO 95.
- The Mab21l2 KO is about 27% deeper per cell than its matched WT E9.5 (95 vs 75): a systematic, non-genotype shift
  that the T3 "systematic component" must carry.
- Coarse composition (lineage map): the heart-setting files differ a lot in scope (E8.25_late is 4% core cardiac, E8.75
  17.5%, E9.5 37%); E8.25_late and E8.75 are dominated by non-cardiac tissue.
- The KO and WT vocabularies do not line up even at the coarse level (pharyngeal mesoderm exists only under WT names),
  so KO-vs-WT composition needs classifier-based relabelling on the shared panel, not label matching.
