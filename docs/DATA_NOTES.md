# Data notes (structure only)

From `vec data inventory` on the DGX, 2026-10-18, re-run on 2026-10-28 with the six validation answers
(`runs/inventory.json` there, the first run kept as `runs/inventory.2026-10-18.json`; QC in `runs/qc.json` and
`runs/qc.csv`; log `runs/logs/intake_2026-10-28.log`). Checksums of all 15 files in `data/raw/SHA256SUMS` (the nine
of 2026-10-18 are unchanged; the first file is kept as `SHA256SUMS.2026-10-18`); the site publishes no checksums to
compare with. Structural facts only; no expression values. No test inputs are distributed.

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

Every file satisfies `expm1(X).sum(1) == 1e4` (log1p CP10k; the validation answers too, median 10,000). Stage cell
counts are exactly 10x the board `ref_cells` in `boards/index.json` (rounded), which confirms the server's 10% working
copies.

## Quirks to handle in loaders

Settled by `experiments/qc/` (2026-10-19; `vec data qc`) and `experiments/coords/`:

- **z is section position x a per-file step**: 10 (E6.75, E7.25, Mab21l2 KO), 14 (E8.75), 15 (E8.0, E8.25_late, E9.5).
  E9.5's 98 raw z values are 40 levels exactly 15 apart; the extra values are floating-point copies (spread at most
  1.1e-13), one per separately placed tissue piece (E6.75 likewise: 65 raw values on 37 levels). Snap with
  `vec.spatial.snap_z`, never by rounding. xy are slide micrometres (each WT section is a pure translation of the slide
  coordinates; 30 of 36 KO sections are rigid copies, 6 were warped). The stored step is most likely the section pitch;
  use native coordinates and never rescale z across stages. Every training file has 30-40 sections (the validation
  answers have 25 to 46 z levels, below), so the z extent is set by acquisition, and the stack cuts the heart in E9.5
  (start) and the KO (end).
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
  - Centring: the three WT validation files are centred on their mean like the WT training files; both Gata4 files
    are not centred, like the Mab21l2 KO (`centred_on` in `runs/qc.csv`).
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
  of the training files, `cm_celltype` exists only in E8.0, E8.25_late, E8.75 and E9.5, where "Unknown" means "no
  subtype given"; the validation answers E8.5 and E10.5 have it too, E7.5 does not). Only CM transfers cleanly between
  E8.75 and E9.5; use the merged lineage level for cross-stage work.
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

## Data Protocol page (organisers, posted 2026-10-29; read through a screening agent on 2026-10-31)

https://virtualembryo.ai/challenge/data-protocol, released stages only (the screen withheld 24 items about how the
held-out stages were collected and one about the test boards' format; the page itself withholds every file property of a
held-out target until it is released).
- **Coordinates:** "the physical unit of the `spatial_3D` coordinates" and the registration method are "not stated here
  yet"; coordinates are per-embryo local frames, not registered across stages. "The z coordinate stored in
  `spatial_3D` steps by the amount in the last column" of a per-file table that matches the steps above (10, 14 or 15),
  and "The files are released as they are, without rescaling, so use the stored coordinate as it is." Sections are
  14 µm apart at every released stage, so stored z is not micrometres in every file. Our handling (native coordinates,
  z never rescaled across stages) is what the page asks for; this settles follow-up question 5 as far as it can be.
- **Labels:** annotated "stage by stage" by clustering and manual naming; "the vocabulary is not harmonised across
  stages", and no cross-stage correspondence table exists for `celltype` or `cm_celltype` ("not used in scoring";
  "Unknown" there means not given a subtype). A label missing from a stage may exist under another name.
- **Sampling:** E8.0 whole embryos from the yolk-sac boundary to where the cardiac crescent ends; E8.25 to E9.5 the
  heart-tube portion of the embryo with the yolk sac; E10.5 hearts dissected and sectioned whole. So E9.5 to E10.5
  mixes development with a change in what was sampled. Multiome E8.5 to E10.5: brain, spinal cord, posterior region
  and most yolk sac removed, ten embryos pooled per replicate; E10.5 five hearts per replicate.
- **Processing:** Cellpose 2 on DAPI; MERFISH cells with fewer than five genes dropped; expression log1p(1e4 c / sum c),
  checked on every released file. Knockouts and their wild-type controls were processed in separate runs; the two
  Gata4 replicates are independent littermates. Duplicated or interpolated records are not mentioned.
- **Staging:** timed mating (plug morning = E0.5), checked against EMAP reference morphology; no somite counts.
  `E8.25_late` is about E8.3, older than another E8.25 collection.
- **Use:** any use beyond the Challenge, publication or public presentation included, needs prior approval from
  Dr. Neil Chi and the data-generating team.

## T1 composition by chamber category (configs/t1_chamber_classes.yaml; measured 2026-11-01)

Share of each released T1 file's cells in the categories that the E12.5 test keeps (the heart chamber: cm,
endocardium, epicardium, vascular endothelium, blood) and in those it leaves out:

| File | Chamber | Left out (pharyngeal, surrounding, progenitor, pericardium) | Notable |
|---|---|---|---|
| E8.5_RNA | 0.31 | 0.69 | surrounding 0.40, progenitor 0.25 |
| E9.5_RNA | 0.48 | 0.52 | surrounding 0.26, progenitor 0.13 |
| E10.5_RNA | 0.69 | 0.31 | pharyngeal 0.12, surrounding 0.11; blood 0.17, endocardium 0.23, CM 0.14 |

The E10.5 truth of the T1 pseudo-task therefore keeps about a third of the tissue the E12.5 test leaves out, and its
blood share is far above B0's planned 2%.

