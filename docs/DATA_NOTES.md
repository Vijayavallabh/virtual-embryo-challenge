# Data notes (structure only)

From `vec data inventory` on the DGX, 2026-10-18 (`runs/inventory.json` there; checksums in `data/raw/SHA256SUMS`).
Structural facts only; no expression values. Re-run after each release (validation answers and test inputs arrive on
2026-10-20).

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

Every file satisfies `expm1(X).sum(1) == 1e4` (log1p CP10k). Stage cell counts are exactly 10x the board `ref_cells`
in `boards/index.json` (rounded), which confirms the server's 10% working copies.

## Quirks to handle in loaders

- **The z convention differs by stage:** median z step 10 (E6.75, E7.25, Mab21l2 KO), 15 (E8.0, E8.25_late), 14
  (E8.75); E9.5 has 98 z levels with no clean step. z is a section index scaled by an unknown per-stage factor, so RMS
  radii mix units across stages. The organisers say not to treat z as a measured distance until the Data Protocol says
  otherwise.
- **Mab21l2 KO:** non-unique obs names (make unique on load); 70 duplicate coordinates (jitter or drop before kNN);
  different cell-type vocabulary from the WT files (vCM, LPM, Paraxial Mesoderm, Branchial Arch...); dense X.
- **Vocabularies are not harmonised across stages** (E8.25_late and E8.75 share label definitions; E9.5 does not;
  `cm_celltype` exists only in E8.0, E8.25_late, E8.75 and E9.5, where "Unknown" means "no subtype given").
- All MERFISH files carry an integer `counts` layer; T1 counts must be recovered from `.X`.
- A layer stored under the literal key `None` appears in every file (an artefact; ignore it).
