# External public data for the Virtual Embryo Challenge (Human Team track): inventory and legality

Compiled 2026-10-18 for the NeurIPS 2026 Virtual Embryo Challenge (VEC), Track 1 (Human Team) only.
Scope: public mouse data that could serve T1 (scRNA, E12.5 extrapolation), T2-heart (3D MERFISH, E12.5 extrapolation; E8.5 interpolation on the val board), T2-embryo (3D MERFISH, E7.75 interpolation) and T3 (Mesp1-Cre conditional KOs at E8.75), judged against the challenge's external-data rules.

How this was built: the challenge site text and rules in `SCRATCH/site_text/`; primary repository listings (GEO FTP and GEO text records, ArrayExpress/BioStudies, CNGB/STOmicsDB, CELLxGENE curation API, JAX/Shendure download server, Marioni-lab servers, ENCODE API, Dryad API); PubMed/PMC full texts; GitHub (`gh api`). Only metadata was downloaded (all files below 15 MB, kept in `SCRATCH/dl_meta/`); no expression matrices were downloaded. Sizes come from HTTP HEAD requests or server directory listings.

Confidence labels used throughout:
- **[H]** verified today against a primary source (repository listing, sample sheet, paper text, API record).
- **[M]** from a secondary source or partial text, consistent with other evidence.
- **[L]** inferred or not verified. Treat as a lead, not a fact.

---

## 0. Bottom line

1. **For T1 and T2-heart, E12.5 can be bracketed instead of extrapolated.** E10.5 challenge data is released at P3 (2026-10-20) as training input; legal external data restarts after E13.5. Usable "late" anchors: JAX/Qiu 2024 embryos E13.75-P0, Feng 2022 hearts E14.5-P9, Feulner 2026 hearts E14.5-P7, MOSTA E14.5-E16.5 sections, Hu 2026 heart Visium E14.5-P2, ENCODE heart bulk RNA-seq E14.5-P0. Nothing in (E9.5, E13.5] may be used, by any route. [H for the data, H for the rule text]
2. **Lab stage labels for the same developmental state differ by up to about half a day.** "E8.5" means 3-5 somites (Sampath Kumar 2023), 6-7 (extended gastrulation atlas 2024), 0-12 (JAX pooled bin), 10-12 (Li 2016). "E9.5" means 15-18, 18-20, 24-26 or 26-29 somites depending on the lab. The rules say somite counts and Theiler stages are placed on the challenge scale but give no table. Section 2 gives a strict "any-convention" rule we recommend, plus a list of questions for the organisers. [H]
3. **T3 hazard found: GEO GSE165300** ([redacted 2026-10-19: B]) is scRNA-seq of *Mesp1-Cre; Wls fl/fl* embryos (Wnt secretion lost in the Mesp1 lineage) at E8.0, E8.5 and E9.5. This is very plausibly a phenocopy of the held-out *Mesp1-Cre; Ctnnb1 F/F* KO at a comparable stage. Exclude it (and any Ctnnb1 or Gata4 data, e.g. GSE208162, a single-cell Gata4 endocardial KO at E12.5). [H for the dataset; the phenocopy call is our judgement]
4. **Same-lab data hazard: GEO GSE176306** (Zhang ... Chi 2021, the lab that generated the challenge data) has Mesp1-lineage scRNA-seq at E7.25, E7.5, E7.75 and E8.25. Its E7.5 and E7.75 samples sit inside the T2-embryo window (E7.75 is the hidden test stage). Never use them for T2-embryo, nor in any shared model that produces T2-embryo predictions. [H]
5. **The CELLxGENE Census is not the main access route.** Mouse embryonic content is essentially JAX/Qiu 2024 plus the Sampath Kumar Slide-seq embryos and the Feulner perinatal hearts. Stages are coded only as Theiler stages (TS12-TS27), which are too coarse for the windows. The author day/somite columns are not in the Census. The one exact route is `donor_id` (equal to JAX `embryo_id`) joined to the JAX sample sheet. Pijuan-Sala 2019, MOCA, the extended atlas and de Soysa 2019 are not on CELLxGENE at all. [H]
6. **The rules conflict on T2/T3.** The Data page says "Task 1 is the only task where training on external public single-cell data is allowed", while Rules s.10 and the FAQ allow external data for every task and define T2-specific windows. No precedence clause exists. Get written confirmation before any T2/T3 submission uses external data. [H]
7. **Pretrained-model route.** The Navigo checkpoints (probably the organisers' `dynode_flow` reference; that identification is the brief's inference) are run on JAX-style data that includes E12.5-E13.5 days. Mouse foundation models trained on CELLxGENE or GEO corpora very likely saw E10-E13.5 cells. Both are prohibited for T1/T2-heart unless retrained on filtered data. The virtualembryo.ai atlas and predictor are built on the benchmark resource itself; do not use them as data. [M]

---

## 1. Governing rules (exact wording) and the open ambiguities

Rules s.10 (from `site_text/challenge_rules.txt`; the Data page and FAQ repeat the same text):
- "External public data, pre-trained models and published code may be used, provided their licences permit it and their use is described in the method summary. For Task 1 this is expressly encouraged."
- "Every external source must be disclosed with the submission. An undisclosed source is a rules violation whether or not it changed the result."
- "No measured data from a held-out stage or genotype may be used, by any route: not directly, not through a model pre-trained on it, and not through a public dataset that contains it."
- "Nearness counts, and the stage label is not the test ... so is a dataset labelled with a somite count or a Theiler stage that lands in the same window."
- Extrapolation (T1, T2-heart): "External data is excluded from after E9.5 up to and including E13.5. E10.5 and E12.5, the validation and test targets themselves, are excluded absolutely. Anything after E13.5 may be used, provided the source and its stages are stated explicitly with the submission." Exactly E9.5 is permitted; exactly E13.5 is not.
- Interpolation (T2): "external data may not be used if it stages strictly between the two released stages that bracket a held-out one. In Task 2 that means heart data strictly between E8.25 and E8.75, and embryo data strictly between E7.25 and E8.0. The bracketing stages themselves are permitted ... Fractional stagings do arise in external data, from somite counts and Theiler stages, and they are placed on this scale and judged by the same boundary; the boundary itself stays on the stages we published." (This restatement applies to submissions from 26 August 2026.)
- Genotype (T3): "Data from the held-out knockouts, from another allele of the same gene at a comparable stage, or from a perturbation that phenocopies them, is treated as the held-out condition."
- "A general-purpose resource that happens to span a held-out stage is not banned outright: the cells inside the window must be removed before training, and the removal described in the method summary."
- "If a source sits near the boundary and you are unsure, ask before submitting. The Organisers will answer, and an answer given in advance is not a finding of misconduct afterwards."

Ambiguities to settle with the organisers (all [H] as to the text; the resolutions are ours):

| # | Issue | Evidence | Working assumption until answered |
|---|---|---|---|
| A1 | Is external data allowed for T2 and T3? | Data page (Task 1 section): "Task 1 is the only task where training on external public single-cell data is allowed, provided the source is disclosed." Data page "Access" paragraph: "Task 1 additionally permits training on external public single-cell data". Versus Rules s.10 and FAQ: allowed for all tasks, with explicit T2 windows. No precedence clause in Rules or Terms. | Probably allowed (the T2 windows would be pointless otherwise), but confirm in writing before relying on it in a T2/T3 submission. |
| A2 | Which somite-to-E-day mapping is used? | None published by the organisers. Their own heart stage is called "E8.25_late". | Use the strict any-convention rule in Section 2.3. |
| A3 | Does the T2-heart E8.5 window still apply in P3, once E8.5 answers are released and no E8.5 board is scored? | The rule text is not phase-dependent. | Keep the exclusion. |
| A4 | Are E13.75 embryos (JAX, staged by limb morphometry, eMOSS) "after E13.5"? | Letter of the rule: yes. | Allowed with disclosure; prefer E14.0 and later if the margin matters. |
| A5 | How broad is "phenocopies"? | Not defined. | Exclude Wnt-pathway loss in mesoderm/heart and Gata4-partner losses; ask before using second-heart-field-loss mutants (Section 3.7). |
| A6 | Cross-species data | Not addressed. Human Carnegie stages CS13-CS20 roughly correspond to mouse E10.5-E13.5. | Avoid human embryonic data in that range for T1/T2-heart unless cleared. Human 12 pcw fetal heart is well past E13.5 by any mapping. |
| A7 | Shared models across tasks | Each task has its own held-out set. | A component used for task X must be trained only on data legal for X. The "global-safe" set (Section 2.4) satisfies every task at once. |

---

## 2. Placing external stagings on the challenge scale

### 2.1 Theiler stages (eMouseAtlas definitions) [H]

| TS | Typical dpc (range) | Somites | Note |
|---|---|---|---|
| TS9 | 6.5 (6.25-7.25) | - | pre/early streak |
| TS10 | 7.0 (6.5-7.75) | - | mid/late streak, early bud |
| TS11 | 7.5 (7.25-8.0) | 0 | OB, EB/LB, early head fold (EHF), late head fold (LHF) |
| TS12 | 8.0 (7.5-8.75) | 1-7 (12a 1-4; 12b 5-7) | heart starting to form |
| TS13 | 8.5 (8.0-9.25) | 8-12 | turning |
| TS14 | 9.0 (8.5-9.75) | 13-20 | |
| TS15 | 9.5 (9.0-10.5) | 21-29 | forelimb bud |
| TS16 | 10.0 (9.5-10.75) | 30-34 | |
| TS17 | 10.5 (10-11.25) | 35-39 | |
| TS18 | 11.0 (10.5-11.25) | 40-44 | |
| TS19 | 11.5 (11-12.25) | 45-47 | |
| TS20 | 12.0 (11.5-13) | 48-51 | |
| TS21 | 13.0 (12.5-14) | 52-55 | |
| TS22 | 14.0 (13.5-15) | 56-~60 | range overlaps the excluded E13.5 |
| TS23 / 24 / 25 / 26 / 27 | 15 / 16 / 17 / 18 / 19 (newborn) | - | |

Consequences: TS11 lies entirely inside the T2-embryo window. TS13 is the T2-heart window. TS15 straddles E9.5 versus E9.75. TS22 straddles E13.5 versus E14.0.

### 2.2 How labs actually label somite counts [H unless noted]

| Source | Label to somite counts |
|---|---|
| JAX, Qiu et al. 2024 (public sample sheet, 74 embryos) | "E8.5" bin = 0, 2-12 somites (E8.0-E8.5 pooled); E8.75 = 14-15; E9.0 = 16-18; E9.25 = 20-23; E9.5 = 24-26; E9.75 = 27-31; E10.0 = 32-34; morphometric (eMOSS) staging from E10.25 |
| Extended gastrulation atlas (Imaz-Rosshandler et al. 2024; per-cell `somite_count`) | E8.5 = 6, 7; E8.75 = 8, 10, 12; E9.0 = 13-16; E9.25 = 21-24; E9.5 = 18, 20 (paper: E9.5 embryos "smaller than expected") |
| Sampath Kumar et al. 2023 (Slide-seq) | E8.5 = 3-5; E9.0 = 10-12; E9.5 = 15-18 |
| de Soysa et al. 2019 | E7.75 = 4; E8.25 = 8; E9.25 = 21 |
| Ibarra-Soria et al. 2018 | E8.25 = 4 (TS12) |
| Li et al. 2016 | E8.5 = 10-12; E9.5 = 26-29; E10.5 = 35-40 |
| Lohoff et al. 2022 (seqFISH) | 8-12 somites, described as E8.5-E8.75 |
| Pijuan-Sala et al. 2019 (authors' metadata script) | E6.5-E6.75 = TS9; E7.0-E7.25 = TS10; E7.5-E7.75 = TS11; E8.0-E8.5 = TS12; "mixed_gastrulation" = TS9-10 |

### 2.3 Recommended removal predicates (strict by default)

An embryo or sample is **in-window (remove)** if any of its own label, its somite count or its Theiler stage places it inside the window.

- **T2-heart interpolation window (E8.25, E8.75)**: remove if the label is strictly between E8.25 and E8.75; or somites are 6-13; or the stage is TS13; or the sample is an unstaged pool labelled E8.5.
  - Why 6-13: TS13 = 8-12 somites; the extended atlas calls 6-7 somites "E8.5"; linear interpolation between JAX anchors (about 1 somite at E8.0, 14-15 at E8.75) puts 6-13 somites strictly inside. That interpolation is our derivation, not a published mapping. 4-5 somites stay allowed unless labelled inside (TS12a = E8.0; de Soysa labels 4 somites E7.75; Ibarra-Soria labels them E8.25). 14+ somites stay allowed (JAX labels 14-15 as exactly E8.75).
  - Moderate alternative, only if the organisers agree: label inside, or 8-12 somites.
- **T1 and T2-heart extrapolation window (E9.5, E13.5]**: remove if the label is in (E9.5, E13.5]; or somites are 27 or more (JAX's E9.75 bin starts at 27); or the stage is TS16-TS21; or the stage is TS22 without a day label of E13.75 or later. Keep labels of E9.5 or earlier with 26 or fewer somites (or unknown somites), and labels of E13.75 or later.
  - Moderate alternative: drop only 30+ somites (TS16), keeping all TS15 embryos labelled E9.5.
- **T2-embryo interpolation window (E7.25, E8.0)**: remove if the label is strictly between E7.25 and E8.0; or the stage is TS11 (OB, EB, LB, EHF or LHF morphology); or the embryo has 0 somites; or the sample is an unstaged pooled gastrulation sample (e.g. Pijuan-Sala `mixed_gastrulation`). Keep E7.25 or earlier, and E8.0 or later with at least 1 somite.
- **T3**: stage is irrelevant; genotype rules in Section 3.7.

### 2.4 "Global-safe" set for a shared backbone (Generality Award designs)

Apply all three stage predicates plus the T3 genotype exclusions. Data surviving this filter is legal for every task, so one pretrained trunk can serve T1, T2 and T3 without per-task bookkeeping. Task-specific heads can then add back data that is legal for their own task only (e.g. the E7.5-E7.75 gastrulation cells for a T1-only or T2-heart-only component; the E8.5 cells for T1 and T2-embryo).

### 2.5 Worked boundary cases

| Case | T1 | T2-heart | T2-embryo | Reason |
|---|---|---|---|---|
| JAX embryo_1 (0 somites, in the "E8.5" bin) | keep | remove | remove | label E8.5 is in the heart window; 0 somites is head-fold (TS11) |
| JAX embryos 13-14 (14-15 somites, E8.75) | keep | keep | keep | bracket stage, exactly E8.75 |
| JAX embryos 22-24 (24-26 somites, E9.5) | keep | keep | keep | exactly E9.5 |
| JAX embryos 25-29 (27-31 somites, E9.75) | remove | remove | keep | after E9.5 |
| JAX embryos 50-52 (E13.75) | keep (disclose) | keep (disclose) | keep | after E13.5 |
| Extended atlas E8.75 (8-12 somites) | keep | remove | keep | label is a bracket but TS13 is inside |
| Extended atlas E9.0, 13 somites | keep | remove (strict) | keep | interpolated about E8.67 |
| Extended atlas E9.5 (18-20 somites) | keep | keep | keep | younger than E9.5 by any convention |
| de Soysa E8.25 (8 somites) | keep | remove (strict); ask | keep | TS13 |
| de Soysa E7.75 (4 somites) | keep | keep | remove | label inside the embryo window |
| Sampath Kumar "E8.5" (3-5 somites) | keep | remove | keep | label inside the heart window; at least 1 somite |
| Sampath Kumar "E9.0" (10-12 somites) | keep | remove | keep | TS13 |
| Sampath Kumar "E9.5" (15-18 somites) | keep | keep | keep | TS14 |
| Li 2016 E9.5 (26-29 somites) | strict: remove (27-29 = JAX E9.75); moderate: keep | same as T1 | keep | ask the organisers |
| Ibarra-Soria E8.25 (4 somites) | keep | keep | keep | TS12a, label E8.25 |
| Pijuan-Sala E8.25 (pooled, TS12) | keep | keep (label-based) | keep | somites per embryo unknown |

---

## 3. Dataset cards

Field order in each card: citation; modality and scope; stages (with conversion); cells; genes; access, formats and sizes; licence and terms; legal use and removal per task; expected value.

### 3.1 Gastrulation and early organogenesis (dissociated)

#### 3.1.1 Pijuan-Sala et al. 2019, mouse gastrulation atlas (with WT and Tal1-/- chimeras)
- **Citation**: Pijuan-Sala B, et al. Nature 566:490-495 (2019). DOI 10.1038/s41586-019-0933-9. [H]
- **Modality / scope**: droplet scRNA-seq (10x Genomics; the MouseGastrulationData vignette says v1 chemistry [M]); whole embryos, pooled per stage, extra-embryonic tissue included.
- **Stages**: E6.5, E6.75, E7.0, E7.25, E7.5, E7.75, E8.0, E8.25, E8.5, plus `mixed_gastrulation` (two pools of several stages, coded TS9-10). Authors' Theiler mapping: E8.0-E8.5 = TS12, i.e. 1-7 somites, so their "E8.5" may be younger than other labs' E8.5. Staged by Downs-Davies morphology to exclude outliers. [H]
- **Cells**: 116,312 [H]. Per stage, from the extended-atlas cell browser [H]: E6.5 3,484; E6.75 2,075; E7.0 14,749; E7.25 13,537; E7.5 10,994; E7.75 14,493; E8.0 16,681; E8.25 15,935; E8.5 16,909; mixed 7,455. Chimeras: WT (10 samples, E7.5 and E8.5) and Tal1-/- (4 samples, E8.5).
- **Genes**: about 29.5k (29,452 in the vignette example) [M].
- **Access**: ArrayExpress E-MTAB-6967 (atlas), E-MTAB-7324 (WT chimeras), E-MTAB-7325 (Tal1-/- chimeras), E-MTAB-6970 (Smart-seq2 endothelium) [H]. Processed: `https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz` (1.59 GB), `chimera_tal1_data.tar.gz` (573 MB), `chimera_wt_data.tar.gz` (197 MB) [H]. Bioconductor `MouseGastrulationData` (`EmbryoAtlasData`, `WTChimeraData`, `Tal1ChimeraData`; release version 1.24.0 [M]).
- **Licence / terms**: EMBL-EBI terms ("imposes no additional restriction ... than those provided by the data owner"; attribution expected) [H]; R package GPL-3 [H].
- **Legal use and removal**:
  - T1: all cells.
  - T2-heart: remove `stage == E8.5`. E8.25 kept (label-based).
  - T2-embryo: remove E7.5, E7.75 and `mixed_gastrulation`.
  - T3: all; the Tal1-/- chimeras are a legal perturbation (Tier C, Section 3.7).
- **Value**: high for T2-embryo (transcriptome at E6.5-E7.25 and E8.0-E8.5 around the hidden E7.75); medium for T1 and T2-heart (early heart progenitors; mostly superseded by 3.1.2, which contains these cells).

#### 3.1.2 Extended mouse gastrulation and organogenesis atlas (Imaz-Rosshandler et al. 2024)
- **Citation**: Imaz-Rosshandler I, Rode C, Guibentif C, et al. Development 151:dev201867 (2024). DOI 10.1242/dev.201867. [H]
- **Modality / scope**: 10x scRNA-seq (chemistry not verified [L]) of whole embryos, some split into anatomical pieces (`anatomy`: EP, YS, anterior/medial/posterior section, ...); merged with Pijuan-Sala 2019 (`embryo_version` = Original / Extension).
- **Stages**: 13 time points E6.5-E9.5 at 6 h spacing; new sampling at E8.5, E8.75, E9.0, E9.25, E9.5. Per-cell `stage` and `somite_count` (from the UCSC cell browser binaries, cross-tabulated today) [H]:

  | stage | somites (cells) |
  |---|---|
  | E8.5 | pooled original 16,909; 6 (24,694); 7 (16,469) |
  | E8.75 | 8 (26,578); 10 (11,001); 12 (31,331) |
  | E9.0 | 13 (28,623); 14 (34,933); 15 (14,788); 16 (14,676) |
  | E9.25 | 21 (16,547); 22 (17,768); 23 (19,443); 24 (17,545) |
  | E9.5 | 18 (22,354); 20 (17,277) |

- **Cells**: 430,339 total (116,312 original plus 314,027 new); 88 cell states including "Cardiopharyngeal progenitors SHF", "Lateral plate mesoderm", "Pharyngeal mesoderm" and cardiomyocyte states [H].
- **Access**: `https://marionilab.github.io/ExtendedMouseAtlas/`. One tarball, `https://bioinformatics.stemcells.cam.ac.uk/rlh60/Supplemental/ExtendedMouseAtlas/ExtendedMouseAtlas.tar.gz` (25 GB), containing `embryo_complete.h5ad`, `embryo_sce.rds`, `embryo_counts.tar.gz`, `embryo.zarr.zip`, `metadata_cells.csv` and others [H]. ArrayExpress E-MTAB-11763 (10x) and E-MTAB-13645 (Smart-seq2) [H]. Cell browser `https://cells-test.gi.ucsc.edu/?ds=ext-mouse-atlas` [H].
- **Licence / terms**: none stated on the site; ArrayExpress terms apply to raw data [M].
- **Legal use and removal**:
  - T1: all cells (labels at most E9.5; at most 24 somites).
  - T2-heart: remove `stage in {E8.5, E8.75}`; strict also removes `somite_count == 13`; also remove `Mixed gastrulation`. Keep E6.5-E8.25, E9.0 at 14-16 somites, E9.25 and E9.5.
  - T2-embryo: remove E7.5, E7.75 and `Mixed gastrulation`.
  - T3: all.
- **Value**: very high for T1 (dense, somite-resolved E8.5-E9.5 ladder with second-heart-field and cardiomyocyte states; directly supports "supervision across the whole stage ladder"). High for T2-heart (E9.0-E9.5 programmes on the 500-gene panel) and T2-embryo (brackets).

#### 3.1.3 Single-embryo time-resolved gastrulation (Mittnenzweig et al. 2021) and the related Weizmann series
- **Citation**: Mittnenzweig M, et al. Cell 184:2825-2842 (2021). DOI 10.1016/j.cell.2021.04.004. [H]
- **Modality / scope**: MARS-seq of 153 individually sampled whole embryos; each embryo has an inferred time stamp (Et6.5-Et8.1) from morphology, transcriptional ranking and embryo size [H]. GEO lists strains C57BL/6, ICR, and chimera/tetraploid hosts [H].
- **Cells**: about 33,700-33,900 [H].
- **Perturbations**: Foxc1/Foxc2 double-KO chimeras and tetraploid embryos [H].
- **Access**: GEO GSE169210 (`RAW.tar` 159 MB plus per-condition cell metadata, 16-270 KB each) [H]; code `github.com/tanaylab/embflow`; viewer `tanaylab.weizmann.ac.il/embflow` (responds today).
- **Legal use and removal**:
  - T2-embryo: keep only embryos with inferred time at or below E7.25 or at or above E8.0. Most of the series falls inside the window, so little remains.
  - T1, T2-heart, T3: all (Foxc1/2 DKO is Tier C).
- **Value**: low for T2-embryo after removal; low elsewhere.
- **Related, same labs (not inspected in detail)**: Mayshar et al. 2023 Cell (rabbit and mouse; DOI 10.1016/j.cell.2023.04.037); Cheng et al. 2022 Cell (Tet TKO; DOI 10.1016/j.cell.2022.06.049); Mukamel et al. 2022 NSMB (Dnmt3a/b; DOI 10.1038/s41594-022-00885-6); Hadas et al. 2024 Nature (Bmp4; DOI 10.1038/s41586-024-07937-5). Same rules apply: remove embryos with Et in (7.25, 8.0) for T2-embryo. These are legal T3 perturbations (Tier C). [M]

#### 3.1.4 Grosswendt et al. 2020, zygotic epigenetic-regulator knockouts
- **Citation**: Grosswendt S, et al. Nature 584:102-108 (2020). DOI 10.1038/s41586-020-2552-x. [H]
- **Content**: 10 KOs (Dnmt1, Dnmt3a, Dnmt3b, G9a, Rnf2, Eed, Kdm2b, L3mbtl2, Kmt2a, Kmt2b) at E8.5, 8-12 embryos each (7,548-25,408 cells per KO); Eed KO also at E6.5 and E7.5; WT reference of 88,779 cells from 50 embryos E6.5-E8.5 [H]. GEO also lists Cdkn2a and Eed/Cdkn2a samples [H].
- **Access**: GEO GSE137337 (per-KO tar.gz 82-593 MB; RAW.tar 6.0 GB) [H].
- **Legal use and removal**:
  - T3: all (no Wnt/Gata genes; Tier C).
  - T1: all.
  - T2-heart: remove E8.5.
  - T2-embryo: remove E7.5 (and any E7.75).
- **Value**: medium for T3 (many KO-vs-WT contrasts at nearly the right stage, for learning response magnitude and direction priors); low elsewhere.

#### 3.1.5 Ibarra-Soria et al. 2018 (E8.25 whole embryos)
- **Citation**: Nat Cell Biol 20:127-134 (2018). DOI 10.1038/s41556-017-0013-z. [H]
- **Content**: 10x, 4-somite embryos (TS12), 19,396 cells from 23 embryos; replication about 7,000 cells (C57BL/6 x CBA) [H].
- **Access**: ArrayExpress E-MTAB-6153 (and E-MTAB-5728, E-MTAB-6155) [H].
- **Legal use**: all tasks. The label E8.25 and 4 somites (TS12a) are both outside the heart window.
- **Value**: low-medium (a clean E8.25 bracket for T2-heart).

#### 3.1.6 Lower-priority gastrulation sets
- Nowotschin et al. 2019 Nature (endoderm; 112,217 cells, E3.5-E8.75; GEO GSE123046; DOI 10.1038/s41586-019-1127-1) [H]. For T2-embryo remove E7.5; for T2-heart, E8.75 is a bracket and allowed.
- Argelaguet et al. 2019 Nature (scNMT-seq, E4.5-E7.5; GEO GSE121708; DOI 10.1038/s41586-019-1825-8) [H]. For T2-embryo remove E7.5.
- Argelaguet multiome (MouseGastrulationData `RAMultiomeData`; 11 samples, WT and T (Brachyury)-KO; stages not verified [L]).
- Guibentif et al. 2021 T-/- chimeras (MouseGastrulationData `TChimeraData`; E7.5 and E8.5 pools) [H]. For T2-embryo remove E7.5; for T2-heart remove E8.5. T3 Tier C (Wnt-adjacent; disclose).

### 3.2 Cardiac progenitors and the early heart

#### 3.2.1 Lescroart et al. 2018 (Mesp1+ progenitors)
- **Citation**: Science 359:1177 (2018). DOI 10.1126/science.aao4174. [H]
- **Content**: 598 single cells [H]: E6.75 Mesp1-GFP+ WT 172, E6.75 Mesp1-KO 85, E7.25 WT 341 (plate-based full-length scRNA-seq [M]).
- **Access**: GEO GSE100471 (`HTSeq_counts.txt.gz`, 9.0 MB) [H].
- **Legal use**: all tasks; E6.75 and E7.25 are released or bracket stages for T2-embryo. Mesp1-KO is Tier C for T3.
- **Value**: low-medium; small, but the closest match to the T3 lineage definition at origin.

#### 3.2.2 Zhang et al. 2021, Chi lab Mesp1-lineage scRNA-seq (same lab as the challenge data)
- **Citation**: Zhang Q, Carlin D, Zhu F, ..., Chi NC. Circ Res 129:474-487 (2021). DOI 10.1161/CIRCRESAHA.121.318943; bioRxiv 10.1101/2021.01.08.425950. [H]
- **Content**: 10x v2 scRNA-seq of *Mesp1-Cre; Rosa26-tdT* embryos, staged by morphology [H]: E7.25 (1 sample), E7.5 (4), E7.75 (2), E8.25 (3; yolk sac removed in one). Cell count not verified [L].
- **Access**: GEO GSE176306 (`combined.Rdata.gz`, 4.1 GB; Seurat 3.1.5 objects) [H].
- **Legal use and removal**:
  - T2-embryo: remove E7.5 and E7.75. Keep E7.25 and E8.25.
  - T2-heart: all kept. E8.25 is the bracket; E7.25-E7.75 are earlier than the window.
  - T1: all.
  - T3: all (WT Mesp1 lineage).
- **Caution**: same lab as the challenge data; ask whether these embryos overlap the withheld E7.75 scRNA. Keep the E7.5/E7.75 cells out of any shared model used for T2-embryo.
- **Value**: high for T3 (WT reference of exactly the Mesp1-Cre lineage, same lab and protocol family); medium for T2-heart and T1.

#### 3.2.3 Tyser et al. 2021 (cardiac crescent to linear heart tube)
- **Citation**: Science 371:eabb2986 (2021). DOI 10.1126/science.abb2986. [H]
- **Content**: Smart-seq2; 3,105 QC cells from the micro-dissected anterior cardiac region over "six defined stages spanning approximately 12 hours (roughly E7.75 to E8.25)": stages -1 to 3 (crescent; stage 3 about E8.0) and the linear heart tube (LHT); plus a small 4-region mesoderm reference at stage 3. No somite counts [H].
- **Access**: ArrayExpress E-MTAB-7403; ENA PRJEB14363 [H]. Processed `.Rds` at `https://content.cruk.cam.ac.uk/jmlab/mouseEmbryonicHeartAtlas/` (7-130 MB each) [H]. The web atlas returned HTTP 503 today.
- **Legal use and removal**:
  - T1: all.
  - T2-heart: keep the crescent stages; LHT is borderline (about E8.25, somites unknown), so strict removes it.
  - T2-embryo: remove the crescent stages (inside (E7.25, E8.0)); only LHT may remain, and it adds little.
  - T3: all.
- **Value**: low-medium (early heart-field states, small n).

#### 3.2.4 de Soysa et al. 2019 (cardiogenic region, WT and Hand2-null)
- **Citation**: Nature 572:120-124 (2019). DOI 10.1038/s41586-019-1414-x. [H]
- **Content**: 10x v2. "The entire cardiogenic region" plus surrounding tissue, with pharyngeal arches at E9.25. E7.75 = 4 somites, E8.25 = 8 somites, E9.25 = 21 somites. 36,777 WT cells and 37,149 Hand2-null cells [H].
- **Access**: GEO GSE126128: `AllTimePoints_WT_10X.csv.gz` (1.3 GB), `E775E825_WTKO_10X.csv.gz` (1.6 GB), `E925_WTKO_10X.csv.gz` (1.0 GB); browser `mouse-cardiac.cells.ucsc.edu` (responds) [H].
- **Legal use and removal**:
  - T1: all WT stages.
  - T2-heart: keep E7.75 and E9.25. E8.25 (8 somites = TS13) is removed under the strict rule; allowed only on a label-only reading, so ask.
  - T2-embryo: remove E7.75 (label inside the window); keep E8.25 and E9.25.
  - T3: WT kept; Hand2-null is Tier B (possible phenocopy of right ventricle/outflow tract loss; exclude unless cleared).
- **Value**: very high for T1 (dissection closest to T1's heart-centred dissection; E9.25 sits just before E9.5); high for T2-heart (E9.25).

#### 3.2.5 Li et al. 2016 (anatomically zoned single heart cells)
- **Citation**: Dev Cell 39:491-507 (2016). DOI 10.1016/j.devcel.2016.10.014. [H]
- **Content**: Fluidigm C1 SMARTer. E8.5 (10-12 somites) 118 cells; E9.5 (26-29 somites) 949; E10.5 (35-40 somites) 1,166; manual dissection zones (2/7/9); plus Nkx2-5-/- (Smart-seq2) and Isl1-lineage cells; CD1 mice [H].
- **Access**: GEO GSE76118 (`RAW.tar` 426 MB) [H].
- **Legal use and removal**:
  - T1 and T2-heart: remove E10.5. E9.5 is strict-remove (27-29 somites fall in the JAX E9.75 bin) or moderate-keep (label E9.5, TS15).
  - T2-heart: also remove E8.5.
  - T2-embryo: all.
  - T3: Nkx2-5-/- is Tier C/B; ask.
- **Value**: low.

#### 3.2.6 Jia et al. 2018 (Nkx2-5 and Isl1 reporter progenitors)
- **Citation**: Jia G, et al. Nat Commun 9:4877 (2018). DOI 10.1038/s41467-018-07307-6; PMC6242939. [H]
- **Content**: Fluidigm C1; 167 Nkx2-5+ and 254 Isl1+ cells, E7.5-E9.5; Isl1 KO at E9.5; scATAC-seq [H].
- **Access**: ENA PRJEB23303 [H].
- **Legal use**:
  - T1: all.
  - T2-heart: remove E8.5.
  - T2-embryo: remove E7.5.
  - T3: Isl1 KO is Tier B.
- **Value**: low.

### 3.3 Whole-embryo organogenesis to birth

#### 3.3.1 JAX single-cell time-lapse, E8 to P0 (Qiu et al. 2024 Nature)
- **Citation**: Qiu C, Martin BK, Welsh IC, et al. Nature 626:1084-1093 (2024). DOI 10.1038/s41586-024-07069-w. [H]
- **Modality / scope**: sci-RNA-seq3 (nuclei), whole embryos. 12.4M nuclei from 83 embryos overall; 11,441,407 nuclei from 74 embryos after QC (E8-P0); about 154,614 nuclei per embryo on average [H]. Additional embryos 76-84 (birth series, P0 natural birth and C-section) and 85-96 (somitogenesis validation, labelled E8.5-E9.25) are in GEO but not in the 74-embryo sample sheet [H].
- **Stages**: 43 distinct day labels across the 74 sheet embryos [H].
  - Somite-staged E8-E10: embryos 1-12 = 0, 2-12 somites (bin "E8.5"); 13-14 = 14-15 (E8.75); 15-17 = 16-18 (E9.0); 18-21 = 20-23 (E9.25); 22-24 = 24-26 (E9.5); 25-29 = 27-31 (E9.75); 30-32 = 32-34 (E10.0).
  - eMOSS-staged: 33-49 = E10.25-E13.5.
  - 50-52 = E13.75; 53 = E14.0; 54-56 = E14.25-E14.75; 57-72 = E15.0-E18.75 at 6 h steps; 73-74 = P0.
- **Nuclei per embryo (GEO annotation files, counted today)** [H]: embryos 1-12 have 1,679-24,551 nuclei each (153,597 total). Later embryos have 33k-480k per sequencing run (e.g. E13.75: 33,187 and 49,200; E14.0: 56,650; P0 embryo_74: 1.61M in run 25). Counts for embryos 13-32 were not checked (other runs).
- **Genes**: GENCODE vM12 (mm10); 45,525 features in the CELLxGENE version [H].
- **Access** [H]:
  - JAX server `https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/download/`: `adata/adata_JAX_dataset_{1..4}.h5ad` (11, 8.9, 9.9, 11 GB); `adata/df_cell.csv` (1.5 GB; columns `cell_id, day, embryo_id, experimental_batch, major_trajectory, celltype_update`); `mtx/` (61 GB per run); `embryo_exon_intron/` (per-embryo RDS, for velocity); `other/exp.tar.gz` (926 MB) and `exp2.tar.gz` (626 MB): per-cell-type mean expression by 6-h day and by somite count; `other/sample_sheet.txt` (embryo, day, somite count).
  - GEO GSE228590 (per-run mtx 0.16-5.7 GB plus annotation CSVs) and GSE186069 (an earlier E8.5 pilot of 0-12-somite embryos; 2.1 GB mtx).
  - CELLxGENE collection `45d5d2c3-bc28-4814-aed6-0bb6f0e11c82` (published 2024-03-11): raw counts, 11,441,407 cells (H5AD 46.3 GB, `https://datasets.cellxgene.cziscience.com/1bc30289-9565-4099-abf9-3326328c11ac.h5ad`); four normalized subsets of about 34 GB each; per-cluster subsets such as Cardiomyocytes (37,205 cells, 0.56 GB, `.../b5bab198-b6bf-40e9-813f-7d9692c49ef0.h5ad`), Endothelium (312,029 cells, 3.9 GB), Mesoderm (3,267,338 cells, 40.0 GB), Hepatocytes (289,569 cells, 3.5 GB). On CELLxGENE, `donor_id` = `embryo_id` and stages are coded TS12-TS27.
- **Licence / terms**: CELLxGENE data are CC BY 4.0 (AWS Open Data registry entry for the Census) [M]; the JAX server states no licence [H]; GEO: "NCBI places no restrictions on the use or distribution of the GEO data", though submitters may claim rights [M; the page was captcha-blocked today, text confirmed via search snippet].
- **Legal use and removal** (computed from the sample sheet with the strict predicates):
  - T1: keep embryos 1-24 and 50-74 (49 embryos); remove 25-49 (E9.75-E13.5). Validation embryos 85-96 are allowed (labels at most E9.25). Optional margin: start the late block at embryo 53 (E14.0).
  - T2-heart: keep embryos 13-24 and 50-74 (37 embryos); remove 1-12 (labelled E8.5) and 25-49; remove validation embryos 85-91 (labelled E8.5).
  - T2-embryo: keep embryos 2-74; remove embryo_1 (0 somites).
  - T3: all.
- **Value**: very high for T1 (brackets E12.5 on both sides in one technology; somite-resolved E8.0-E9.5; heart-region trajectories). High for T2-heart (expression programmes on the panel genes; no spatial data). Medium for T3 (WT programmes of Mesp1-lineage derivatives).

#### 3.3.2 MOCA, mouse organogenesis cell atlas (Cao et al. 2019)
- **Citation**: Nature 566:496-502 (2019). DOI 10.1038/s41586-019-0969-x. [H]
- **Content**: sci-RNA-seq3; about 2 million cells from 61 embryos at E9.5, E10.5, E11.5, E12.5 and E13.5, staged by gestational day only (no somite counts) [H].
- **Access**: GEO GSE119945: `gene_count.txt.gz` (3.7 GB); `cell_annotate.csv.gz` (101 MB; columns `sample, id, sex, day, ..., Main_Cluster, Sub_Cluster, detected_doublet`); `gene_annotate.csv.gz` [H].
- **Legal use and removal**:
  - T1 and T2-heart: keep `day == 9.5` only. It is moderate-risk, since actual somite stages are unknown.
  - T2-embryo and T3: all, but of low relevance.
- **Value**: low-medium; superseded by the somite-staged JAX E9.25-E9.5 embryos.

#### 3.3.3 Huang et al. 2023 (22 mutants at E13.5)
- **Citation**: Nature 623:772-781 (2023). DOI 10.1038/s41586-023-06548-w. [H]
- **Content**: sci-RNA-seq3; 101 embryos at E13.5, 22 mutant and 4 WT strains; about 1.63-1.67M nuclei. Mutants include Ttc21b, Carm1, Gli2, Sox9 regulatory alleles, Scn11a, Ror2, Gorab, Cdkl5, Scn10a/11a, Atp6v0a2, Fat1-TAD, and regulatory deletions near Smad3, Tbx5, Neurog2, Sim1, Smad7, Dmrt1, Tbx3, Twist1 and the ZRS. No Ctnnb1, Gata4 or Gata6; Ror2 is a non-canonical Wnt receptor [H]. Browser `atlas.gs.washington.edu/mmca_v2/`. GEO accession not verified [L].
- **Legal use**: T1 and T2-heart: none (E13.5 is excluded). T2-embryo: legal but irrelevant. T3: legal (Tier C).
- **Value**: low-medium for T3 (many matched KO/WT contrasts, but far from E8.75).

### 3.4 Spatial embryo atlases (2D and 3D)

#### 3.4.1 Spateo 3D Stereo-seq embryos and hearts (Qiu et al. 2024 Cell; also "MOSTA3D")
- **Citation**: Qiu X, Zhu DY, Lu Y, et al. "Spatiotemporal modeling of molecular holograms." Cell 187:7351-7373 (2024). DOI 10.1016/j.cell.2024.10.011 (erratum 10.1016/j.cell.2025.02.026). Companion preprint: Cheng M, et al. bioRxiv 10.1101/2024.08.17.608366 (MOSTA3D; 94 slices at E9.5 and 91 at E11.5). [H]
- **Content**: 3D whole-embryo Stereo-seq (serial sections) at E9.5 and E11.5, "capturing eight million cells"; whole transcriptome [H]. The data availability text could not be read (publisher returned 403).
- **Files** [H]: `https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/mouse_E9.5_embryo.h5ad` (822 MB), `mouse_E9.5_heart.h5ad` (35.8 MB), `mouse_E11.5_embryo.h5ad` (6.69 GB), `mouse_E11.5_heart.h5ad` (138 MB), plus `E11.5_embryo_heart_3D.h5ad` (222 MB, no link). The organisers' lab GitHub (`aristoteleo/PantheonOS`) describes `mouse_E9.5_embryo.h5ad` as "a 3D Stereo-seq dataset", and the Spateo datasets page lists the same four files.
- **Licence / terms**: CNGB terms (the page did not render its text today) [L]; the Spateo code is BSD-2 [H].
- **Legal use and removal**:
  - T1 and T2-heart: E9.5 files only; remove every E11.5 file.
  - T2-embryo and T3: legal, low relevance.
- **Value**: high for T2-heart. It is the only public whole-transcriptome 3D mouse heart at a released stage (E9.5), useful for spatial priors (neighbourhood structure, shape, gene-gene structure on the 500-gene panel) and for calibrating 3D growth. Check its physical units before using it for scale.

#### 3.4.2 MOSTA, Stereo-seq sagittal sections E9.5-E16.5 (Chen et al. 2022)
- **Citation**: Cell 185:1777-1792 (2022). DOI 10.1016/j.cell.2022.04.003. [H]
- **Content**: 53 sagittal sections by the download page count [H]: E9.5 (5), E10.5 (4), E11.5 (4), E12.5 (6), E13.5 (4), E14.5 (7), E15.5 (5), E16.5 (18). Bin50 (about 25 µm; [M]) and some cell-bin h5ad files; whole transcriptome.
- **Access**: `https://db.cngb.org/stomics/mosta/download/` (STOmicsDB STDS0000058). URL pattern `https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/{STAGE}_{SECTION}.MOSTA.h5ad`. Sizes: E9.5 0.29-3.8 GB; E14.5 2.1-5.9 GB (E14.5_E1S1 = 5.93 GB); E15.5 3.0-4.6 GB; E16.5 2.1-5.3 GB (E16.5_E1S1 = 4.54 GB). `Mouse_embryo_all_stage.h5ad` (7.2 GB) mixes all stages. GEM and barcode files are also offered [H].
- **Licence / terms**: CNGB terms (not readable today) [L].
- **Legal use and removal**:
  - T1 and T2-heart: keep E9.5, E14.5, E15.5 and E16.5 sections; remove E10.5-E13.5. Filter the all-stage file by stage, or do not use it.
  - T2-embryo and T3: all, low relevance.
- **Value**: medium for T1/T2-heart as a late bracket. Heart regions at E14.5-E16.5 give expression plus physical size in 2D. Heart annotation in the h5ad is likely but not verified [L].

#### 3.4.3 Sampath Kumar et al. 2023, Slide-seq whole embryos with sc3D reconstructions
- **Citation**: Nat Genet 55:1176-1185 (2023). DOI 10.1038/s41588-023-01435-6. [H]
- **Content** [H]: Slide-seqV2.
  - E8.5: 2 embryos at 3-5 somites; 15 and 17 sagittal sections of 10 µm about 30 µm apart.
  - E9.0: 10-12 somites; 26 sections about 20 µm apart.
  - E9.5: 15-18 somites; partial (13 midline slices across 3 embryos).
  - Tbx6 KO at E9.5: trunk transverse sections, 1 WT and 2 KO.
  - 533,116 beads in total (median 1,798 UMIs and 1,224 genes per bead). The 3D virtual embryos include the heart.
- **Access**: GEO GSE197353 (`RAW.tar` 1.56 GB) [H]; figshare 3D inputs (E8.5 embryo 1 `figshare.com/s/1c29d867bc8b90d754d2`, E8.5 embryo 2 `21695849`, E9.0 `21695879`) [H]; CELLxGENE collection `d74b6979-efba-47cd-990a-9d80ccf29055` (75 datasets: one per puck, plus "All integrated" with 533,116 beads at 2.77 GB, `https://datasets.cellxgene.cziscience.com/e38f5610-2487-448e-a0aa-f5dc7c3b5850.h5ad`, and a Tbx6 WT/KO set of 499 beads). Every dataset is coded TS14, which is inconsistent with the paper's stages [H].
- **Legal use and removal**:
  - T1: all.
  - T2-heart: remove E8.5 (label) and E9.0 (10-12 somites = TS13); E9.5 (15-18 somites) is legal but partial.
  - T2-embryo: all (at least 3 somites; at or after E8.0 by somites).
  - T3: WT kept; Tbx6 KO is Tier B ([redacted 2026-10-19: B]; ask).
- **Value**: medium-high for T2-embryo. It is the only public 3D whole-embryo transcriptomic reconstruction just after the E8.0 bracket (3-5-somite "E8.5" embryos are about E8.0-E8.25 by somites). Low for T2-heart.

#### 3.4.4 Lohoff et al. 2022, seqFISH of 8-12-somite embryos
- **Citation**: Nat Biotechnol 40:74-85 (2022). DOI 10.1038/s41587-021-01006-2. [H]
- **Content**: 3 embryos at 8-12 somites (E8.5-E8.75), 20 µm sagittal sections, 387 targeted genes (351 barcoded seqFISH plus 36 smFISH), 57,536 cells; transcriptome-wide imputation available. The heart tube is captured [H].
- **Access**: Bioconductor `MouseGastrulationData::LohoffSeqFISHData` (observed, imputed and molecule positions); code `github.com/MarioniLab/SpatialMouseAtlas2020` [H].
- **Legal use**: T2-heart: none (TS13). T1, T2-embryo and T3: all.
- **Value**: low-medium (2D; after the embryo bracket).

#### 3.4.5 Harland et al. 2025, seqFISH spatiotemporal gastrulation atlas
- **Citation**: Cell Rep 44:116047 (2025). DOI 10.1016/j.celrep.2025.116047. [H]
- **Content**: 4 new embryos (2 labelled E7.25 at late streak, 2 at E7.5), 20 optical sections, 351-gene seqFISH (387-gene panel), 14,794 cells; integrated with the Lohoff E8.5 data and the E6.5-E9.5 scRNA atlas (more than 150k cells, 82 types) [H].
- **Access**: images Zenodo 13977985; processed data and portal `shiny.maths.usyd.edu.au/SpatiotemporalMouseAtlas/` (timed out today); code `github.com/ltgharland/Spatiotemporal-Atlas-of-Mouse-Gastrulation` [H].
- **Legal use**: T2-embryo: keep the E7.25 embryos and remove the E7.5 embryos. T1, T2-heart and T3: all.
- **Value**: medium for T2-embryo (spatial expression at the E7.25 bracket, 2D only).

#### 3.4.6 Peng et al. 2019, Geo-seq spatial transcriptome E5.5-E7.5
- **Citation**: Nature 572:528-532 (2019). DOI 10.1038/s41586-019-1469-8. [H]
- **Content**: laser-capture "zipcode" samples of 20-40 cells; 954 GEO samples covering E5.5, E6.0, E6.5, E7.0 and E7.5 [H]; 3D "corn-plot" maps; an E7.0 single-cell TPM table.
- **Access**: GEO GSE120963 (per-embryo FPKM tables 0.4-3.9 MB; E7.0 single-cell 16 MB); STOmicsDB STDS0000077; portal `egastrulation.sibcb.ac.cn` (responds) [H].
- **Legal use**: T2-embryo: remove E7.5. Everything else is legal for all tasks.
- **Value**: medium for T2-embryo (3D-indexed spatial expression at E6.5-E7.0, close to the E6.75/E7.25 training stages).

#### 3.4.7 Visium whole-embryo sections E9.5-E12.5 (GSE178636)
- **Citation**: Chen Y, et al. Brief Bioinform 24:bbad068 (2023), STEEL. DOI 10.1093/bib/bbad068. [H]
- **Content**: one Visium section each at E9.5, E10.5, E11.5 and E12.5 [H].
- **Legal use**: T1 and T2-heart: E9.5 only. **Value**: low.

### 3.5 Later-stage hearts (E14-P0) to bracket E12.5 (T1, T2-heart)

#### 3.5.1 Feng et al. 2022, chamber-resolved heart atlas E9.5-P9
- **Citation**: Nat Commun 13:7960 (2022). DOI 10.1038/s41467-022-35691-7. [H]
- **Content**: 10x v3 with MULTI-seq. CD1: 18 stages, 72 samples, 29,001 cells after filtering. C57BL/6: 17 stages, 68 samples, 25,605 cells. Hearts split into LA, RA, LV and RV (AV canal and septum included with LV). Mutants: Wt1 (E10.5, E11.5, E13.5, E14.5) and Tbx18 (E14.5, E15.5, E17.5) [H]. Metadata columns include `stage`, `zone`, `strain` and cell-type labels [H].
- **Access**: GEO GSE193346: Seurat `.rds.gz` files (CD1 5.1 GB, C57BL/6 4.7 GB, tbx18 2.9 GB, wt1 1.1 GB) and metadata CSVs (0.8-3.2 MB) [H]. RDS needs R; use R 4.4 or later in an isolated environment (older R had an RDS deserialization vulnerability).
- **Legal use and removal**:
  - T1 and T2-heart: keep WT at `stage` E9.5 and E14.5 or later; remove E10.5-E13.5.
  - Remove mutants for WT modelling. They are legal at E14.5 or later, but they are mutants.
- **Value**: high for T1 (heart-focused late bracket plus an E9.5 anchor on one platform). Medium for T2-heart (expression only).

#### 3.5.2 Feulner et al. 2026, perinatal whole-heart Drop-seq (with Wünnemann et al. 2020 E14.5)
- **Citation**: Feulner L, Wünnemann F, Rafatian N, Liang J, et al. Biol Open (2026; PMID 42725465). DOI 10.1242/bio.062493. E14.5 hearts from Wünnemann F, et al. Nat Genet 52:40 (2020), DOI 10.1038/s41588-019-0536-2. [H]
- **Content**: Drop-seq of whole hearts, 6 biological replicates each at E16.5, E18.5, P0, P4 and P7 (C57BL/6NJ), plus E14.5 (Adamts19 WT n=7 and KO n=8, GSE109247) [H]. The CELLxGENE version has 49,769 cells (E14.5-P7; donors `HeartMat_E14.5_r1` ...) [H].
- **Access**: GEO GSE150817 (`raw_umi_counts.final_cells.tsv.gz` 101 MB; metadata 3.2 MB) and GSE109247; CELLxGENE collection `613f5480-4957-4f80-b804-0e2b85ac454c`, H5AD 0.32 GB at `https://datasets.cellxgene.cziscience.com/525401aa-ac2c-4e45-9a21-c820b6964e0d.h5ad` [H].
- **Licence**: CELLxGENE, CC BY 4.0 [M].
- **Legal use**: T1 and T2-heart: everything (remove Adamts19 KO cells for WT modelling).
- **Value**: high for T1 (whole-heart composition from E14.5; small and easy to handle).

#### 3.5.3 DeLaughter et al. 2016, heart time course E9.5-P21
- **Citation**: Dev Cell 39:480-490 (2016). DOI 10.1016/j.devcel.2016.10.001. [H]
- **Content**: Fluidigm IFC SMARTer; more than 1,200 cells at E9.5 (common atria and primitive ventricle), E11.5, E14.5, E18.5, P0, P3, P7 and P21 (LV, RV and LA free walls from E14.5); Nkx2.5+/- cells at E14.5, P0, P7 and P21 [H].
- **Access**: GNomEx only (`b2b.hci.utah.edu/gnomex`, guest login, accessions 272R, 274R, 275-292R, 439R, 440R). No GEO record for the new data [H].
- **Legal use**: T1 and T2-heart: remove E11.5.
- **Value**: low (small, awkward access).

#### 3.5.4 Goodyer et al. 2019, cardiac conduction system at E16.5
- **Citation**: Circ Res 125:379 (2019). DOI 10.1161/CIRCRESAHA.118.314578. [H]
- **Content**: 10x; micro-dissected E16.5 CD1 hearts in 4 zones (SAN; AVN/His; BB/PF left and right) [H]. Cell count not verified [L].
- **Access**: GEO GSE132658 (`RAW.tar` 166 MB) [H].
- **Legal use**: T1 and T2-heart: legal. **Value**: low-medium (conduction subtypes).

#### 3.5.5 Hu et al. 2026, heart Visium and Slide-seq atlas
- **Citation**: Commun Biol (2026). DOI 10.1038/s42003-026-10259-z. [H]
- **Content**: Visium at E11.5, E14.5, E16.5, E17.5, P0 and P2. Curio Slide-seq at E10.5, E11.5, E12.5, E13.5, E14.5 and E15.5. Plus Wt1-CreERT2 lineage tracing (E14.5) and EPDC ablation (E12.5) [H].
- **Access**: GEO GSE282547 (`RAW.tar` 24.1 GB) [H].
- **Legal use and removal**: T1 and T2-heart: keep Visium E14.5-P2 and Slide-seq E14.5-E15.5; remove E10.5-E13.5 and the E12.5 ablation.
- **Value**: medium (2D spatial heart sections after E13.5, physical scale).

#### 3.5.6 Other late-heart sources
- **Kang et al. 2025 GigaScience**, Stereo-seq heart sections at "E20", P1, P4 and P14; 330,857 bins (bin50); SRA PRJNA1148773; CNGB STT0000062 (PMC11831923) [H]. Legal for T1/T2-heart; low value (perinatal and postnatal).
- **GSE154715**: Wt1-lineage epicardium-derived cells at E12.5 and E16.5, plus MRTF deletion (PMID 34230480) [H]. Keep E16.5 only.
- **ENCODE bulk heart RNA-seq (polyA plus)** [H]:
  - Legal: E14.5 (ENCSR727FHP, ENCSR000CHF), E15.5 (ENCSR597UZW), E16.5 (ENCSR020DGG), P0 (ENCSR526SEX).
  - Excluded: E10.5 (ENCSR049UJU), E11.5 (ENCSR691OPQ), E12.5 (ENCSR150CUE), E13.5 (ENCSR284YKY).
  - Use the legal stages as pseudobulk trajectory anchors for T1's de_score and de_direction. ENCODE release policy: freely available [M].
- **JAX embryos 50-74** and **MOSTA E14.5-E16.5**: see above.
- **3D heart transcriptomics after E13.5 in mouse: none found** (absence not proven).
- **Cross-species option**: MERFISH+ 3D human developing heart at 12 post-conception weeks; 1,835-gene panel; about 3.1M cells; 53 transverse sections 112 µm apart; Kern C, Zhang Q, ..., Chi NC, Qiu X, Zhu Q, Bintu B; bioRxiv 10.1101/2025.11.02.686137; Dryad doi:10.5061/dryad.k98sf7mkx, CC0, 68.9 GB [H]. At 12 pcw this is fetal, well past any mouse E13.5 equivalent. Legality is probable but cross-species staging is not addressed by the rules (A6). Value: low-medium for T2-heart 3D architecture priors.

### 3.6 Morphology and imaging (shape and scale only; same windows apply)
- **McDole et al. 2018 Cell**, in-toto light-sheet imaging and a dynamic atlas, E6.5-E8.5 (DOI 10.1016/j.cell.2018.09.031). IDR `idr0044` (CC BY 4.0) [H].
  - T2-embryo: only frames staged at E7.25 or earlier, or E8.0 or later.
  - T2-heart: drop frames in (E8.25, E8.75).
  - Frames must be staged morphologically. Value: low-medium (growth/shape priors).
- **Esteban et al. 2022 Nat Cardiovasc Res**, 52 segmented embryos, nominal E7.75-E8.5, staged by the myocardial d1/d2 ratio (10 groups, no somite counts) [H]. Mendeley Data `t828xhg66k` (licence not stated) [H].
  - T2-heart: crescent groups only. Do not use the pseudo-dynamic 3D+t models, which pool in-window specimens.
  - T2-embryo: crescent specimens are inside the window.
  - Value: low.
- **DMDD HREM**, E14.5 embryos (607 datasets, 198 WT, per Geyer et al. 2021) at `dmdd.org.uk` (timed out today; licence not found) [M]. Legal for T1/T2-heart. A possible heart-size anchor at E14.5. Value: low.
- **eMouseAtlas EMA 3D models** (Theiler-staged). Use TS15 or earlier, or TS23 or later, for T1/T2-heart; avoid TS11 for T2-embryo [M].

### 3.7 T3: perturbation datasets, legal versus held-out

Held-out conditions: Gata4 KO (val) and Ctnnb1 (beta-catenin) KO (test), both Mesp1-Cre driven, at E8.75. The rule extends to other alleles of the same gene at a comparable stage and to phenocopies.

**Tier A: exclude (same gene, or a direct Wnt/beta-catenin loss in the Mesp1/heart lineage)** [H for each dataset's content]:

| Accession | Content | Why |
|---|---|---|
| GSE165300 | [redacted 2026-10-19: B]. *Mesp1-Cre; Ai9; Wls fl/fl* vs *fl/+*, Mesp1+ cells at E8.0, E8.5, E9.5, 10x v3 with MULTI-seq | Loss of Wnt secretion in the same lineage at the same stages; likely phenocopy |
| [redacted 2026-10-19: B] | | |
| [redacted 2026-10-19: B] | | |
| [redacted 2026-10-19: B] | | |
| [redacted 2026-10-19: B] | | |
| [redacted 2026-10-19: B] | | |
| GSE208162 | single-cell E12.5 heart with inducible endocardial Gata4 KO | same gene |
| [redacted 2026-10-19: B] | | |
| [redacted 2026-10-19: B] | | |
| [redacted 2026-10-19: B] | | |
| [redacted 2026-10-19: B] | | |
| [redacted 2026-10-19: B] | | |

Search strings for re-screening (GEO, mouse, series type): `(Ctnnb1 OR "beta-catenin") AND (heart OR cardiac OR mesoderm OR "heart field") AND (embryo OR embryonic)` returned 32 hits; `Gata4 AND (heart OR cardiac OR endoderm) AND (embryo OR embryonic OR fetal)` returned 121 hits. Re-run before the final submission.

**Tier B: ask before using** (second-heart-field-loss [redacted 2026-10-19: B] mutants that could be called phenocopies): de Soysa Hand2-null (E7.75-E9.25); Jia Isl1-/- (E9.5); Sampath Kumar Tbx6 KO (E9.5); Gata6, Gata5/6 or Zfpm2 (Fog2) mutants; Tbx1, Mef2c, Fgf8/Fgf10 SHF mutants (no single-cell sets checked).

**Tier C: legal (disclose)**: Grosswendt 2020 epigenetic KOs at E8.5 (GSE137337); Pijuan-Sala Tal1-/- chimeras at E8.5; Guibentif T-/- chimeras at E7.5/E8.5 (Wnt-adjacent; disclose); Argelaguet T-KO multiome; Mittnenzweig Foxc1/2 DKO; Weizmann Dnmt3a/b, Tet TKO and Bmp4 series; Li 2016 Nkx2-5-/- (arguably Tier B for Gata4; ask); Feng 2022 Wt1 and Tbx18 mutants; Lescroart Mesp1-/-; Huang 2023 (22 mutants at E13.5).

**Literature-derived gene sets**: DE lists from beta-catenin or Gata4 LOF experiments at E8-E9.5 are "data from another allele at a comparable stage" and should be treated as held-out. Generic pathway gene sets are probably fine; disclose them and ask (open question Q8).

### 3.8 Leakage routes to avoid
- **virtualembryo.ai atlas and predictor**: built on "the same multi-modal mouse embryo resource the benchmark is drawn from" (proposal), with a predictor covering E7.8 to E12.5 [H]. The Data page says using it at target time points "is using the answer, however it is routed". Recommendation: orientation only, never training data.
- **Navigo** (`aristoteleo/Navigo-release`, BSD-2): its tutorials run pretrained checkpoints on JAX-style data with E12.5-E13.5 days (KO notebook subsets E12.5-E13.5) [H]. The checkpoints very likely saw held-out-window data [M]. Do not use them for T1/T2-heart. Retraining the code on filtered data is fine.
- **Mouse or multi-species foundation models pretrained on CELLxGENE or GEO corpora** (e.g. UCE, TranscriptFormer, GeneCompass, Mouse-Geneformer, scPRINT, Nicheformer): if the corpus includes JAX (published on CELLxGENE 2024-03-11), MOCA or MOSTA E10.5-E13.5, the model is "pre-trained on it", so prohibited for T1/T2-heart. Corpus membership for each model was not verified [L]. Human-only models avoid the mouse stage issue but raise A6 if they include human CS13-CS20 tissue.
- **Same-lab data**: GSE176306 E7.5/E7.75 (Section 3.2.2), and any future Chi-lab mouse heart MERFISH or multiome release. None was found today; the MERFISH+ release is human only.

---

## 4. CELLxGENE Discover and Census: what is there and is a stage query practical?

**What exists** (curation API, 2,243 public datasets queried today) [H]. Mouse datasets with embryonic stage terms come from 12 collections. The relevant ones are:
- JAX/Qiu 2024: 31 datasets (one primary raw-counts set of 11,441,407 cells, normalized subsets, 25 major-cluster subsets); TS12-TS27; `donor_id` = `embryo_id`.
- Sampath Kumar 2023 Slide-seqV2: 75 datasets, all coded TS14.
- Feulner 2026 perinatal heart: 49,769 cells; TS22, TS24, TS26, TS27 plus P4 and P7. This equals E14.5, E16.5, E18.5, P0, P4 and P7.

The rest are tissue-specific (limb, gonad, face, thalamus, cerebellum, kidney, lung, islets). Not present: Pijuan-Sala 2019, the extended atlas, MOCA, de Soysa, Tyser, Feng.

**Census facts**:
- Releases: stable LTS `2025-11-08`; latest `2025-11-17` [H].
- Accepted assays include sci-RNA-seq3 (EFO:0030028), Drop-seq and Slide-seqV2 (spatial data sit in `census_spatial_sequencing`) [H].
- Obs keeps only standardized columns: `dataset_id`, `donor_id`, `development_stage(_ontology_term_id)`, `cell_type`, `is_primary_data`, and so on. Author columns such as day, somite count or `celltype_update` are not kept [H].
- JAX inclusion in the current Census build is expected (published 2024-03-11, accepted assay) but was not checked directly [M].
- Client: `cellxgene-census` 1.18.0 (Python 3.10 or later), `tiledbsoma` 2.3.0 [H].

**Practicality verdict**: a stage query is technically easy but semantically unsafe.
- TS13 sits inside the T2-heart window.
- TS15 contains JAX E9.75 embryos (27-29 somites), which are in the extrapolation window.
- TS22 spans E13.5-E15.
- The JAX day-to-TS mapping used by curators is undocumented.

So filter by `donor_id` instead, and add `is_primary_data == True` to avoid duplicates from the cluster subsets. An exact recipe for the DGX (not run here):

```python
import cellxgene_census  # 1.18.0
JAX_RAW = "dcfd4feb-18a3-4b30-81d7-1b0c544a8ab3"  # "Whole dataset: Raw counts only", 11,441,407 cells
T1_OK = [f"embryo_{i}" for i in list(range(1, 25)) + list(range(50, 75))]     # E8.0-E9.5 (0-26 somites) and E13.75-P0
T2H_OK = [f"embryo_{i}" for i in list(range(13, 25)) + list(range(50, 75))]   # E8.75-E9.5 and E13.75-P0
with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    adata = cellxgene_census.get_anndata(
        census, organism="Mus musculus",
        obs_value_filter=f"dataset_id == '{JAX_RAW}' and is_primary_data == True and donor_id in {T1_OK}",
        obs_column_names=["donor_id", "development_stage", "cell_type", "assay"])
# Join donor_id to SCRATCH/dl_meta/jax/sample_sheet.txt (day, somite_count) for disclosure tables.
```

For fine cell types (e.g. "First heart field", "Second heart field", "Epicardium"), the JAX server's `df_cell.csv` (`celltype_update`, `major_trajectory`) is the better source. Pre-cut H5ADs, such as the Cardiomyocytes subset (0.56 GB), are the quickest legal entry point; filter them by `donor_id`.

---

## 5. Legality matrix (strict rules; "partial" means the listed removals are required)

| Dataset | T1 | T2-heart | T2-embryo | T3 |
|---|---|---|---|---|
| Pijuan-Sala 2019 | all | partial (remove E8.5) | partial (remove E7.5, E7.75, mixed) | all; Tal1 chimeras Tier C |
| Extended atlas 2024 | all | partial (remove E8.5, E8.75, 13-somite E9.0, mixed) | partial (remove E7.5, E7.75, mixed) | all |
| JAX Qiu 2024 | partial (embryos 1-24, 50-74) | partial (13-24, 50-74) | partial (2-74) | all |
| MOCA 2019 | partial (E9.5 only) | partial (E9.5 only) | all | all |
| de Soysa 2019 | all WT | partial (E7.75, E9.25; E8.25: ask) | partial (E8.25, E9.25) | WT; Hand2-null Tier B |
| Tyser 2021 | all | partial (crescent stages) | partial (LHT only) | all |
| Lescroart 2018 | all | all | all | all (Mesp1-KO Tier C) |
| Zhang (Chi) 2021, GSE176306 | all | all | partial (E7.25, E8.25) | all |
| Li 2016 | partial (E8.5; E9.5 ask) | partial (E9.5 ask) | all | WT; Nkx2-5-/- ask |
| Jia 2018 | all | partial (remove E8.5) | partial (remove E7.5) | Isl1-/- Tier B |
| DeLaughter 2016 | partial (remove E11.5) | partial (remove E11.5) | all | all |
| Feng 2022 | partial (E9.5, E14.5+) | partial (E9.5, E14.5+) | all | Wt1/Tbx18 Tier C |
| Feulner 2026 / Wünnemann 2020 | all | all | all | all |
| Goodyer 2019 | all | all | all | all |
| Hu 2026 | partial (E14.5+) | partial (E14.5+) | all | all |
| ENCODE heart bulk | partial (E14.5-P0) | partial (E14.5-P0) | n/a | n/a |
| MOSTA 2022 | partial (E9.5, E14.5-E16.5) | partial (same) | all | all |
| Spateo 3D (E9.5, E11.5) | partial (E9.5) | partial (E9.5) | all | all |
| Sampath Kumar 2023 | all | partial (E9.5 only) | all | WT; Tbx6 KO Tier B |
| Lohoff 2022 | all | none | all | all |
| Harland 2025 | all | all | partial (E7.25 embryos) | all |
| Peng 2019 | all | all | partial (remove E7.5) | all |
| Mittnenzweig 2021 | all | all | partial (Et at most 7.25 or at least 8.0) | Foxc1/2 Tier C |
| Grosswendt 2020 | all | partial (remove E8.5) | partial (remove E7.5) | all (Tier C) |
| Huang 2023 | none | none | all | all (Tier C) |
| GSE165300 (Wls, Mesp1) | legal by stage, but avoid | partial at best | all by stage | **excluded** |
| Ctnnb1 / Gata4 datasets (Tier A) | stage rules apply | stage rules apply | stage rules apply | **excluded** |
| virtualembryo.ai atlas / predictor; Navigo checkpoints | avoid | avoid | avoid | avoid |

---

## 6. Ranked shortlist per task, with exact cell-removal rules

Effort figures assume the DGX (H100 box, about 628 GB free disk); all downloads go there, not to the local machine.

### T1, scRNA, predict E12.5 (heart-centred dissection)
1. **JAX Qiu 2024.** Keep `embryo_id` in 1-24 and 50-74 (optionally start at 53 for an E14.0 margin); remove 25-49. Select heart-region trajectories to match T1's dissection: cardiomyocytes, lateral plate and second heart field mesoderm, endothelium/endocardium, epicardium and septum transversum, hepatocytes, foregut epithelium, neural crest, blood. About 40-60 GB via the JAX h5ads plus `df_cell.csv`, or 0.56-44 GB via the CELLxGENE subsets. **Value: highest.**
2. **Feng 2022 (GSE193346).** Keep `stage == E9.5` and `stage >= E14.5`, WT only, both strains; remove E10.5-E13.5 and all mutant objects. About 10 GB of RDS.
3. **Feulner 2026 + Wünnemann 2020 (GSE150817/GSE109247, or the CELLxGENE H5AD).** All stages E14.5-P7; remove Adamts19 KO. 0.3 GB.
4. **Extended atlas 2024.** All cells legal for T1 (E6.5-E9.5). Use E8.5-E9.5 heart lineages with `somite_count` as a fine time axis. 25 GB tarball.
5. **de Soysa 2019 (GSE126128).** WT E7.75, E8.25, E9.25. 3.9 GB CSV.
6. **ENCODE heart bulk.** E14.5, E15.5, E16.5, P0 only (pseudobulk anchors); remove E10.5-E13.5.
7. **MOSTA heart regions** (E9.5; E14.5, E15.5, E16.5) and **Spateo** `mouse_E9.5_heart.h5ad`.
8. **Optional**: Hu 2026 Visium E14.5-P2; Goodyer 2019 E16.5; DeLaughter 2016 (E9.5, E14.5+).
- **Never**: anything in (E9.5, E13.5]: MOCA E10.5-E13.5, MOSTA E10.5-E13.5, Spateo E11.5, ENCODE E10.5-E13.5, Li 2016 E10.5, Feng E10.5-E13.5, Hu Slide-seq E10.5-E13.5, GSE178636 E10.5-E12.5, Huang 2023, JAX embryos 25-49. Also Navigo checkpoints and CELLxGENE-pretrained mouse models.

### T2-heart, 3D MERFISH, predict E12.5 (the E8.5 interpolation board is validation only)
1. **Spateo 3D E9.5** (`mouse_E9.5_heart.h5ad`, `mouse_E9.5_embryo.h5ad`). Keep all; remove every E11.5 file. 0.86 GB.
2. **JAX.** Keep embryos 13-24 and 50-74; remove 1-12 (E8.5 bin) and 25-49 (and validation embryos 85-91). Subset to the 500-gene panel.
3. **Extended atlas.** Keep E6.5-E8.25, E9.0 at 14-16 somites, E9.25 and E9.5; remove E8.5, E8.75, E9.0 at 13 somites, and Mixed gastrulation.
4. **Late anchors for expression and size.** Feng 2022 (E9.5, E14.5+), Feulner 2026, Hu 2026 Visium (E14.5-P2; physical scale), MOSTA E14.5-E16.5 heart regions.
5. **de Soysa** E7.75 and E9.25; remove E8.25 (strict).
6. **Optional morphology**: DMDD E14.5 hearts (size anchor). Esteban 2022 crescent-group meshes only.
- **Never**: Lohoff 2022; Sampath Kumar E8.5 and E9.0; Li 2016 E8.5/E10.5; Tyser LHT (strict); any 6-13-somite embryo; anything in (E9.5, E13.5].

### T2-embryo, 3D MERFISH whole embryo, predict E7.75 (interpolation between E7.25 and E8.0)
1. **Pijuan-Sala 2019 / extended atlas.** Keep E6.5, E6.75, E7.0, E7.25, E8.0, E8.25, E8.5 (and later if wanted); remove E7.5, E7.75 and mixed_gastrulation.
2. **Sampath Kumar 2023 3D Slide-seq.** E8.5 (3-5 somites) and E9.0 virtual embryos; all legal for this task.
3. **Peng 2019 Geo-seq.** E5.5-E7.0; remove E7.5.
4. **Harland 2025 seqFISH.** E7.25 embryos; remove the E7.5 embryos. **Lohoff 2022** E8.5 is legal here.
5. **JAX embryos 2-12** (1-12 somites); remove embryo_1. **Zhang/Chi 2021 GSE176306**: E7.25 and E8.25 only.
6. **Lescroart 2018** (E6.75, E7.25). **Grosswendt 2020** WT at E6.5 and E8.5. **Mittnenzweig 2021** embryos with Et at most 7.25 or at least 8.0.
7. **McDole 2018 light-sheet** frames at E7.25 or earlier, or E8.0 or later (morphology and growth only).
- **Never**: any label E7.5/E7.75, TS11 or head-fold-stage data; de Soysa E7.75; Tyser crescent stages; GSE176306 E7.5/E7.75; Argelaguet E7.5; Nowotschin E7.5; Grosswendt Eed E7.5; T chimera E7.5 pools.

### T3, predict the Mesp1-Cre; Ctnnb1 F/F mutant at E8.75 (Gata4 on validation)
1. **WT Mesp1-lineage references**: GSE176306 (all four stages legal for T3), Lescroart 2018, the extended atlas E8.5-E9.0, JAX E8.5-E9.0. Use them to define lineage programmes and response "receivers".
2. **Legal perturbation priors (Tier C)**: Grosswendt 2020 (10 KOs, E8.5), Tal1-/- and T-/- chimeras, Foxc1/2 DKO, Weizmann Tet/Dnmt3/Bmp4 series, Huang 2023 (22 mutants), Feng Wt1/Tbx18. Use them to learn response size and direction relative to WT, e.g. for severity_slope and de_direction priors.
3. **Tier B only after written clearance**: Hand2-null, Isl1-/-, Tbx6 KO, Nkx2-5-/-.
- **Never (Tier A)**: GSE165300 (Wls in the Mesp1 lineage) and every Ctnnb1 or Gata4 dataset listed in 3.7, including their littermate controls, to keep provenance clean. No literature DE lists from beta-catenin or Gata4 loss of function at E8-E9.5.

---

## 7. Questions to send the organisers before using any of this (ask early in P3)
1. Is training on external public data allowed for T2 and T3? The Data page says T1 only; Rules s.10 and the FAQ say all tasks.
2. Which somite-to-E-day convention applies? What are the somite counts of the challenge's E8.25_late, E8.5, E8.75 and E9.5 heart embryos and the E7.25/E8.0 embryos? Is our strict rule (remove 6-13 somites for T2-heart; remove 27+ somites for T1/T2-heart) acceptable?
3. Does the T2-heart (E8.25, E8.75) exclusion still apply in P3, after E8.5 answers are released and no E8.5 board is scored?
4. Are JAX E13.75 (eMOSS-staged) embryos acceptable as "after E13.5"?
5. Are de Soysa E8.25 (8 somites) and Li 2016 E9.5 (26-29 somites) allowed for T2-heart and T1?
6. Phenocopy scope: confirm that Wls (GSE165300) is held out; rule on Hand2-, Isl1-, Tbx6-, Nkx2-5- and Gata6-perturbation data.
7. Does GSE176306 (Chi lab) share embryos with the withheld E7.75 scRNA or the E7.5/E7.75 MERFISH stages?
8. May literature-derived gene sets (e.g. Wnt-target lists from beta-catenin studies in other tissues or stages) be used for T3?
9. Confirm that pretrained models whose corpora include E9.5-E13.5 mouse data (Navigo checkpoints, CELLxGENE-trained models) are prohibited for T1/T2-heart.
10. Cross-species: are human embryonic heart data at Carnegie stages equivalent to E10-E13.5 treated as held-out? Is the 12-pcw MERFISH+ heart acceptable?

---

## 8. References and access links

Rules and challenge
- Virtual Embryo Challenge rules, data and FAQ pages (local copies `SCRATCH/site_text/challenge_rules.txt`, `challenge_data.txt`, `challenge_faq.txt`); https://virtualembryo.ai/challenge
- eMouseAtlas Theiler stage definitions: https://www.emouseatlas.org/emap/ema/theiler_stages/StageDefinition/stagedefinition.html

Datasets and papers (DOI; main access)
- Pijuan-Sala 2019, https://doi.org/10.1038/s41586-019-0933-9; ArrayExpress E-MTAB-6967/7324/7325: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6967; MouseGastrulationData: https://bioconductor.org/packages/MouseGastrulationData; processed data: https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz; scripts: https://github.com/MarioniLab/EmbryoTimecourse2018
- Imaz-Rosshandler 2024, https://doi.org/10.1242/dev.201867; https://marionilab.github.io/ExtendedMouseAtlas/; https://bioinformatics.stemcells.cam.ac.uk/rlh60/Supplemental/ExtendedMouseAtlas/; https://cells-test.gi.ucsc.edu/?ds=ext-mouse-atlas; E-MTAB-11763, E-MTAB-13645
- Qiu C 2024 (JAX), https://doi.org/10.1038/s41586-024-07069-w; https://omg.gs.washington.edu/jax/public/download.html; https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/download/; GEO GSE228590, GSE186069; CELLxGENE https://cellxgene.cziscience.com/collections/45d5d2c3-bc28-4814-aed6-0bb6f0e11c82; code https://github.com/ChengxiangQiu/JAX_code
- Cao 2019 (MOCA), https://doi.org/10.1038/s41586-019-0969-x; GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119945
- Chen A 2022 (MOSTA), https://doi.org/10.1016/j.cell.2022.04.003; https://db.cngb.org/stomics/mosta/download/
- Qiu X 2024 (Spateo, "molecular holograms"), https://doi.org/10.1016/j.cell.2024.10.011; files under https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/; MOSTA3D preprint https://doi.org/10.1101/2024.08.17.608366; code https://github.com/aristoteleo/spateo-release
- de Soysa 2019, https://doi.org/10.1038/s41586-019-1414-x; GEO GSE126128; https://mouse-cardiac.cells.ucsc.edu/
- Tyser 2021, https://doi.org/10.1126/science.abb2986; E-MTAB-7403; https://content.cruk.cam.ac.uk/jmlab/mouseEmbryonicHeartAtlas/
- Lescroart 2018, https://doi.org/10.1126/science.aao4174; GEO GSE100471
- Zhang Q (Chi lab) 2021, https://doi.org/10.1161/CIRCRESAHA.121.318943; GEO GSE176306
- Li G 2016, https://doi.org/10.1016/j.devcel.2016.10.014; GEO GSE76118
- DeLaughter 2016, https://doi.org/10.1016/j.devcel.2016.10.001; GNomEx https://b2b.hci.utah.edu/gnomex
- Jia 2018, https://doi.org/10.1038/s41467-018-07307-6; ENA PRJEB23303
- Lohoff 2022, https://doi.org/10.1038/s41587-021-01006-2
- Harland 2025, https://doi.org/10.1016/j.celrep.2025.116047; Zenodo https://zenodo.org/records/13977985
- Peng 2019, https://doi.org/10.1038/s41586-019-1469-8; GEO GSE120963; http://egastrulation.sibcb.ac.cn/
- Sampath Kumar 2023, https://doi.org/10.1038/s41588-023-01435-6; GEO GSE197353; CELLxGENE https://cellxgene.cziscience.com/collections/d74b6979-efba-47cd-990a-9d80ccf29055
- Mittnenzweig 2021, https://doi.org/10.1016/j.cell.2021.04.004; GEO GSE169210
- Grosswendt 2020, https://doi.org/10.1038/s41586-020-2552-x; GEO GSE137337
- Ibarra-Soria 2018, https://doi.org/10.1038/s41556-017-0013-z; E-MTAB-6153
- Nowotschin 2019, https://doi.org/10.1038/s41586-019-1127-1; Argelaguet 2019, https://doi.org/10.1038/s41586-019-1825-8
- Huang 2023, https://doi.org/10.1038/s41586-023-06548-w
- Feng W 2022, https://doi.org/10.1038/s41467-022-35691-7; GEO GSE193346
- Feulner 2026, https://doi.org/10.1242/bio.062493; GEO GSE150817; Wünnemann 2020, https://doi.org/10.1038/s41588-019-0536-2, GEO GSE109247; CELLxGENE https://cellxgene.cziscience.com/collections/613f5480-4957-4f80-b804-0e2b85ac454c
- Goodyer 2019, https://doi.org/10.1161/CIRCRESAHA.118.314578; GEO GSE132658
- Hu J 2026, https://doi.org/10.1038/s42003-026-10259-z; GEO GSE282547
- Kang 2025 (GigaScience), https://pmc.ncbi.nlm.nih.gov/articles/PMC11831923/
- Chen Y 2023 (STEEL; Visium E9.5-E12.5), https://doi.org/10.1093/bib/bbad068; GEO GSE178636
- ENCODE heart RNA-seq, e.g. https://www.encodeproject.org/experiments/ENCSR727FHP/ (E14.5)
- MERFISH+ human heart, https://doi.org/10.1101/2025.11.02.686137; Dryad https://doi.org/10.5061/dryad.k98sf7mkx
- McDole 2018, https://doi.org/10.1016/j.cell.2018.09.031; IDR https://idr.openmicroscopy.org/study/idr0044/
- Esteban 2022, https://pmc.ncbi.nlm.nih.gov/articles/PMC11357989/; Mendeley https://data.mendeley.com/datasets/t828xhg66k/1
- [redacted 2026-10-19: B]
- [redacted 2026-10-19: B]
- [redacted 2026-10-19: B]

Infrastructure and terms
- CELLxGENE curation API: https://api.cellxgene.cziscience.com/curation/v1/datasets; Census schema: https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md; Census releases: https://census.cellxgene.cziscience.com/cellxgene-census/v1/release.json; Census licence (CC BY 4.0): https://registry.opendata.aws/biohub-cellxgene-census/
- EMBL-EBI terms of use: https://www.ebi.ac.uk/about/terms-of-use/
- GEO disclaimer ("no restrictions on the use or distribution"): https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html
- Navigo: https://github.com/aristoteleo/Navigo-release

Local evidence files (metadata only): `SCRATCH/dl_meta/jax/sample_sheet.txt` (JAX embryo staging), `SCRATCH/dl_meta/ema/*.bin.gz` and `dataset.json` (extended-atlas stage and somite fields), `SCRATCH/dl_meta/cxg/all_datasets.json` (CELLxGENE index), `SCRATCH/dl_meta/geo/*.txt` (GEO records), `SCRATCH/scripts_extdata/jax_rules.py` (per-task JAX embryo lists). SCRATCH = `SCRATCH`.
