# Task 3 biology: what we may use (redacted 2026-10-18)

**Redaction note.** The original version of this report (2026-10-18) contained gene-level findings from published
beta-catenin (Ctnnb1) and Gata4 loss-of-function studies in cardiac mesoderm: which genes go up or down, and by how
much. On 2026-10-15 the organisers ruled that such gene-level results count as data from the held-out condition even
when read off a paper (`docs/DATA_POLICY.md`). Those findings, the gene-by-gene direction tables derived from them and
the paper titles that state them were removed from this file and from every other note in the repo. What remains is
allowed: facts about the target and panel, broad phenotype background, Mesp1-lineage scope, the biology of the
released training knockout (Mab21l2), and the rules analysis. **Task 3 gene directions must come only from released
data, by the pre-registered procedure in `docs/PLAN.md` section 4.4.**

## 1. The target

| Item | Value |
|---|---|
| Test condition | Mesp1-Cre; beta-catenin F/F at E8.75, 2 biological replicates (300 + 381 MB) |
| Validation condition | Mesp1-Cre; Gata4 F/F; Gata6 F/+ at E8.75, 2 replicates; answers released at the start of the final phase (2026-10-20) |
| Training condition | Mab21l2 KO at E9.5 (458 MB); allele not stated on the site |
| Reference | matched WT at E8.75 (shared with the T2 heart setting) and at E9.5 |
| Panel | the 500-gene MERFISH panel, identical for T2 heart and T3 |
| Not on the panel | Ctnnb1 and [redacted 2026-10-19: B] |
| Scored | de_score 0.30, de_direction 0.25, abs(severity_slope) 0.25, mutant-state distribution 0.20 (mmd_u 0.6 / variogram 0.4); all DE terms use the pseudobulk over every cell in the file, differenced against the matched WT |

Because the scored quantities are per-gene signs and ranks of the pseudobulk change, a hand-written gene-direction
table is close to the scored answer itself. That is why the organisers' ruling matters so much here, and why the
procedure must derive every gene direction from released data.

## 2. Broad phenotype background (allowed under the 2026-10-15 ruling)

- **Mesodermal beta-catenin loss, around E8.5-E9.5:** a second-heart-field expansion failure. The outflow tract and
  right ventricle are hypoplastic and cardiac looping is disrupted, while the first-heart-field / left-ventricle
  component is comparatively preserved. Progenitor proliferation is reduced; this reads as a growth and allocation
  defect rather than a failure of specification or differentiation. At E8.75 (about 8-12 somites) the mutant heart is
  expected to look visibly abnormal.
- **Beta-catenin is also an adherens-junction component.** Loss in endothelium alters vascular patterning while early
  vasculogenesis proceeds; [redacted 2026-10-19: B].
- **Mesp1-Cre; Gata4 F/F; Gata6 F/+ (validation, unpublished genotype):** by analogy with related Gata4/Gata6 losses,
  expect thin myocardium, right-ventricular hypoplasia and reduced cardiomyocyte proliferation; complete loss of both
  genes blocks cardiomyocyte differentiation. Both knockouts are Mesp1-Cre at E8.75 [redacted 2026-10-19: B],
  which is why the released Gata4 data is the best legal calibration for the test board.

## 3. Mesp1-lineage scope (lineage biology, not knockout data)

- Mesp1-Cre labels early cardiac and cranial mesoderm and their derivatives: heart-tube myocardium, endocardium,
  pericardium/proepicardium and pharyngeal mesoderm, plus extra-embryonic mesoderm. It labels trunk paraxial mesoderm
  poorly. Endoderm, surface ectoderm, neural tube and neural crest are not Mesp1-derived (Saga et al. 1999,
  doi:10.1242/dev.126.15.3437; Kitajima et al. 2000, doi:10.1242/dev.127.15.3215).
- Practical consequence: non-Mesp1 compartments in a knockout sample should not change cell-autonomously, so their
  changes estimate nuisance (depth, section coverage, staging, batch) plus non-autonomous effects. `configs/lineage_map.yaml`
  carries a `mesp1` flag per coarse group for this purpose.
- The common Mesp1-Cre allele is a knock-in (Mesp1^tm2(cre)Ysa), so Cre carriers are Mesp1 heterozygous; whether the
  matched WT carries Cre is not stated.

## 4. The released training knockout: Mab21l2 (not a held-out genotype)

- Yamada et al. 2004 (doi:10.1016/j.ydbio.2004.07.016): Mab21l2-deficient mice have a rudimentary retina, defective
  ventral body wall formation and die in utero at mid-gestation. Probably a germline null; the released file's
  `genotype` column is `mab21l2` for every cell and the probed Mab21l2 transcript should be checked in the data.
- Saito et al. 2012 (doi:10.1371/journal.pone.0032991): Mab21l2 is expressed in trabecular and compact myocardium and in
  the septum transversum mesenchyme including the proepicardium. No overt heart defect at E8.5-E9.5; at E9.5 the
  septum transversum / proepicardium is reduced (smaller Tbx18+/Wt1+ proepicardium, reduced Itga4), with myocardial
  thinning appearing from E10.5-E11.5. Nkx2-5, Gata4 and Bmp10 were reported unchanged.
- So the Mab21l2 E9.5 signal is expected to be mostly a proepicardium/septum-transversum deficit plus nuisance (the
  released KO is about 27% deeper per cell than its WT; `docs/DATA_NOTES.md`). It is a different stage, possibly a
  different allele type, and a different tissue from the test knockout, which is why transferring its pseudobulk delta
  fails.

## 5. Pseudobulk algebra (method, no knockout data)

The DE reference is the matched WT; the pseudobulk is the mean over all cells:

    dpb_g = sum_t [ f_t^KO * m_gt^KO - f_t^WT * m_gt^WT ]

This splits into a composition term (cell-type fractions change) and a within-type term (expression changes inside a
type, from biology, depth or batch). On sparse MERFISH (median ~5 log units per detected count), only highly
expressed, compartment-specific genes cross the |lfc| >= 0.25 DE threshold from a modest composition change; many
within-type changes are invisible to de_score but visible to de_direction. A composition change and its consequences
for each gene must therefore be computed from the released WT expression profiles of the affected cell states, not
written down by hand.

## 6. Rules outcome

- Official Rules s.10 treats "data from the held-out knockouts, from another allele of the same gene at a comparable
  stage, or from a perturbation that phenocopies them" as held out.
- Organiser ruling, 2026-10-15 (Slack, Weize Xu): gene-level results of losing the held-out gene in embryos at a
  comparable stage are not allowed even when read off a paper; general background on what the phenotype looks like in
  broad terms is allowed; generic pathway and gene-set databases are allowed except gene sets built from a knockout or
  perturbation of the held-out gene.
- Consequences for us: no gene-level beta-catenin or Gata4 knockout findings, no curated target lists derived from such
  experiments, no LLM-elicited gene priors (models trained on that literature), no phenocopy datasets (GSE165300,
  Mesp1-Cre; Wls, is excluded). Broad phenotype (section 2), lineage scope (section 3), Mab21l2 biology (section 4) and
  released data are the allowed inputs.

## 7. Excluded sources (listed so they are avoided, not used)

Beta-catenin loss-of-function studies in cardiac mesoderm ([redacted 2026-10-19: B]
[redacted 2026-10-19: B]
[redacted 2026-10-19: B]) and Gata4/Gata6 loss-of-function studies ([redacted 2026-10-19: B]
[redacted 2026-10-19: B]). Their
findings are not reproduced anywhere in this repo.
