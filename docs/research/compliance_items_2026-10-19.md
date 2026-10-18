# Compliance items and flagged uncertainties, 2026-10-19

Companion to `digest.md`. Every protected-window (W) item the workflow produced, every size/leak (A) exposure, and
every judgement call I was unsure about. **No banned content appears here**: where an item exists only because of
material you must not see, it is described without its content.

Source labels: `dim:<name>` = research dimension, `ver:<name>` = its adversarial verifier, `gap:<n>` = gap-fill,
`critic` = completeness critic.

---

## 1. Window items (W): external data in or near a protected window

The windows (from `docs/DATA_POLICY.md`): T1 and T2-heart exclude **after E9.5 up to and including E13.5**
(exactly E9.5 allowed; after E13.5 allowed with disclosure). T2-heart interpolation excludes **strictly between E8.25
and E8.75**. T2-embryo excludes **strictly between E7.25 and E8.0**. A shared model respects the union.

### W-1. JAX / Qiu 2024 atlas — which embryos were used, and is that enough to confirm compliance?

**What the reports say.** The external-direction gap-fill states at the top of its report and in its compliance log
that it used **JAX embryos 1-22 and 54-74 only**, plus T1 E8.5/E9.5 and external heart data at E9.5 and from E14.5,
and that **no E10.5 or E12.5 data was used**. Specifically:
- `df_cell.csv` was **stream-filtered by embryo_id**, so in-window rows (embryos 25-49) were never written to disk.
- A per-embryo count-file route was used for embryos 1-22 and 54-74, so no in-window expression was downloaded.
- **Embryos 23-24 were excluded**: they are staged E9.5 by somite count but come from a litter harvested as "E9.75",
  and the report applies the strict any-label rule.
- **Embryos 50-53 were not used.** Their morphometric stage is E13.75-E14.0, and the external-data verifier showed
  separately that embryos 50-51 come from a **gestational E12.5** litter and 52-53 from a **gestational E13.5** litter.
  The corrected recipes are `T1_OK = 1-24 + 54-74` and `T2-heart_OK = 13-24 + 54-74`.
- A one-off read of an atlas sub-analysis metadata file covering all stages was inner-joined to legal cells only, and
  the raw file was then deleted.
- The T1 composition gap-fill read the same metadata for **embryos 1-24 and 50-74** (labels only, no expression), and
  its own rule-risk note says embryos 50-53 were excluded for margin; its reported numbers use embryos 18-24 and
  55-72. **This is an internal inconsistency between its action list ("register embryos 1-24 and 50-74") and its
  rule-risk note ("50-53 excluded"), and it must be resolved before anything is registered.**

**Do the stages respect the windows?** For the clearly in-window block, yes: embryos 25-49 were never read, and the
late block starts at 54. Three things are **not** confirmed by the reports:
1. **The boundary at E9.5.** Embryos up to 22 are at or before E9.5 by somite staging, but the organisers' somite-to-
   E-day convention is unpublished, and the external-data verifier showed the same somite count maps to different
   labelled days in different datasets. Whether embryos 18-22 (and 23-24, currently excluded) are "exactly E9.5 or
   earlier" on the organisers' scale is an open question that only they can settle.
2. **The T2-heart interpolation window.** The gap-fill's own rule-risk note states that **JAX embryos 1-12 (the E8.5
   bin) fall inside (E8.25, E8.75)**. Its results are therefore valid for T1 only and **must not feed any T2-heart
   decision or any model shared with T2-heart**. Nothing in the reports verifies that this separation is enforced in
   code.
3. **Gestational versus morphometric age.** The verifier's finding that the atlas's "day" is a morphometric estimate
   that can differ from the litter's gestational age applies in principle to the kept embryos too; it was checked
   only for embryos 50-53.

**What the lead must do.** (a) Fix the inconsistency above and register exactly one embryo list per task in
`EXTERNAL_DATA.md`, with the committed filter code, before the data enters any model. (b) Enforce the T1-only
restriction on embryos 1-12 in code, not by convention. (c) Put the somite-boundary question and the "does E13.75
count as after E13.5" question to the organisers in writing before the first submission. (d) Never validate or
evaluate on in-window embryos either — evaluation on held-out-window data is also use.

### W-2. External heart dataset (10x, chamber-dissected) — in-window rows present in the raw files
`gap:3`, `gap:1`, `ver:external_data`. The raw object downloaded to the DGX contains E10.5-E13.5 cells; only a
filtering script that drops them before writing was used, and the full metadata files contain in-window rows that
were dropped at read time. **Check:** quarantine or delete the raw object, keep the filter script in the repo, and
make sure an older unfiltered extract in the scratch tree is never reprocessed. Also note the per-cell stage column
must be used because the GEO samples are multiplexed pools.

### W-3. Disclosed slip: in-window QC metadata printed once
`gap:3` rule-risk log. One exploratory print showed per-stage medians of library size and mitochondrial percentage
for **all** stages of an external dataset, including E10.5-E13.5. These are QC metadata, not expression, and the
report states they were not used anywhere. **Check:** record the exposure in the provenance log and confirm no
parameter was set from it.

### W-4. External sources used but not yet registered
`gap:1`, `gap:3`. Label-only metadata from three external sources (two heart datasets at E9.5 and E14.5+, and the
atlas labels) are on the DGX and were used to build the T1 composition prior and the direction comparison, but no
`EXTERNAL_DATA.md` rows exist yet. The project policy requires a row and committed filter code **before** a source
enters a model. **Check:** register them, move the filter scripts into the repo, and disclose the author-stated ages
(one source's cells are coded at a Theiler stage that spans the window boundary while the author age is legal).

### W-5. A public Visium/Slide-seq heart dataset mixes stages on one capture area
`ver:external_data`. One file labelled as a late stage physically contains an in-window heart section on the same
capture area, and another mixes a postnatal stage with E14.5. **Check:** remove those spots by tissue position or skip
the sample, and describe the removal. Its E14.5 Slide-seq samples are tamoxifen-treated reporter animals, not plain
wild type.

### W-6. Atlases that span a window
`dim:external_data`, `dim:t1_methods`, `dim:landscape`. The large mouse atlases and spatial resources used or
considered (the time-lapse nuclei atlas, the Stereo-seq resource, the extended gastrulation atlas, the 2019
gastrulation atlas, the 3D Stereo-seq embryo/heart files, ENCODE heart bulk) each contain in-window material.
**Check:** in-window cells must be removed **before any processing** (normalisation, HVG selection, integration,
training), by code committed to the repo, and the removal described in the method summary. Exactly E13.5 is excluded.

### W-7. In-window developmental literature used for E12.5 priors — **removed from the digest**
`dim:t2_spatial` section 4.1 and key finding 7; `ver:t2_spatial` verdict 7 and correction 4; `gap:1` rule-risk 0;
`dim:t1_methods` (emergent-state class list). The research used published developmental timing from **after E9.5
through E13.5** to set the E12.5 composition prior, the emergent-population list and the anatomical niches for
placement. The verifier flagged this as inconsistent with the same report's own rule note, and recommended taking
quantitative priors only from sources at or before E9.5 or after E13.5 until the organisers answer. A question on
whether qualitative in-window background may inform priors was raised in the organisers' channel on 2026-10-16 and is
**still unanswered**. **Check:** decide whether the E12.5 class list may rest on this at all; if yes, disclose it; if
no, derive the class list from released stages plus legal post-E13.5 data only. The removed material is the only
reason this item exists; its content is not reproduced anywhere in the screened outputs.

### W-8. Contaminated pretrained models and checkpoints
`dim:t1_methods`, `dim:external_data`, `ver:external_data`, `ver:t1_methods`, `dim:t2_spatial`, `dim:architecture`.
Named model families whose corpora contain mouse embryo cells inside a window, or which were trained on an atlas
spanning it, are listed in the digest (section 2.1) as **do not use**. One flow-matching checkpoint family is
confirmed contaminated (a knockout checkpoint run on in-window data; an interpolation checkpoint trained across
21 time points spanning the window). A benchmark model whose evaluation design the team wants to copy is also
contaminated — copy the design, never the weights. **Check:** audit corpus **and training date** for every pretrained
model before use (anything trained on the public cell atlas after 2024-03-11 may contain the in-window embryos);
extend the audit to T2-embryo (gastrulation atlases carry E7.5/E7.75) and to T3 (perturbation corpora); and remember
a shared model must respect the union of windows.

### W-9. The organisers' own atlas and predictor
`dim:external_data` rule-risk 6. Explicitly prohibited at target time points. **Check:** keep it out of every
pipeline and out of any "related work" code path.

### W-10. Cross-task use of released data that falls inside another task's window
`dim:architecture` key finding 6 and rule-risks 2-3; `ver:architecture`. Before 2026-10-20, T1's released E8.5 scRNA
is measured data at a stage held out for T2-heart and must not feed a T2-heart pipeline. In the test phase, cross-modal
use of E10.5 data (T1 scRNA for T2-heart E12.5, heart MERFISH for T1 E12.5) is ambiguous because E10.5 sits inside the
stated exclusion window. **Check:** keep two variants until the organisers answer; this is organiser question Q-A/Q-B
in the architecture report.

### W-11. A public dataset from the challenge's own data-generating lab
`dim:external_data` key finding 3 and rule-risk 4; `ver:external_data`. Its E7.5 and E7.75 samples sit inside the
T2-embryo window (E7.75 is the hidden test stage). Literal cell overlap with the withheld challenge data is unlikely
(sorted lineage cells vs whole embryos) but batch overlap cannot be excluded. **Check:** never use those samples for
T2-embryo or in a shared model that serves it; ask the organisers whether the dataset shares embryos or protocol
batches with the withheld material. The other stages are a legal wild-type lineage reference for T3.

### W-12. "Head-fold stage" is effectively the hidden T2-embryo test stage
`ver:external_data` correction 4. The data-generating lab stages embryos morphologically, and its convention places
the head-fold stage at E7.75. **Check:** drop head-fold-staged data from **any** lab from T2-embryo work, not just
data labelled E7.75.

### W-13. Somite and Theiler staging make label-only filtering unsafe
`dim:external_data` key finding 1; `ver:external_data` verdict 1 and correction 4; `dim:t1_methods` rule-risk 2.
The same somite count maps to different labelled days across datasets, and Theiler stages straddle window boundaries.
The strict predicates the report proposes are **not convention-proof**. **Check:** ask the organisers for the somite
counts of their own stages and their mapping; until then filter on original somite/bin metadata, not atlas-curated
Theiler labels, and keep the mapping in the repo.

### W-14. Cross-species data at equivalent stages
`dim:external_data` open question 8; `dim:t1_methods` open question 2. The rules are silent on human embryonic data at
Carnegie stages equivalent to a protected window, including human data inside human-only foundation models.
**Check:** ask before using any such source.

### W-15. The Rules / Data-page conflict on external data for T2 and T3
`dim:external_data` key finding 5; `ver:external_data` verdict 4; `ver:t3_biology` verdict 4 (**refuted** the stronger
"three places" version: only one sentence carries the exclusivity wording and it is about training on external
single-cell data, not about reading papers). The binding text is Rules s.10, which permits external data for all tasks
with disclosure and windows. **Check:** get written confirmation before any T2/T3 submission uses external data, and
keep a fully internal fallback pipeline that uses released stages only.

### W-16. Does the genotype exclusion for external data survive the validation release?
`ver:external_data` correction 9; `ver:t3_methods` correction 7; `critic` contradiction 2. Rules s.10 still lists the
Gata4 knockout at E8.75 as held out, while the Data and timeline pages say validation answers become training
material at P3, and s.10 also treats "a perturbation that phenocopies them" as the held-out condition. The T3 plan
depends on using the released Gata4 cells inside the beta-catenin prediction. The drafted question was **omitted from
the email that was actually sent**. **Check:** send it and get a written answer before the 16 November freeze; mark
that pipeline component "pending confirmation" until then.

### W-17. Disclosure is itself a rule
`dim:external_data` rule-risk 9; `dim:t3_biology`; `dim:landscape` rule-risk 5. An undisclosed external source is a
violation even if it changed no result, and the disclosure must cover sources that shaped the method without changing
a number. **Check:** maintain `EXTERNAL_DATA.md` as you go and have the disclosure text ready before the first test
submission.

---

## 2. Size and leak items (A)

### A-1. The scorer's pre-release artefacts contain held-out geometry
`dim:t2_spatial` section 1.5; `dim:architecture` key findings 4-5; `dim:metrics` correction 0; `dim:engineering`
correction 0; `dim:landscape` correction 1; `gap:4` rule-risk 1; `gap:5`. Already recorded in `docs/DATA_POLICY.md`.
**Check:** (a) never `pip install` the PyPI release of the scorer anywhere, including the DGX — it is the unmerged
pre-release code; pin the git commit; (b) nobody building the size model may open the pull request or the wheel's
shape-metrics file; (c) the private report to the organisers should be widened to cover the wheel and the pull-request
description as well, and should also ask which build and metric set the test boards and the final evaluation will use
and whether anchors will be re-baked (the sent email covers only the metric change).

### A-2. Copies of the leaked file inside this session's scratch tree
`ver:architecture` correction 1 lists four local paths (three copies of the shape-metrics file from the wheel and the
pull request, plus a saved copy of the pull-request body) under
`/tmp/claude-1001/-mnt-md0-IITM-BackUp-Home-vijayavallabh/07d440b9-3ae8-4421-8a08-a046b42f307d/scratchpad/` in
`dl/veckit012/`, `factcheck/veckit_cmp/pr6/`, `factcheck/veckit_cmp/whl012/x/` and `landscape_work/gh/`.
**Check:** delete or quarantine them. I did not touch them (the brief forbids editing inputs).

### A-3. A held-out-derived spacing bound entered a research note
`ver:t2_spatial` verdict 3 and correction 0. One number quoted in `research/t2_spatial.md` section 2.2 (and in its
summary) does not match any released stage and matches a held-out column of the leaked table. It is removed from the
screened copy. **Check:** confirm that no size, spacing or density parameter anywhere in the plan or code was set from
it, and record the exposure in the data-provenance log.

### A-4. Board anchors reveal target size relations
`dim:t2_spatial` rule-risk 1; `dim:metrics` rule-risk 0; `dim:architecture` rule-risk 1; `dim:engineering`
correction 0; `gap:4` rule-risk 2. The published validation anchors let the target size be recovered directly, and
several reports used anchor-derived previews. The Evaluation page says nothing about the test boards is published.
**Check:** predict scale only from released files; never derive it from an anchor or from score inversion; state in
the method summary that anchors were not used for this; and re-derive every preview from the released files once they
are out.

### A-5. Conclusions whose rationale came partly from the leaked table
`ver:architecture` correction 2. Two qualitative claims in the T2-heart plan ("cell spacing is constant across
stages", "the imaged field of view varies by stage") were partly derived from it, and one report's per-stage geometry
table was quoted from it. The qualitative field-of-view conclusion is also independently supported by released-stage
measurements in `gap:4`. **Check:** re-derive any spacing, size or field-of-view prior from released stages only,
commit that derivation before the first test submission, and keep the size-model author away from the affected report
sections (they are redacted in the screened copies).

### A-6. Score-based recovery of target properties
`dim:metrics` rule-risk 0; `dim:landscape` rule-risk 3; `dim:engineering` rule-risk 1; `ver:t2_spatial`
correction 10. Iterative try-and-score is explicitly allowed by Rules s.10; what is prohibited is solving for a target
property — inverting the skill map, differencing submissions, or deriving a size from a published floor. One public
community repository did exactly this. **Check:** keep "no bisection" as a team policy but describe it correctly in
the method summary, and make sure the second submission is not built from a property read off the first submission's
per-metric scores (`ver:t3_methods` correction 8; `critic` contradiction 8).

---

## 3. Judgement calls I made, and what you may want to revisit

Where I was unsure, I removed the material and recorded it here, as the brief requires.

- **U-1. Canonical pathway readout genes (B).** Several reports name the classic target genes of the held-out gene's
  pathway, as an off-panel activity signature to be learned from wild-type data. These are generic pathway knowledge,
  not knockout results, and the 2026-10-15 ruling explicitly allows generic pathway and gene-set databases — but they
  are also "named target genes tied to" the held-out gene, so I removed every such list. **You may rebuild the list
  from a generic pathway database**, provided it was not built from perturbing the held-out genes, and record its
  provenance. The method around it (build the activity score in the 32k-gene scRNA space at cluster or neighbourhood
  level, regress panel genes on it within cell type, then project) is kept in full.
- **U-2. Excluded-source lists (B).** Two reports keep explicit lists of publications and dataset accessions to
  avoid for T3. A list of citations that state the banned results is itself banned material under the brief, so every
  publication list and the extended accession list are removed from the digest and from the screened copies. The two
  accessions that `docs/DATA_POLICY.md` already names (the approved-for-T3 endocardial dataset and the lineage
  Wnt-secretion phenocopy) are left where the screened text simply restates the policy, since you read that file at
  startup anyway. The blanket rule covers the rest without naming it. If you want a machine-checkable blacklist for
  the provenance filter, have it generated into a file that the T3 modeller does not read.
- **U-3. A reported relationship between the two knockouts' responses (B).** One community source and one verifier
  correction state how the released training knockout's response relates to the Gata4 knockout's. It is not a paper
  result, but it is a prior about a knockout of a held-out gene, so I removed it. **It can be measured directly from
  released data after 2026-10-20**, which is the compliant route and is already in the pre-registered gap-fill 2
  script.
- **U-4. Broad phenotype statements I kept.** The digest keeps compartment-level phenotype statements for both
  evaluation genotypes (second-heart-field/outflow-tract expansion failure, right-ventricular hypoplasia, thin
  myocardium, reduced progenitor proliferation, first-heart-field relatively spared). The 2026-10-15 ruling allows
  broad background of exactly this kind, and the brief says broad phenotype statements stay. If you read the ruling
  more strictly, treat these paragraphs as removable; nothing downstream depends on them numerically.
- **U-5. Scope factors measured on a released stage (A).** Gap-fill 4 measured how much the radius shrinks when the
  cloud is restricted to cardiac compartments, **on a released stage**. I kept those factors (released data) but
  removed every value that applies them to a held-out stage. If you would rather treat the factor itself as a
  held-out-stage scale factor, drop it and re-measure on E10.5 after the release — which the report recommends anyway.
- **U-6. Shape-metric anchors kept.** I kept the published `d2_shape` and `occupancy_dice` floor/ceiling anchors and
  skill arithmetic: they are shape-overlap statistics, not sizes, and they are needed to reason about the scorer.
  I removed the `scale_log_ratio` anchors and everything derived from them. Flagging this in case you want the
  stricter reading.
- **U-7. Released-stage geometry kept.** Per-stage cell counts, section steps, extents, principal SDs and
  nearest-neighbour spacing measured from the released files are kept in the screened `gap_t2_units_scale_shape.md`.
  Values for the same released stages that were quoted from the leaked pre-release table were removed, because the
  brief bars radius values from that source for **any** stage.
- **U-8. The parity hedge (gap-fill 5).** Not a data-policy item, but an ethics one the report itself raises: the
  hedge exploits a known scorer defect, uses no score returns and only orients our own prediction. The report
  recommends disclosing it or asking the organisers in the same message as the build question. **Your call.**
- **U-9. Competitor-derived numbers.** T3 and T2 method lessons taken from public community repositories are kept
  (they are scoring results, not target properties), with the verifier's caution that those repositories are
  unlicensed: ideas only, no code, and disclose any reuse. One clause was removed under U-3.

---

## 4. State of the inputs while I was screening

- **Another agent was redacting `scratchpad/research/*.md` in place during this screen.** Between my first read and
  my last write, `community_digest.md`, `metrics.md`, `t2_spatial.md` and `gap_t2_units_scale_shape.md` changed on
  disk (their redaction markers are dated 2026-10-19 and sit on the same lines mine do). My screened copies were
  produced from the state of 2026-10-19 00:44-00:50 IST and each redaction was verified to land on the intended text.
  **Treat the screened copies, not the originals, as current**, and re-run a screen if you need a later version.
- **Seven reports were not screened**: `t1_methods.md`, `external_data.md`, `metrics.md`, `landscape.md`,
  `architecture.md`, `engineering.md`, `community_digest.md`. Their substance is in `digest.md` section 2. The
  originals contain A and B material (held-out-stage size values; gene-level knockout findings and the citations that
  state them) and must not be opened by anyone building the T3 or heart-size methods.
- The workflow's raw structured output and the per-agent journal also contain that material. They are not safe to
  read directly.
