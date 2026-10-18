# Data policy

Condensed from the Official Rules s.10 and the Data page (https://virtualembryo.ai/challenge/rules,
https://virtualembryo.ai/challenge/data), as read on 2026-10-18. When in doubt, the live rules win; re-read them before
each submission. Any change to this file needs a dated note at the bottom.

## Competition data

- Released for research and education within the Challenge. Never commit it, share it with other teams or post it.
- Unpublished (Qingquan Zhang and Neil Chi, UCSD). Any publication or public presentation needs prior written approval
  from Dr. Neil Chi's group and must credit Dr. Qingquan Zhang and Dr. Neil Chi, plus acknowledge the Virtual Embryo
  Challenge (Terms of Use clause 14).
- Validation answers are released at the start of the final phase (2026-10-20) and may then be trained on. Test ground
  truth is never released.

## Held-out stages and genotypes (no measured data from these, by any route)

| Task | Held out |
|---|---|
| T1 (scRNA) | E10.5 (val), E12.5 (test) |
| T2 embryo | E7.5 (val), E7.75 (test) |
| T2 heart | E8.5 (val interp), E10.5 (val extrap), E12.5 (test) |
| T3 | Gata4 KO and beta-catenin KO at E8.75 |

"By any route" includes pretrained models whose training data contains the window, and public atlases that contain it.

## External data windows

- **Extrapolation targets (T1, T2 heart):** external data excluded **after E9.5 up to and including E13.5**. Exactly E9.5
  is allowed. **After E13.5 is allowed**, with the source and stages disclosed.
- **Interpolation targets:** external data strictly between the bracketing released stages is excluded.
  - T2 heart: strictly between E8.25 and E8.75.
  - T2 embryo: strictly between E7.25 and E8.0.
  - The bracketing stages themselves are allowed.
- **Nearness counts.** Somite counts and Theiler stages are converted to E-days and judged by the same boundaries.
- **Genotype.** Data from the held-out knockouts, from another allele of the same gene at a comparable stage, or from a
  perturbation that phenocopies them counts as held out (beta-catenin/Ctnnb1 and Gata4 at around E8.75).
- **General resources** that span a window may be used only after removing the in-window cells; the removal must be
  described in the method summary.
- A model shared across tasks must respect the **union** of the windows of every task it serves.

## Process for adding an external source

1. Check its stages and genotypes against every window above, for every task the model serves.
2. Add a row to `EXTERNAL_DATA.md` (source, version or accession, licence, stages used, cells removed and how, which
   tasks, how it enters the method) **before** using it.
3. Keep the filtering code in the repo so the removal is reproducible.
4. If a source sits near a boundary, ask the organisers first. An answer given in advance protects us.

## Scores

Returned scores may be used to choose among predictions, not to compute them. Trying a setting and seeing its score is
fine; solving for a property of a hidden target from scores (inverting the skill map, differencing submissions) is
prohibited.

## Organiser rulings on record (Slack #virtual-embryo-challenge; details in `docs/research/slack_digest.md`)

- **2026-10-15, T3 priors (Weize Xu):** generic pathway and gene-set databases are allowed, except gene sets built
  from a knockout or perturbation of the held-out gene. **Gene-level results of losing the held-out gene in embryos at
  a comparable stage (which genes go down, by how much) are not allowed, even when read off a paper.** Broad
  background on what the phenotype looks like is allowed. For us: no gene-level beta-catenin or Gata4 knockout findings
  from the literature; phenotype-level knowledge (e.g. which compartments fail to expand) is fine.
- **2026-10-15 and 2026-09-04, T2 heart window:** the extrapolation window (after E9.5 up to and including E13.5
  excluded; after E13.5 allowed with disclosure) applies to the T2 heart targets as well as T1.
- **2026-08-18, external perturbation data far from the test stage:** GSE208162 (E12.5 endocardial Gata4 inactivation)
  approved for T3 "as it is substantially far away for the perturbation time point". External data is therefore usable
  for T3 within the windows.
- **2026-10-16, scoring conventions:** submissions are never renormalised (submit log1p CP10k); cell limits are inclusive
  and enforced at both ends.

## Open questions for the organisers

Tracked in `docs/ORGANIZER_QUESTIONS.md` until answered; record answers here with the date.
