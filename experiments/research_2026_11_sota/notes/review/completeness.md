# Review of the draft REPORT.md (SOTA sweep, 2026-11-06): completeness and recency

Reviewer: completeness and recency critic. Object:
`experiments/research_2026_11_sota/REPORT.md` (draft, section 9 pending). AS_OF 2026-11-06. Read-only; this file is
the only write. No web.

## Verdict

The bottom line holds: no 2025-2026 item in the audited notes moves our total materially under our rules, and the
corpus is recent (125 of 163 registry entries, 77%, are dated 2025-2026). Revise before publication on four points:
- the first screen does not answer the captain in his own terms (issue 2);
- the embryo s2 threshold given to the captain drops the embryo critic's calibrated break-even (issue 1);
- the largest open tail ruling, T1's sheet C, is missing from the action lists (issue 3);
- three passages state a designer's reading where the critic disagreed or the source says less (issues 4, 5, 6).

None of the coverage gaps (issue 7, item 1 below) is likely to change the bottom line. One could flip a ruling: legal
evidence on per-embryo staging spread decides the sign of the embryo s2 bet (range -1.0 to +0.8 Overall).

## Issues, most severe first

### 1. High. The embryo s2 threshold and the critics' disagreement are misstated

**Where:**
- Executive summary, "What to do" item 4: "the embryo s2 at mixture fraction 1/4: about +0.15, and only if he holds
  P(E7.75 is E7.5-like) at or above about 0.2;".
- Section 8, Recommendations item 4: "worth it only if he himself puts P(E7.75 is E7.5-like) at about 0.2 or more",
  and "it moves cell state and neighbourhoods only".
- Section 8, Confidence: "The disagreement is about the captain's own prior, which only he holds."

**What is wrong.**
- The 0.2 threshold and the +0.15 are the cross critic's figures, built on sheet B's pure-mixture toy.
- The embryo critic's -0.05 (range -1.0 to +0.8) comes from a different skill map, not a different prior. It adds
  the non-mixture residual that committed skills show (matched-mixture skill 0.7-0.8, rho 0.083-0.19). With it:
  - break-even rises to P(t' < 0.375) of about 0.3-0.6;
  - the near-E7.5 gain falls to +0.7 to +2.9 Overall.
- The embryo critic recommends adoption only if the captain states a staging spread sigma of at least 0.15 day
  (gestational labels) and accepts the low-residual end (`notes/design/embryo_critique.md` sections 1, 3 and 4 C1).
- The cross critic notes the residual too ("the winning-world gain shrinks"), yet keeps sheet B's 0.15-0.2.
- The embryo critic also corrected "cell state and neighbourhoods only". The scorer's shared reference-noise term
  makes the expression group partly hedgeable (+0.0 to +0.15 Overall; embryo critique, finding 1).

**Fix.**
- Executive summary item: "the embryo s2 at mixture fraction 1/4: between about -0.05 and +0.15. Break-even is
  P(E7.75 is E7.5-like) of about 0.15-0.2 on sheet B's pure-mixture toy, and about 0.3-0.6 once the residual that
  committed skills show is included. The embryo critic advises it only if the captain holds a per-embryo staging
  spread of 0.15 day or more and accepts the low-residual end;".
- Section 8, Confidence: replace the quoted sentence with "They disagree on the skill map, not only on the prior. The
  embryo critic adds the non-mixture residual that committed skills show (rho 0.08-0.19), which raises break-even from
  0.15-0.2 to 0.3-0.6. The staging spread is the captain's prior; the residual is a stated modelling assumption."
- Recommendation 4: replace "it moves cell state and neighbourhoods only" with "it moves cell state and
  neighbourhoods, and the expression group slightly, through the shared reference noise".

### 2. High. The first screen does not answer the captain's question in his terms

**Where:**
- Executive summary, first paragraph: "That lifts the chance of reaching our own 95th percentile from 5% to about 6-7%
  (5.9% at +0.35, 7.3% at +0.8). 209 needs about +16 Overall of sure gain to become an even bet (section 2)."
- Lines 13-19, the compliance list that sits between the title and the answer.

**What is wrong.**
- (a) Every probability in the summary is reading A's, the upper reading. The committed outlook is 184 (170-198),
  about 25 points short of 209 (`docs/PLAN.md` section 8). The summary does not say so; section 2 does.
- (b) The captain's own number, P(total >= 209), is not in the summary. It is 0.42% today, 0.59% with the package and
  0.89% with R1 and sheet C (section 2 table; `experiments/gate/results/path_to_target_2026-11-06.json`). The "+16"
  is reading A's too; on the outlook the gap is about +25.
- (c) Task h's measurement-based T3 payoffs put 209, and even reading A's q95, out of reach. That reading appears only
  as section 2's "Opposing reading".
- (d) The title promises "our top-3 chance", but the report never states one.
  - The team's established statement is "a few percent at most" (2026-11-05 report, summary; `docs/PLAN.md`
    decision log 2026-11-05 (f)).
  - 209 was the leading Human total when the gate was set (`docs/PLAN.md` section 8), a first-place figure.
  - The top-3 line moves with the field ("Rank, not level", 2026-11-05 report section 1).
  - A reader can take 0.42% as the top-3 chance, or 5-7% as a central figure.
- (e) The compliance list pushes the answer down the page.

**Fix.** Insert after the bold bottom line:

> **In the captain's terms.** On reading A (the upper reading; mean 191.5), P(total >= 209) is 0.42% today, 0.59%
> with this package and 0.89% with R1 and sheet C. An even chance needs a sure +16 on T1 or T3. On the committed
> outlook (184, range 170-198; `docs/PLAN.md` section 8) the gap is about 25. Under task h's measurement-based T3
> payoffs, 209 is out of reach. 209 was the leading Human total when the gate was set. The top-3 line moves with the
> field, and the top-3 chance stays a few percent at most (2026-11-05 report). The draws that reach 209 are our
> existing bets landing together (T3 about 75, T1 about 73), not anything found here.

Then move "No agent saw a test target ..." and its list to the end of section 1. Add the ruling deadline to "What to
do" item 4 ("before about 24 Nov, so builds can be re-hashed before the freeze", from section 11).

### 3. Medium. The largest open tail ruling (T1 sheet C) is missing from the action lists; options are ranked by mean only

**Where:**
- Executive summary, "What to do" item 4 ("Optional captain's rulings, each small").
- Section 8, Recommendations item 4.
- The section 8 tables, which have only a "Corrected gain" column.

**What is wrong.**
- Sheet C appears only inside the summary's "+0.8" sentence. It is the T1 strong-null down-tail prior of the
  2026-11-05 report (section 7, sheet C):
  - +1.5 to +2 T1 if right, and -0.2 to -0.45 if not;
  - break-even at a prior of 0.2-0.3;
  - deadline: the s1 freeze, 2026-11-27 00:00 UTC.
- The T1 critic flags it as "the one skewed T1 ruling on record", and says every T1 s1 change competes with it and
  sends T1's s2 back to the captain (`notes/design/t1_critique.md` section 3).
- By the cross design's ranking table it adds about +1.17 points of P(>= our q95) (+0.29 expected). Its critic did not
  check that figure. Even so, it exceeds every item in recommendation 4: the embryo s2 adds +0.4-0.5 points after the
  cross critic, and A15 hygiene about +0.3.
- The captain's objective is the upper tail, so the per-option tail effect is the decision quantity. Section 8 drops
  the ΔP(>= q90/q95) figures that the cross design and its critic computed.
- The cycling-quota line omits its critic's caveats:
  - its range (-1.0 to +0.4) is skewed down;
  - its sign comes from a cross-lab, depth-sensitive cycling call on Leonard's nuclei, the objection that killed A and
    B (t1 critique, D).

**Fix.**
- Add to both lists: "Still open from 2026-11-05: T1's strong-null down-tail prior (sheet C). +1.5 to +2 T1 if right,
  -0.2 to -0.45 if not; break-even prior 0.2-0.3; deadline the s1 freeze. The cross design puts it at about +1.2
  points of P(>= our q95); its critic did not check that figure. Any T1 s1 change competes with it and sends T1's s2
  back to the captain."
- Add a column "ΔP(>= our q95), points" to section 8's tables:
  - embryo s2 +0.4-0.5 (cross critic; the embryo critic puts its tail value lower);
  - A15 hygiene about +0.3 (cross critic);
  - every other row "not computed".
- Extend the cycling-quota line: "range -1.0 to +0.4; its sign rests on a cross-lab cycling call, to be measured
  depth-matched first".

### 4. Medium. Heart: the report credits a "new point" that its critic rejected, and keeps the variant the critic flagged

**Where:**
- Section 8, Heart, first bullet: "The designer's new point stands: the part of B0's momentum shared across lineages
  carries the reference run's per-gene offset, which enters the test with the wrong sign. That is the reason for H-B."
- Recommendation 4: "Hold the heart's momentum at test (or keep it for cardiomyocytes only): about +0.03."

**What is wrong.** `notes/design/heart_critique.md` says the opposite on both counts.
- The premise "cannot be identified in these files" (sections 1 and 5). The lineage-shared part also carries:
  - segmentation admixture, which tracks each stage's imaged field;
  - the monotone depth drift (233 down to 56 counts per cell);
  - shared biology;
  - at test only, the dissection change.
- H-B's mechanism "is not new". The 2026-11-05 dynamics note already put the run term in the across-break increment,
  and sheet E was priced on it. H-A repeats that note's closed L3.
- The critic offers only the full hold. The CM-only form rests on an uncommitted carry-over line that may be
  scRNA-derived, which A15's working answer keeps out of uploads.

**Fix.**
- Replace the bullet with: "The designer's premise, that the lineage-shared part of B0's momentum is the reference
  run's offset, cannot be identified in these files: that part also carries admixture, the depth drift, shared
  biology and, at test, the dissection change. It is not new: the 2026-11-05 dynamics note priced sheet E on the same
  run term, and H-A repeats that note's closed L3. H-B is sheet E, still pending, at about +0.03."
- Recommendation 4: "Hold the heart's momentum at test, in the full-hold form only (the CM-only form rests on an
  uncommitted line and may engage A15): about +0.03; default no change."

### 5. Medium. T3 phenotype: the report says the phenotype "agrees" where its note says the main part is least supported

**Where:** Section 4, "Phenotype only": "This agrees with the depletion s1 and s2 already make."

**What is wrong.**
- s1 depletes five Mesp1-lineage states: OFT/RV-CM, aSHF, pSHF, aPHM and pPHM.
- Task c (finding 15; mapping 11) reads [78] as supporting the sparing of FHF/LV states and the depletion of OFT/RV.
- It also says "depleting the aSHF/pSHF/aPHM/pPHM progenitor states is the least supported part of the scenario",
  since pharyngeal mesoderm is reported grossly normal.
- The reading rests on one 2007 paper. Related second-heart-field Cre losses were deliberately left unopened to limit
  class-B exposure (task c, Gaps).

**Fix.** Replace with: "It agrees with sparing FHF/LV states and depleting OFT/RV-CM. It gives the weakest support to
depleting the four progenitor states aSHF, pSHF, aPHM and pPHM, since pharyngeal mesoderm is reported grossly normal
(task c, mapping 11). It rests on one 2007 paper; related second-heart-field Cre losses were left unopened to limit
class-B exposure. T3's composition is closed, so this is risk information only."

### 6. Medium. Key finding 2 rests on old, off-organ sources, and the protocol argument does not apply to Leonard

**Where:**
- Section 10, finding 2: "but each changes protocol, lab or dissection by more than the effects we would measure
  [45][46][47]."
- Section 3, "The critics' main objection to building on these resources:", where "Protocol moves shares by tens of
  points" is listed first.

**What is wrong.**
- [45] is a 2022 study of human fetal hearts.
- [46] is a 2020 study of adult mouse hearts, read at search-snippet level (task b, source 16, authority 6).
- [47] is a 2024 study of lung adenocarcinoma.
- None measures a per-type temporal direction in developing mouse heart. Task b's Gaps record that no study does, so
  "by more than the effects we would measure" is unsupported.
- Leonard is in T1's own chemistry (task b, M1: "same assay as T1"), so the protocol objection does not touch it.
- The T1 critic's decisive reason is the one listed last. Leonard is a single stage from another lab, so no
  within-dataset time contrast cancels its lab and dissection offset, and the backtest reproduces its dissection
  confound (t1 critique, A.1, A.3, B.2, and section 5, "same chemistry is not same scale").

**Fix.**
- Finding 2: "New legal anchors exist (Leonard 2025, mdCLA 2026) [50][33]. Leonard is one lab's single stage, so no
  within-dataset time contrast cancels its lab and dissection offset. mdCLA changes protocol and dissection. No study
  measures how far such offsets move a per-type temporal direction (task b, Gaps). Older off-organ comparisons show
  protocol moving cell-type shares by tens of points [45][46][47] ([46] at snippet level). They bear on share priors,
  not on Leonard."
- Section 3: put the lab-offset objection first, and label [45]-[47] "(2020-2024; human fetal heart, adult mouse heart
  at snippet level, lung tumour)".

### 7. Medium. The limitations list one venue gap; the notes record more, and some bear on open rulings

**Where:**
- Section 11, Limitations, first bullet: "ICLR 2027 submissions could not be read (an OpenReview browser check), and no
  search engine indexed them yet."
- Section 1, Confidence ("medium for coverage").

**What is wrong.** The notes record further gaps that the report omits:
- The NeurIPS 2026 accepted-paper list was not scanned; only PACE is confirmed, from its repository (task a, Gaps;
  task a2, Gaps). NeurIPS 2026 is the challenge's own venue.
- ICML 2026 is named in no note's venue list.
- A legal lead was not assessed: a 2026 Cell single-cell time series of whole-embryo endothelial cells at 0.5-day
  sampling (PubMed 41702400; task a, Gaps).
- No 2025-2026 comparison of single-cell and single-nucleus data was found (task b, Gaps). That gap leaves C4a's gene
  classes resting on 2017-2024 sources.
- The per-embryo staging spread that decides the embryo s2 has no source ("sigma itself is sourced nowhere",
  embryo critique).
- No mouse in vivo knockout single-cell study at E8-E9.5 was found (task c, Gaps).

**Fix.** Add these as Limitations bullets, and add item 1's follow-ups below to "Future directions", the staging-spread
diagnostic first, because it feeds a ruling due before the freeze.

### 8. Medium. The untracked research folder also trips the clean-tree guard, and a cited process record does not exist yet

**Where:**
- Line 19: "Process notes are in `docs/DATA_POLICY.md` 2026-11-06 (a)."
- Executive summary item 1, "**Done, operational:** commit `experiments/gate/path_to_target.py` ...".
- Section 8, Recommendation 1, "(done with this report)".

**What is wrong.**
- `docs/DATA_POLICY.md` is unmodified and has no 2026-11-06 entry.
- At review time `experiments/gate/path_to_target.py` and its result are still untracked.
- So are `REPORT.md` and all of `notes/`. vec.provenance's clean-tree rule exempts only `results/`, `runs/` and
  `experiments/**/results/`, so the untracked notes folder trips the same upload guard that recommendation 1 names.

**Fix.**
- Write the 2026-11-06 (a) process entry before the report is committed, or reword line 19 to "to be recorded".
- Reword recommendation 1: "Commit the script, its result, this report and `notes/` together. The clean-tree guard
  exempts only `results/`, `runs/` and `experiments/**/results/`."
- Mark the item "done" only after the commit.

### 9. Low. Section 6 says s1 answers both organiser pitfalls, but is silent on the file that recommendation 2 would upload

**Where:** Section 6, "The T3 page": "s1 measures its factor in compartments the knockout cannot touch and at the test
stage's WT, which answers both."

**What is wrong.**
- Under the ruling of 2026-11-06 (b), a "zero" answer to the DM of recommendation 2 makes T3's s2 the barred s1. That
  file is the Mab21l2 (E9.5) factor against WT E8.75. Its credit on Gata4 comes mostly from that E9.5-against-E8.75
  offset (cross design C3, citing t3_swap_review). That offset is the stage difference the T3 page calls a pitfall.
- Task e's lever 1 is exactly this counter-claim. It prices dropping the fallback at 0 to +0.8, on the organisers'
  7.5-point P2 loss. The committed readings do not reproduce that loss: D11b scored 59.40 / 62.15, against
  52.39 / 50.61 for WT cells.
- N11 found the swap's gain sensitive to how the barred file performs (+1.38 falls to +0.42).

**Fix.** After the quoted sentence, add: "The barred s1, which a 'zero' answer would make T3's s2, does not: it
applies the E9.5 Mab21l2 factor at E8.75, and its Gata4 credit is mostly that stage offset. Task e's lever 1 prices
this risk; the committed readings do not show the organisers' 7.5-point loss, and N11 already found the swap's gain
sensitive to it (+1.38 to +0.42)."

### 10. Low. Section 3's opposing reading picks the weakest form

**Where:** Section 3, Opposing reading: "Even then, a 27-36% W2 gain at hour-scale horizons in zebrafish says little
about two days of mouse organogenesis."

**What is wrong.**
- scFM's reprogramming arm extrapolates 24-48 h (task a, finding 1 audit), a day-scale horizon, with the same 27-36%.
- Task a2 (Gaps) gives a stronger opposing reading. mmd_u, T1's largest-weight and largest-headroom metric, is computed
  in a 30-PC space fitted on the truth, which projects much per-cell noise out. Denoising may then cost little there,
  and only the variogram penalty is certain.

**Fix.** Replace the sentence with: "Even then, the gains come from dense ladders of 6-13 snapshots (hours in zebrafish
and Drosophila, 24-48 h in reprogramming), against our three stages. A second reading (task a2, Gaps): mmd_u's
truth-fitted 30-PC space projects much noise out, so denoising may cost little on T1's largest-weight metric. Against
it, the organisers warn that a model collapsing onto the mean profile 'forfeits most of the score', and an oracle
per-type shift moved mmd_u only +0.12 and +0.30 (`notes/design/t1.md` section 1, citing TRIES)."

### 11. Low. Section 7's rebuttal contradicts the summary

**Where:** Section 7, Opposing reading: "... still holds in principle, but no surviving bold option has positive
expectation."

**What is wrong.** Two bold options on record have positive expectation:
- sheet C, at +0.29 by the cross design, which the summary counts in its +0.8;
- the embryo s2 at 1/4, at +0.15 by the cross critic.

**Fix.** Replace the clause with: "and the bold options still open are small and rest on the captain's priors: T1's
sheet C (about +0.3 at a prior of 0.3) and the embryo s2 at 1/4 (-0.05 to +0.15)."

### 12. Low. The organisers' models are not below the floor on every practice board

**Where:** Executive summary: "**The organisers' own models** sit at or below the floor on the practice boards."

**What is wrong.** pseudobulk_shift scored 55.22 on embryo interpolation (section 6). Task e's spot-check corrected this
same overstatement.

**Fix.** Replace with "sit below the floor on the T1, heart-extrapolation and Gata4 practice boards (pseudobulk_shift
scored 55.22 on embryo interpolation)".

### 13. Low. Weak or summary-level sources are stated without their caveats

**Where:**
- Section 7: "On the released test split, a deterministic centred transfer scored 0.3413 against the winning entry's
  0.3401 [144]."
- Section 3, Blood: "... more than 95% enucleated by E15.5 [40];" and "... more than 94% of circulating cells at E16.5
  [39];".
- Section 5: "SEU-3D only at E8.0 [98], Harland et al. only at or before E7.25 [99]."

**What is wrong.**
- [144] ran four frozen versions in sequence on TEST, so its figure is not a single blind test (task f, deep read [9]).
  It carries the summary's "sitting close to simple pseudobulk statistics".
- [141]'s "15th of 337" comes from a search summary. [159] was read as a search summary only (authority 3). [151] is an
  unchecked community mirror. The M5 example in [158] is from a 2020 contest.
- The audit could not verify [40] (no full text) and did not check [39] (task b, spot-check 4).
- [98] was not opened and [99] was read at metadata level, with an unresolved stage list (task d, Gaps).
- Decision relevance is low in each case: the conclusions they support (keep the pairs, keep today's T1 pair, A15)
  stand on other grounds.

**Fix.**
- Add "after four frozen versions were run in sequence on that split, so not a single blind test" to the [144]
  sentence.
- Extend section 7's Confidence: "[141]'s rank and [159] were read as search summaries, [151] is an unchecked mirror,
  and the M5 example in [158] is a 2020 contest".
- Add "(not verifiable in audit)" after [40] and "(not checked)" after [39].
- Add "(metadata only; Harland's stage list unresolved)" after [99].

### 14. Low. Source accounting double-counts four works, and one cited source has no reference entry

**Where:**
- Header: "Sources: 163 registry entries from 41 domains".
- Executive summary: "with 163 sources mostly from 2025-2026".
- Section 5: "[13 in task d, cited for method only]".

**What is wrong.**
- The registry gives two numbers to each of four works: scTimeBench [1]/[26], Chreode [13]/[127], mdCLA [33]/[102]
  and MERFISH+ [104]/[126]. So 163 entries are 159 distinct works. The 2026-11-05 report reported distinct works
  ("131 registry entries, 128 distinct works").
- Raiola et al. (eLife 2026), task d's source 13, has a screened title. The registry skipped it ("removed source lines
  skipped: 1"), so it has no global number and no entry in section 12.

**Fix.**
- Header: "163 registry entries (159 distinct works)"; summary: "159 distinct sources".
- Give Raiola a global entry with its screened title and add it to section 12. Otherwise cite it as "(task d, source
  13; method only; not in the registry)".

## Item 1. State-of-the-art areas not covered, or covered only at abstract level (follow-ups; no web used)

| # | Area | Status | Board and group | Could it change a conclusion? |
|---|---|---|---|---|
| 1 | NeurIPS 2026 accepted papers (single-cell temporal forecasting, perturbation transfer, spatial generation) | Not scanned; only PACE confirmed, from its repository (task a and a2, Gaps) | T1 direction, DE and mmd_u (.25/.25/.30); T3 identity and direction; T2 | Unlikely. A new method still needs a +0.3 counted gain before 27 Nov and meets the same denoising and momentum objections. At most about +0.3 Overall. Follow-up: a screened scan of the list for persistence or nearest-stage baselines by about 20 Nov |
| 2 | ICML 2026 proceedings (July 2026) | Named in no note | T1, T2 (flow matching, OT) | As row 1 |
| 3 | ICLR 2027 submissions | Disclosed (OpenReview check) | T1, T2, T3 | As row 1; results land after the freeze |
| 4 | Legal per-embryo staging spread at stages outside every window. Candidates: mdCLA, which records each embryo's gestational day with a staging photo (task b, source 3); JAX/Qiu 2024 somite counts in 2-h bins at E8.75-E9.5, if plug age is recorded per embryo; the DMDD E14.5 sub-stage spread (read only partly, task d) | Not gathered; the embryo critic says sigma is "sourced nowhere" | Embryo: cell state and neighbourhoods (50 of 100 board points) | **Yes, for one ruling.** It is the input that decides the embryo s2. At sigma 0.05 day the bet loses in every model; at 0.15 it is +0.33 to -0.02 Overall on the calibrated maps (+1.50 on the toy; embryo critique, section 3); range -1.0 to +0.8 Overall. Follow-up: a target-free diagnostic before the ruling, on legal stages only (remove in-window embryos by code), with the source registered first (an X003 amendment for JAX) and the staging rule checked for the DMDD sub-stages |
| 5 | 2026 Cell whole-embryo endothelial single-cell time series, 0.5-day sampling (PubMed 41702400) | "Not assessed" (task a, Gaps) | T1 direction and DE for endothelium, one of the four Feng-directed types | Possibly small. A within-dataset long-range bracket (E9.5 or earlier against E14.25 or later) is the design the T1 critic says cancels lab offsets. Ruling-only, at most about +0.1-0.2 Overall. Follow-up: a screened read of its stages, chemistry and licence; in-window removal by our code; a register row before any use |
| 6 | 2025-2026 comparisons of single-cell and single-nucleus data (per-type temporal fold changes) in developing tissue | None found (task b, Gaps); the report relies on 2017-2024 sources | T1: C4a hygiene (+0.10, -0.3 to +0.4); Key finding 2 | It could flip the hygiene bet's sign; it would not reopen T1's A or B |
| 7 | Mouse in vivo knockout single-cell studies near E8-E9.5, including pre-2025 multi-knockout embryo atlases | Not found (task c, Gaps); the pre-2025 sets were not searched | T3: how often an embryo knockout response carries a generic delay or stress component shared across genes (p_anchor) | Evidence only (scRNA; A15 bars it as an input). It could shift the p_anchor pricing, not an upload before the freeze |
| 8 | Imaging-ST segmentation and molecule-misassignment correction, 2025-2026 (e.g. Proseg, resolVI; from my own knowledge, unverified) | Read only through [70] | Heart (admixture is why H-A's premise is unidentifiable); T3 (spillover; C4) | Unlikely. Correcting our cells moves them away from an uncorrected truth; about 0 |
| 9 | VCC 2026 winners (mid-late November), the Cell 2026 "lessons" commentary (403), Myllia and VCPI results | Task f, Gaps | Section 7's lessons; T3 transfer | They land at or after the freeze; about 0 for uploads (already under "Watch") |
| 10 | Organiser Slack, Discord and X | Not searched by this research (task e, Gaps); the last screened sweep was 2026-11-05 (j) | Rules: A20, A24, A15 | Only through an answer, which the runbook handles. Section 6's "The scorer and the rules have not changed" should name its scope: Rules, FAQ, Timeline and veckit read on 2026-11-06, Slack swept on 2026-11-05 |

## Item 2. Stale or weak sources in decision-relevant claims

- **Issue:** Key finding 2 and section 3 ([45] 2022, [46] 2020 at snippet level, [47] 2024; off-organ). See issue 6.
- **Issue:** [144]'s TEST figure, the summary-level [141] and [159], the unchecked mirror [151], the 2020 M5 example
  in [158], the unverified blood anchors [39] and [40], and the metadata-only [98] and [99]. See issue 13.
- **Issue:** [78] (2007) is the only phenotype source and is read beyond what task c supports. See issue 5.
- **No issue:**
  - [2] scFM is one unreviewed preprint, but the report marks its venue as unconfirmed. Its conclusion ("no
    naive-baseline win at our metrics") would stand without it.
  - The template figures from [55] are one preprint, but [53], [65], [66] and [72] corroborate the structure.
  - The corpus is recent (77% from 2025-2026).
  - No other pre-2025 source is presented as current state of the art.
  - The organiser record is current: Rules, FAQ, Timeline and veckit were read on 2026-11-06.

## Item 3. Opposing readings, section by section

| Section | Given | Strongest? |
|---|---|---|
| 2 | Task h's payoffs | Strong. It belongs in the summary (issue 2). Missing: rank, not level (the top-3 line moves with the field), and the committed-outlook reading (184) |
| 3 | scFM sampling with noise | Weakened by "hour-scale"; misses task a2's truth-PC reading (issue 10) |
| 4 | Cell lines against two embryo runs | Adequate; no issue. The phenotype over-read is in the body text (issue 5) |
| 5 | A learned heart model | Covers the heart only. The embryo half's opposing reading (E7.75 away from midway) sits in section 8; a pointer would do |
| 6 | Births may win at test | Plausible. Misses the stage-difference pitfall against the barred s1 (issue 9) |
| 7 | Kaggle-scale picks against two uploads | The right reading; the rebuttal contradicts the summary (issue 11) |
| 8 | None (Confidence only) | Missing. The embryo critic's calibrated break-even (issue 1) and the heart critic's verdict (issue 4) are the strongest counter-claims |

## Item 4. Structure

- **Bottom line.** It is plain and on the first screen ("holds no method, dataset or tactic that moves our total
  materially"). No issue.
- **Missing from the first screen.** The captain's own number (P >= 209), the reading-A caveat, task h's reading, and
  the top-3 chance in the team's terms (issue 2).
- **Buried:**
  - sheet C, the largest open tail ruling (issue 3);
  - the embryo critic's break-even (issue 1);
  - the per-option tail effects that the objective needs (issue 3).
- **Order.** The compliance list precedes the answer; move it to section 1 (issue 2).

## Checks that passed (no issue)

- Section 2's figures match `experiments/gate/results/path_to_target_2026-11-06.json`, which is untracked at review
  time. Checked: mean 191.5; q90/q95/q99 of 203.95/205.54/208.07; P(>= 209) 0.42%; the +0.35, +0.8, +1, +3, +5 and +10
  rows; +16 for an even chance; the board means given q95.
- The T1, T3 and cross-board verdict tables match their critiques. The cross critic's kill of front-loading is carried
  into the summary and recommendation 3.
- 209 is reported for information only and in a table of its own, as `docs/PLAN.md` section 8 and the cross critic
  ask.
- The ICLR 2027 gap is disclosed, and the VCC 2026 results are on the watch list.

## Process

- **Read:**
  - REPORT.md;
  - `notes/brief.md`, `notes/registry.md`, `notes/task-a.md`, `task-a2.md` and `task-b.md` to `task-f.md`;
  - all ten files in `notes/design/`;
  - `docs/PLAN.md` sections 7 and 8 and the decision log of 2026-11-04 (e) to 2026-11-06 (b);
  - `docs/DATA_POLICY.md`: Scores, and a grep for 2026-11-06 entries;
  - `docs/CHALLENGE.md`: the prize lines;
  - `docs/ORGANIZER_QUESTIONS.md`: sweep dates;
  - `experiments/research_2026_11/REPORT.md`: front matter, sections 1 and 7;
  - `docs/research/digest_2026-10-19.md`: sections 2.3, 3 and 5 (screened);
  - the clean-tree lines of `src/vec/submit.py` and `src/vec/provenance.py`;
  - `experiments/gate/results/path_to_target_2026-11-06.json`.
- **Commands:** git reads with `GIT_OPTIONAL_LOCKS=0` (log, status, ls-files). One `python3 -I -B` stdlib read of the
  JSON and one awk count of the registry dates. No repository code ran.
- **Not opened:**
  - `boards/`, `data/`, `runs/`, `docs/research/` (apart from the screened digest);
  - the workflow's `raw/` and `screened/` folders;
  - anything under `~/.claude-account` or `~/.claude`.
  - The harness saved one oversized registry listing under `~/.claude-account/.../tool-results/`. It was not opened;
    the listing was redone in smaller pieces.
  - The harness placed the auto-memory index in this agent's context. It was not opened or used.
- **Barred classes:** none met. No per-board leaderboard figure or hidden-target property is used, and no leaderboard
  total other than the captain's 209 is used.
