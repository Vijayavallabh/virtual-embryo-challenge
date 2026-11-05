# Data policy

Condensed from the Official Rules s.10 and the Data page (https://virtualembryo.ai/challenge/rules,
https://virtualembryo.ai/challenge/data), as read on 2026-10-18. When in doubt, the live rules win; re-read them before
each submission. Any change to this file needs a dated note at the bottom.

## Competition data

- Released for research and education within the Challenge. Never commit it, share it with other teams or post it.
- Unpublished (Qingquan Zhang and Neil Chi, UCSD). Any publication or public presentation needs prior written approval
  from Dr. Neil Chi's group and must credit Dr. Qingquan Zhang and Dr. Neil Chi, plus acknowledge the Virtual Embryo
  Challenge (Terms of Use clause 14).
- Validation answers are served in the final phase (opened 2026-10-27, postponed from 2026-10-20) only after a team
  member unlocks them on the Data page, which permanently ends P2 submissions for the whole team; the user unlocked
  them for our team on 2026-10-28 at 09:32 UTC. From then they may be used "solely for preparing P3 predictions"
  (Rules s.10; the exception below the held-out table). Until then they are held out like any other target. Test
  ground truth is never released.
- P3 methods, prediction files and detailed evaluation feedback stay private from other participants until the
  2026-12-11 final score announcement (organisers' email of 2026-10-20 22:03 IST; partner rules s.5, in Yukon's section
  but not limited to Yukon by its wording). Ours stay in this private repo and the DGX workspace; nothing about them
  goes into Slack or any other public place before then. Official P3 scores themselves are published when each
  submission is made (Rules s.11).

## Held-out stages and genotypes (no measured data from these, by any route)

| Task | Held out |
|---|---|
| T1 (scRNA) | E10.5 (val), E12.5 (test) |
| T2 embryo | E7.5 (val), E7.75 (test) |
| T2 heart | E8.5 (val interp), E10.5 (val extrap), E12.5 (test) |
| T3 | Gata4 KO and beta-catenin KO at E8.75 |

"By any route" includes pretrained models whose training data contains the window, and public atlases that contain it.

**Final phase, 2026-10-28.** The official validation files (T1 E10.5, T2 E8.5, E10.5 and E7.5, the two T3 Gata4
replicates), unlocked by our team on the Data page, are the one exception (Rules s.10): "solely for preparing P3
predictions", and the Data page: "you may use those released validation files for P3 preparation. This does not
permit access to hidden test data or other restricted external sources." The user decided on 2026-10-28 that this
covers using them as model inputs when producing the test predictions (A20 passed; organiser confirmation asked in the
fourth email). The exception covers those files only: external data from the validation stages and genotypes is still
prohibited, and the test windows (T1 and T2 heart E12.5, T2 embryo E7.75, the beta-catenin KO) are unchanged.

## External data windows

- **Extrapolation targets (T1, T2 heart):** external data excluded **after E9.5 up to and including E13.5**. Exactly E9.5
  is allowed. **After E13.5 is allowed**, with the source and stages disclosed.
- **Interpolation targets:** external data strictly between the bracketing released stages is excluded.
  - T2 heart: strictly between E8.25 and E8.75.
  - T2 embryo: strictly between E7.25 and E8.0.
  - The bracketing stages themselves are allowed.
- **Nearness counts.** Somite counts and Theiler stages are converted to E-days and judged by the same boundaries.
  Rules s.10 gives the general rule, read word for word by a screening agent on 2026-11-05: "External data is
  treated as a held-out stage when it sits closer to that stage than to the nearest stage we released for the same
  task: the midpoint between them is the boundary." For the extrapolation targets, the same section says the window
  "is stated outright rather than derived, and it relaxes the midpoint rule above rather than following it". So the
  explicit window governs there:
  - "Anything after E13.5 may be used, provided the source and its stages are stated explicitly with the
    submission."
  - Data at exactly E13.5 may not be used, and data at exactly E9.5 may.
  No source the method uses changes status.
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

A screened read of the leaderboard hands the main agent overall totals and overall ranks only. Per-board or per-task
figures of any team (scores, ranges, ranks, gaps, correlations) stay with the screening record and out of every file an
agent reads, the screened digest included (2026-11-05 (h), (i)).

Every screening brief quotes this list. Removed: test-board limits, anchors and truth counts, and any other property of
a hidden target; per-board leaderboard figures; gene-level knockout findings from outside the released data (the ruling
of 2026-10-15 below). Kept: overall totals and ranks; phenotype-level knockout knowledge, which that ruling allows; our
own results; the organisers' published baselines.

## Do-not-use list (known traps)

These leak properties of held-out targets. Do not read them into any model, prior or hyperparameter unless the
organisers approve in writing.

- **The PyPI veckit 0.1.2 docstring quotes an RMS radius for the held-out E12.5 heart stage.** That is a measured property
  of a test target. Use veckit from git `main` (pinned in `pyproject.toml`) and never copy numbers from 0.1.2 docs.
- **Target sizes or other properties recovered from scores** (inverting the skill map with published anchors, differencing
  submissions). At least one public community repo did this; it is prohibited by s.10.
- **Test-board anchors in the final-phase index.json** (floor and ceiling metric values). A floor value is the distance
  between the reference stage and the hidden truth, so it encodes the hidden target (for scale_log_ratio, its size; for
  mmd_u on T3, how far the knockout is from wild type). The intake strips them before anyone looks; they are never read,
  printed, committed or passed to the emulator. Validation-board anchors may be used to calibrate the emulator once
  the validation answers are released (for our team, the unlock of 2026-10-28). Before then they describe held-out
  targets, and nothing is calibrated, tuned or checked against them (the scorer-build anchor check refuses to read
  them before the final phase).
- **Our own size forecasts from incomplete stage sets** (dry runs of the size procedure before the release) are code
  checks, not predictions, and stay out of the docs.
- **The public board file and board endpoints themselves.** `panels/index.json` on the site (and the public
  `/challenge/boards` endpoint) may carry test-board anchors at any time now (the final phase opened on 2026-10-27;
  since then the file answers HTTP 403 to `vec boards update`). Nobody opens them in a browser or fetches them by
  hand or by script; boards are refreshed only through `vec boards update`, which strips non-validation anchors in
  memory. Its 403 is not worked around (the user's rule: "no workarounds for a 403 (no other user agents, no other
  endpoints)"), and the fourth organiser email (2026-10-28) told the organisers we had not tried to get around it and
  asked whether the block is intended. At the user's request ("build the intake", 2026-10-28) a second intake exists
  for a copy the user would save with the browser's "Save link as", unopened, outside the repository:
  `vec boards update --from-file PATH` (`vec.boards.update_boards_from_file`). It strips the anchors in memory the
  same way, refuses a file inside any git work tree before reading it (on this machine the repository's parent folder
  is one), and with `--delete-source` deletes the saved file afterwards (otherwise it says to). **It is not cleared
  for use.** Saving the file with a browser is a fetch by hand and a workaround for the 403, so the route may be used
  only after the user lifts both rules for it, recorded here with the date and the user's words. If it is ever used,
  nobody else reads the saved file. 2026-10-31 (c): told to act on the organisers' most rational answers, AI agents
  predicted no objection (0.8; the 403 is Cloudflare's browser-signature block, error code 1010); the judge and the red
  team held that this general instruction does not by itself lift the two rules, so the route waits for the captain's
  explicit yes, then a dated notice on the fourth-email thread and about 72 hours without objection before the save.
  Test-board truth_cells, min_cells and max_cells can encode a held-out stage's cell count: code reads them from
  `boards/`, nobody quotes them. **2026-11-01: the captain lifted both rules for this one route** ("yes for all",
  answering the plan: the dated notice first, the save not before 2026-11-04 UTC and not if an objection or a fix
  arrives, the link checked by hovering, "Save link as" into a private folder, the intake, the folder deleted on any
  refusal). The two rules stand for every other route: no other user agent, endpoint or script for the board file.
  Team members skip Slack threads that discuss sizes or other values worked out from anchors, and agents reading Slack
  must not quote such threads.
- **The organisers' Virtual Embryo atlas** (virtualembryo.ai, its predictor and API) lists stages inside protected
  windows (e.g. E7.5 and E11.5 spatial data). The rules name it explicitly: no held-out-stage data from it.
- **Public atlases spanning a window** (MOCA E9.5-E13.5, MOSTA E10.5-E13.5, the JAX time-lapse between E9.5 and E13.5)
  until the in-window cells are removed by code in this repo.
- **GSE208162** (E12.5 endocardial Gata4 KO) is approved for T3 only; it is inside the T1/T2-heart window.
- **Gene-level beta-catenin or Gata4 knockout findings from papers** (see the 2026-10-15 ruling below).
- **Code from public repos without a licence** (e.g. tubager/vec-t2, wqty123): ideas only, no copying.
- **The Yukon partner team's material** (Eigen Labs' collaborative reference team, announced 2026-10-20 and detailed
  in the organisers' email of 2026-10-20 22:03 IST: shown on the leaderboard as an unranked reference, not eligible
  for the Challenge's official prizes; its authorised members share executable models, fitted weights, methods and
  evaluation feedback with each other). Because many people can join it, anything the Yukon team produces
  (discussions, code, models, weights, predictions, posts about them) may carry anchor-derived values, leaked numbers
  or banned knockout findings. Agents do not open yukon.org, the Yukon interest form or any Yukon team material;
  agents reading Slack or email only note that such material was posted, without opening or quoting it. A team member
  who looks at what Yukon offers should skip its challenge discussions, code, models and predictions, as with anchor
  threads in Slack. Nothing from the Yukon team enters a method, and no number, threshold or submission choice comes
  from its leaderboard scores. The rules agree: the partner exceptions "apply to the partner category only", Yukon
  members share "within the team", and Rules s.4 bars attempts to access another team's work. Moving from an official
  team to the Yukon team is one-way for whoever moves (no return to an ordinary team this edition, no joining or
  creating one, even after leaving Yukon; earlier submissions stay with the original team) and ends that person's
  eligibility for the official prizes; whether any of us moves is the user's decision alone. Our team's code and
  methods are not any one member's to contribute there (code sharing across teams is prohibited, and the email limits
  contributions to work the person is authorised to share).

### Redaction record (gene-level knockout findings)

- 2026-10-18: the research notes copied into the repo contained gene-level findings from beta-catenin and Gata4
  loss-of-function papers (direction tables and paper titles stating them), which the 2026-10-15 ruling bars. They
  were removed from `docs/research/t3_biology.md`, `architecture.md`, `t3_methods.md` and `external_data.md`, and the
  session scratch copies were overwritten with the redacted versions. Git history was then rewritten (user's
  decision, 2026-10-18) so that no commit contains the original text: main force-pushed as 9fed750, old objects
  purged locally, the DGX mirror's .git refreshed. The GitHub repository was then deleted and recreated from the clean
  history (same name, private), so none of the old commits exist in it; the API returns "No commit found" for the old
  IDs. A local backup bundle of the clean history was written before the deletion.
- Mitigation: T3 gene directions come only from the pre-registered provenance rule in `docs/PLAN.md` section 4.4
  (released data and committed code; scenarios in `configs/t3_scenarios.yaml` name cell states, never genes).

### Redaction record (held-out sizes, 2026-10-19)

- A second redaction pass, by an agent that reported locations but not values, removed 132 spans from
  `docs/research/t2_spatial.md`, `community_digest.md`, `metrics.md`, `docs/CHALLENGE.md` and the size-procedure dry-run
  outputs in `experiments/coords/results/`. Each span is now `[redacted 2026-10-19: A]`. Removed: sizes and growth rates
  of the validation stages (E8.5, E10.5, E7.5) back-calculated from published anchors and community score inversions;
  veckit PR #6 radii for released stages; two figures derived from the veckit 0.1.2 docstring table, which includes
  the E12.5 test stage (see the correction in the leak record below); an in-window literature measurement (heart wall
  thickness at a Theiler stage inside the excluded window); an organiser-site volume factor for interpolation targets
  (removed for caution); and our own dry-run forecasts for E12.5 and E7.75.
- A separate screening agent reads the research workflow's raw results and hands the lead only a cleaned digest
  (`docs/research/` integration uses that digest).
- Exposure: the research agents, the coords experiment agent and its reviewer read `t2_spatial.md`. The lead saw a few
  validation-stage values (E10.5, E7.5) in the reviewer's redaction list. Nobody building the size procedure saw a
  measured or leaked E12.5 or E7.75 value.
- Mitigation: the heart and embryo size procedures contain no size numbers and choose their rules with a rate-free
  argument fixed before the release (`docs/PLAN.md` section 4.3). The measured validation stages released on
  2026-10-20 supersede anything anchor-derived. Test-board anchors are stripped at intake (do-not-use list).
- Git history: rewritten on 2026-10-19 (user's decision) so that every commit carries the redacted versions; the
  GitHub repository was deleted and recreated (private) from the clean history, and the DGX .git refreshed. Commit
  IDs quoted in docs and result files before this date map to the new history as follows (old -> new):
  - `810dc11` -> `eb0c4db`
  - `a419024` -> `67ff77f`
  - `49eef7c` -> `bdc7dd1`
  - `00c14ba` -> `f496896`
  - `ac5a19a` -> `a13552e`
  - `e3b5d2c` -> `049ae66`
  - `5078c82` -> `b768a38`
  - `1b64676` -> `83e85e7`
  - `18bb455` -> `33d6468`
  - `3ac797e` -> `a2f3750`
  - `9fed750` -> `042e595`
  - `eddbbee` -> `428b654`
  - `f87cf53` -> `41e8518`
  - `dbb6eaf` -> `4d9048b`
  - `da5692d` -> `5a9b2a6`
  - `d8c0d2a` -> `977af14`
  - `c511e4c` -> `59d266e`
  - `171695a` -> `daca416`
  - `9f83c59` -> `40c9664`
  - `5044494` -> `7c8941a`
- Git history rewritten again on 2026-10-27 (user's decision) to drop the AI co-author trailer from every commit
  message; trees and dates unchanged. Old -> new commit IDs: `docs/COMMIT_MAP.md`.
- Related findings from the screened research digest: one spacing bound in `t2_spatial.md` matched a held-out column of
  the leaked table (since redacted); no spacing, size or density parameter in the plan or code comes from it (the
  jitter and count rules use spacing measured on released files, `experiments/density/REPORT.md`). One research agent
  printed per-stage QC medians (library size, mitochondrial fraction) for all stages of an external heart dataset,
  including in-window stages, once; nothing was set from them.
- Clean-up (user's decision, 2026-10-19): deleted the session's copies of the veckit 0.1.2 wheel and sdist, the PR #6
  files and PR bodies, and the anchor-derived research outputs (scratchpad `gap_units/`, DGX `~/v/vec-gap-units`). On the
  DGX also deleted the raw external heart object and the two full metadata tables that still held E10.5-E13.5 cells,
  and an atlas cell table covering in-window embryos; only window-filtered subsets remain. External data will be
  re-ingested through committed filter code (feat-007).

### Leak handling record (veckit 0.1.2 held-out radii)

- 2026-10-18: research agents reading the scorer's PyPI release saw the E12.5 heart radius in
  `common/shape_metrics.py` (line ~245 of the 0.1.2 wheel; also in open veckit PR #6). The numbers were redacted
  from every research note before anything entered this repo (`git grep` for them is clean).
- Correction, 2026-10-19: that check covered the radius values only. Two figures derived from the same docstring table
  (a nearest-neighbour spread and a radius swing factor across its stages, which include E12.5) were still in
  `docs/research/community_digest.md`, committed on 2026-10-18. They were redacted on 2026-10-19 (record above).
- Mitigation: the E12.5 heart size comes only from the pre-registered procedure in `docs/PLAN.md` (T2 heart), fitted on
  released training stages by committed code, with its output recorded in git before any test submission. Nobody
  building the size model opens veckit 0.1.2 or PR #6.
- Disclosure: private report sent to the organisers (Weize Xu, Zehua Zeng) in Slack group DM `C0CA7QTEMDF` on
  2026-10-18 22:39 IST. Reply: _pending_ (no reply or reaction as of 2026-10-20 13:49 IST; no organiser has posted on
  Slack since 2026-10-16). Still none at the screened Slack sweep of 2026-10-23, the last one (no organiser post or
  reply anywhere).

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

- **Pending public question** (a participant, #virtual-embryo-challenge, 2026-10-19 19:41 IST): may released validation
  answers be used as training data and model inputs for the test predictions? Rules s.10 still lists the validation
  stages as held out, while the Data page says they become training data in the final phase. Our plan follows the Data
  page (assumption A20); if a ruling says otherwise, T1 and T2 heart are rebuilt from E9.5, T2 embryo from E7.25 and
  E8.0 only, and T3 uses the barred branch. No answer as of 2026-10-20 22:43 IST; the Rules and Data pages were
  unchanged on this point when checked at 17:10 UTC that day. We asked the organisers the same privately (question 2
  of the third email, sent 2026-10-21 01:24 IST; `docs/ORGANIZER_QUESTIONS.md`); no answer yet.
  2026-10-28: the final-phase Rules s.10 make the unlocked validation files an exception "solely for preparing P3
  predictions" (held-out section above), and the user read this as allowing them as model inputs too (A20 passed).
  The fourth email (sent 2026-10-28 15:41 IST) asks the organisers to confirm. No reply to any of the four emails as
  of then, and no public answer at the last Slack sweep (2026-10-23).
- **2026-10-31: no reply to any email, and no organiser answer to any team's rule question on Slack since 2026-10-16.**
  At the captain's instruction the team proceeds on predicted answers, recorded as working answers (never rulings) in
  `docs/ORGANIZER_QUESTIONS.md`. A real answer replaces a prediction the day it arrives.

## Change notes

- 2026-10-20: the organisers postponed the final phase by one week to 2026-10-27 (email to registered teams,
  2026-10-20 07:41 IST; scoring, submission limits, the 2026-12-02 deadline and the 2026-12-11 announcement unchanged).
  Release dates in this file now read 2026-10-27; records written earlier ("released on 2026-10-20") mean the same
  release. The validation answers stay held out until the final phase opens, so two release gates moved with it:
  `vec boards update` writes the T1 composition file with E10.5 in it only once the public phase endpoint reports the
  test split, and the scorer-build anchor check refuses to run before 2026-10-27. Added the Yukon partner team to the
  do-not-use list.
- 2026-10-20 (b): after a screened sweep of Slack and the organiser email: the anchor check also waits for the public
  phase endpoint; validation-board anchors stay unused until the release; the leak-report status and the pending public
  question on validation answers recorded.
- 2026-10-20 (c): the organisers' email of 2026-10-20 22:03 IST (collaborative Yukon track) and the new partner rules
  page: P3 methods, prediction files and evaluation feedback stay private until 2026-12-11 (added under Competition
  data); the Yukon entry now covers the collaborative track, the interest form, the rules behind the quarantine and
  the rule that our team's work is not one member's to contribute. Screened Slack and mail check to 22:43 IST: no
  organiser post and no reply to any of our messages.
- 2026-10-20 (d): while drafting `docs/RELEASE_RULES_SUMMARY.md`, an agent's repo-wide search printed one line each of
  `docs/research/t3_methods.md` and the screened digest (its exclusion filter failed). Neither line held a value or a
  gene-level knockout finding, and nothing from them entered the summary or any method.
- 2026-10-27: the history rewrite of that day (AI co-author trailer dropped from every commit message; old -> new IDs
  in `docs/COMMIT_MAP.md`) added to the held-out sizes record.
- 2026-10-28: the validation answers turned out to be opt-in: a team member unlocks them on the Data page, which
  permanently ends P2 submissions for the whole team (the user unlocked them at 09:32 UTC). Added the final-phase
  exception below the held-out table (Rules s.10 and the Data page; A20 passed on the user's reading, organiser
  confirmation asked in the fourth email). Updated the Competition data line on validation answers, the release date
  in the calibration condition and the board-file entry in the do-not-use list (panel files HTTP 403 since 2026-10-27),
  the leak-report status and the open question on validation answers. Records that say the validation stages were
  released on 2026-10-20 or 2026-10-27 mean, for our team, this unlock.
- 2026-10-28 (b): at the user's request ("build the intake"), the board-file entry names a second route for while the
  panel files answer HTTP 403 to scripts: a board index the user saves with the browser, never opened, goes through
  `vec boards update --from-file`, which strips test-board anchors in memory as the online refresh does. Its tests
  feed a sentinel anchor value through the success path and every refusal path and check that it never reaches the
  terminal or `boards/`. Using the route is the user's call. Before it was committed, an independent review by three
  AI agents found ways a mistaken file in the saved folder could carry anchors into `boards/` (a copy of the index
  saved as the composition or as a gene list, a kept field holding nested values); both routes now refuse those,
  each with a test, and the online refresh had the composition gap too.
- 2026-10-31: a side agent pointed out that saving the board file with a browser breaks two standing rules (nobody
  fetches it by hand; the user's "no workarounds for a 403"). The user had asked for the intake to be built, which
  is not the same as lifting those rules, and the fourth organiser email says we have not tried to get around the
  block. The board-file entry now says the saved-file route is not cleared for use until the user lifts both rules.
  The same day the online refresh still answered 403 (16:52 UTC), and none of the four organiser emails had a reply.
- 2026-10-31 (b), exposure: searching the repository for the organisers' reference rows (feat-005), the AI assistant
  ran two recursive greps meant to skip `docs/research/`, but the filter did not match the output paths (this shell's
  grep is ugrep, which prints paths without a leading "./"), so about 30 lines from `docs/research/` (landscape.md,
  t1_methods.md, t3_methods.md, engineering.md, t2_spatial.md and README.md) reached the main agent. The assistant
  checked each line: published validation-board reference rows and notes on why the organisers' baselines score as
  they do, none of the banned kinds (held-out sizes, test anchors, values from veckit 0.1.2 or PR #6, gene-level
  knockout findings from papers, Yukon material). No decision depends on them: the pre-registered step-4 outputs,
  the T3 plan among them, were committed on 2026-10-28, and the calibration verdict follows the pre-registered A2
  test. The reference rows the calibration uses came through a screening agent
  (`experiments/calibration/reference_rows.json`); searches now exclude the folder (`--exclude-dir=research`).
- 2026-10-31 (c): working answers. The captain told the AI assistant "Now assume that they don't reply to the emails and you then come up with their rational responses and act on it and proceed".
  A screened Slack sweep (16 to 31 Oct: no organiser reply to any rule question, no 403 report from any team, one
  public post linking the Data Protocol page) and a screened read of that page (released-stage facts in
  `docs/DATA_NOTES.md`; 24 items about the held-out stages withheld by the screen) fed seven AI agents that predicted
  the organisers' answers (`docs/ORGANIZER_QUESTIONS.md`). What changes: the T3 build selects the allowed branch;
  nothing else in the data rules changes. Still the captain's: the browser route above, and the red team's proposal
  to make T3's s2 the refused branch's s1.
- 2026-11-01: the captain answered "yes for all" to the open proposals. The browser route for the test-board specs is
  cleared for that one route with the conditions in the board-file entry; T3's s2 becomes the refused branch's s1, so
  one T3 submission carries no number derived from the Gata4 files (a hedge in case a ruling restricts them); and
  model changes are adopted only under the adoption rule (`docs/PLAN.md` section 3, item 8).
- 2026-11-01 (b): the captain allowed the screen to release the Data Protocol page's qualitative description of how
  the Task 1 E12.5 samples were collected (what was dissected and pooled, any processing that removes cell types),
  with every number replaced; collection protocols are not measurements, and the organisers published them. The text
  about the Task 2 heart E12.5 and embryo E7.75 stacks and the test knockout stays withheld, because the T2 size and
  scope choices are pre-registered. The same day the captain approved how the adoption rule counts pseudo-tasks
  (`docs/PLAN.md` section 3, item 8).

- 2026-11-04, exposure: the board index displayed and pasted. At 06:56 UTC `vec boards update` still got the 403
  (Cloudflare 1010). Asked which API was blocked, the AI assistant named the board index's address without repeating
  that it must not be opened; the captain then opened that address in a browser, which displays the file, and pasted
  all of it into the AI assistant's session (about 07:05 UTC). It held each board's key, label, gene-list name and
  published hash, the validation boards' full specs (as already in `boards/`), and for the four test boards no anchors
  and no truth_cells, but their min_cells and max_cells. So the AI assistant has read the test boards' cell limits,
  which the board-file entry says can encode a held-out stage's cell count; read against the validation boards' own
  limits, some of them look as if they do. No value or inference from the paste is written anywhere (repository,
  notes, commit messages, prompts to other agents, session summaries) and none is used: every use of a test board's
  limits is code fixed before today (the count rule through approved rule 7's committed decision, the heart and
  embryo sizes from the procedure pre-registered at 171695a before the release, the s2 rules frozen and hash-pinned on
  2026-11-03, the pre-registered T3 rules), and the builders never print a test-board count. From now on no try
  changes a cell number, a size or a scope on any board, or T1's composition, unless the captain rules otherwise after
  reading this entry (decided by the AI assistant under the captain's delegation of 2026-11-02). The paste stays in
  Claude Code's local log of this session, outside every repository; deleting it is the captain's call, and searches
  of that log stop before 07:00 UTC on 2026-11-04. The route is otherwise unchanged: the captain saves the index and
  the embryo test board's gene list (the one list not already in `boards/` under its published hash) with "Save link
  as" into `~/vec-panels/`; the intake reads them only after the screened Slack and email sweep started at about
  07:00 UTC finds no objection and no fix; nobody opens the saved files.
- 2026-11-04 (b), exposure: the captain then pasted the text of the Data page (as seen signed in) into the AI
  assistant's session, to show where its links are. Besides what the notes already hold, the page shows every board's
  cell limits in its panel table (the test boards' are the index's), says how the upper limit is set (read with the
  limits, that rule can give a hidden truth half's size away), and lists hidden targets' file sizes: the T2 heart
  E12.5 file, which the redaction of 2026-10-19 removed from `docs/CHALLENGE.md` (category A), and the two
  beta-catenin replicates (the T1 E12.5 file size was already in the notes). So every visit to the Data page shows
  these to the captain; the route's step 1 now says never to copy the page's text to anyone, the AI assistant
  included. The containment of the entry above is widened to T3: no try changes a cell number, a size, a scope, or
  T1's or T3's composition, unless the captain rules otherwise. Nothing from the page is written anywhere else; the
  organiser follow-up (drafted by the AI assistant, sent by the captain) names the limit rule and the listed file
  sizes as possible leaks, with no value.
- 2026-11-04 (c), the captain's ruling on the block: a side agent pointed out that the block closes T1 births (the
  plan's ~1.5-5 points); asked whether to keep it or narrow it so births by the legal route stay open (the AI
  assistant's recommendation), the captain answered "Go ahead with whatever maximizes our chances of winning this
  challenge". The AI assistant narrowed it as recommended: a T1 try may add births by the legal route (cell states
  present in legal post-E14.5 chamber data and absent at E10.5, their states and shares taken only from E10.5 and
  those data by a rule written before any run). Everything else in the block stands: cell numbers on every board,
  sizes and scope on every board, T1's other composition (the cardiomyocyte and blood shares stay as the s1 and s2
  rules set them) and T3's composition. Why this is safe: what the exposure showed for T1 is a cell count of the
  scorer's E12.5 working copy, which says nothing about which states appear or their shares; and the T1 E12.5 file
  size, a proxy for the same count, has been in `docs/CHALLENGE.md` since 2026-10-18 without redaction, because Task 1
  scores no size. A births try still goes through the tries ledger, registered before its run.
- 2026-11-04 (d), decided by the AI assistant under the captain's ruling (c) ("Go ahead with whatever maximizes our
  chances of winning this challenge") and the delegation of 2026-11-02 (the captain can overrule): T1's composition
  reopens within the fixed chamber scope by the legal route. A T1 try may change the shares of the chamber categories
  and of the states within them, births included, taken only from the released stages, the unlocked E10.5 and legal
  external data (nothing after E9.5 up to and including E13.5), by a rule committed in its ledger row before any run.
  Why this is safe: (c)'s reason covers every share, not only births. What the exposure showed for T1 is a cell count
  of the scorer's E12.5 working copy, and the T1 E12.5 file size, a proxy for the same count, has been in
  `docs/CHALLENGE.md` since 2026-10-18; a count says nothing about which states appear or in what shares. No exposed
  value is in any file, prompt or summary, and since this session's context was summarised the AI assistant's working
  context holds none. Why now: births were pruned the same day (`experiments/harness/TRIES.md`), leaving no open lever
  on T1's weakest groups, and a screened read of the public leaderboards (board and overall totals only) put a Human
  top-3 place near 202-210 overall, against about 165 from our practice-task results, with
  [per-board leaderboard figure removed 2026-11-05]. Still closed: cell numbers, sizes and scope on every board, and
  T3's composition (also fixed by the frozen T3 rules). The frozen s2 rules (`docs/S2_RULES.md`) apply to whatever s1
  the ledger keeps; an s1 that changes sends its s2 back to the captain.
- 2026-11-04 (e), process slip, no exposure: in the principle workflow (wf_880d37d6-91b), one recursive grep over its
  worktree (`*.md` and `*.json`) opened `boards/index.json`; no line matched and none was printed. The same grep
  printed three citation lines from `docs/research/` (landscape.md, external_data.md, t2_spatial.md), all giving the
  DOIs of errata to Qiu et al., Cell 2024; nothing from them was used or written. Recursive searches now skip both
  folders (`AGENTS.md`, quarantined material).
- 2026-11-04 (f), clarifying (d) (proposed by the rebuild design's synthesis, wf_237f8d8f-e7c, amended by its
  compliance audit; decided by the AI assistant under the captain's ruling (c) and the delegation of 2026-11-02; the
  captain can overrule). (1) A rule that changes which real cells are drawn inside a chamber category, by quotas or by
  weights, changes the shares of the states within it, so it is composition under (d): its rule, every constant
  included, is committed in its ledger row before any run, and before any diagnostic that reads a T1 truth's
  within-category shares as a target (the births entry's precedent). (2) Shares are set through `vec.core`'s quotas
  at the count the frozen count rule gives, never through importance weights: what a weighted draw realises depends on
  that count's ratio to the pool, and cell numbers stay closed. (3) T1 depth stays as built (every drawn cell keeps its
  own counts) and is now closed, as the depth line already is on the T2 boards. A cell count read with the T1 E12.5
  file size gives bytes per cell, which bears on detection depth, so (c)'s reason for leaving that size in the notes
  ("Task 1 scores no size") does not cover depth. The size is now redacted from `docs/CHALLENGE.md` at HEAD as
  `[redacted 2026-11-04: A]`; a check that printed locations only found no other copy in a tracked file outside
  `boards/` and `docs/research/` (those two were not searched). Earlier commits still hold it, and rewriting history
  is the captain's call. A share rule's row reports the drawn population's median library size against E10.5's, and
  a shift beyond 5% sends the rule to the captain. (4) Why shares stay open: a count, the upper-limit rule or a file
  size says how many cells or bytes a file holds, not which states appear or in what shares, as (c) found for births.
  Unchanged: cell numbers, sizes and scope on every board, T3's composition and the frozen s2 rules.
- 2026-11-04 (g), process slips in the rebuild design workflow (wf_237f8d8f-e7c), found by its compliance audit; no
  test-board limit, anchor or truth count and no other hidden value was seen, and the one hidden-target figure shown is
  the file size below: (i) four agents' first reads (the T1-populations, embryo, T3 and heart designers) displayed the
  `docs/CHALLENGE.md` line holding the T1 E12.5 file size (in the notes since 2026-10-18; not quoted, written or used;
  redacted at HEAD under (f)); (ii) the heart critic opened once a harness-saved tool output in Claude Code's session
  store, whose content was `TRIES.md` as its own command had printed it (no session log or transcript was opened);
  (iii) the heart designer started a filesystem-wide `find` for an installed package, which walked `boards/`,
  `docs/research/` and the session store by name only, and stopped it unread; (iv) some agents ran the repository's
  Python with `-I` but without `-B` (timestamps show no bytecode written in the tree); (v) design and critique texts
  used leaderboard figures beyond calibration: decision thresholds pinned to the field's validation totals and
  test-board statistics as evidence for or against a gap hypothesis (the T1-populations design and its critic, one of
  whose figures is in no committed record), the test top 3 turned into a per-group skill to say where the embryo's gap
  lies (the embryo design; its critic used a test-board statistic the same way), a sentence attributing the field's
  heart lead to closed levers (the heart design), and a proposed per-metric leaderboard read outside the screened
  route (the embryo design). All are struck: none enters a ledger row, a decision or a prompt, and the design texts are
  not committed. (vi) An embryo candidate proposed an external atlas spanning E8.5, a held-out heart validation stage
  strictly inside a heart interpolation bracket, for a check of a step in the shared core, before any window check;
  nothing was fetched or run, and it is pruned (`experiments/harness/TRIES.md`). From now on agents never open session
  logs, transcripts or harness-saved tool outputs (under `~/.claude-account` or `~/.claude`): a long output is split,
  or written to the scratchpad and read there. Packages are located through the venv's `*.dist-info` (`RECORD`,
  `direct_url.json`), never by a filesystem-wide search. The repository's Python runs as `python -I -B` (on the DGX
  `uv run --frozen --no-sync python -I -B`), and read-only git uses `GIT_OPTIONAL_LOCKS=0` (`AGENTS.md`).
- 2026-11-04 (h), process notes from the first-wave workflow (wf_b4b013b3-a46), from its agents' own reports and its
  review; no test-board limit, anchor or truth count and no hidden-target size was seen: (i) one oversized `sed` of
  `experiments/harness/TRIES.md` by the run agent was saved by the harness into the session store (its content was
  that file's text); it was not opened, and the rest was read in smaller chunks; (ii) the K2 builder read once the
  harness's output file of its own background wait command (empty; under `/tmp`, outside the session store, but a
  harness-saved output all the same); (iii) the K2 builder's capped read of `baselines_v1`'s T3 decision text showed
  the baseline scores that row K2a already quotes, before any K2 score existed; (iv) the N2 builder displayed
  `experiments/t1_m1/results/diag_v1/README.md`, which holds that diagnostic's practice-task scores and cell counts
  by Feng group (released and validation files, no test board); N2 uses none of it and stays target-free; (v) the
  builders left three repository copies on the DGX, each with its own `boards/index.json`, which a later recursive
  search could walk; after the K2 smoke evidence was archived (`experiments/harness/results/t3_k2_v1/factor_check.txt`;
  the full logs in the ignored `runs/smoke/`), the copies, the local worktrees and their branches were removed the
  same day. (g)'s rules stand: a long output is capped, or written to the scratchpad and read there. In addition, a
  workflow's last agent removes its scratch copies and worktrees once their evidence is archived, and the DGX sync
  excludes `.claude` (`docs/COMPUTE.md`).
- 2026-11-04 (i), the second screened leaderboard read (wf_010b5501-f2a, 14:40 UTC; the captain's choice of a screening
  agent for the leaderboards). The screening agent read the public leaderboard pages. Its audit was not clean: it
  flagged an item in the site's JavaScript and the boards' label strings, and withheld both before handing over. The
  main agent received a redacted table: overall and board totals, ranks and team names, nothing else. No withheld
  item is in any file, prompt or summary, and the main agent saw no test-board limit, anchor, truth count or
  hidden-target size. The figures serve as calibration only (`docs/PLAN.md` decision log 2026-11-04 (e)).
- 2026-11-05, the captain's upload gate: "Until our estimated score in leaderboard crosses 209, don't submit our 1st
  submission." Its figure is the leading Human total on the test boards in the screened read of 2026-11-04 (209.3).
  It is the captain's own upload rule and gates official uploads only. It is not a threshold or evidence for any try,
  adoption or method choice, which a leaderboard figure never is. The estimate it is compared with uses our own
  committed evidence only, never a test score or a leaderboard figure, and later leaderboard reads do not move the
  gate (`docs/PLAN.md` section 8). His deadline exception of the same day lets the first round go in from 2026-11-27
  00:00 UTC whatever the estimate, each upload only on his go-ahead at the time; it changes no data rule.
- 2026-11-05 (b), process notes from the DGX campaign design (wf_7d1d7fe1-e6e), from its compliance audit and its
  agents' own reports; no test-board limit, anchor or truth count and no hidden-target size was seen, and no
  leaderboard page, board file or board endpoint was opened: (i) before any row existed, the T1 composition designer
  applied Feng's E9.5 -> E14.5 composition rate one day forward from E9.5 and compared it with E10.5's committed
  non-cardiomyocyte split, the target of both counted T1 tasks; it was disclosed and marked not to be cited, its
  figures are carried nowhere, it counts as T1's seventh try (rule 4), and E10.5 is never cited for that rule; (ii)
  that designer's probes read Feng and JAX tables in the DGX repository's `data/external` directly with read-only
  stdlib one-liners, without a scratch copy (the modification times of `data/`, `data/raw` and `data/external` are
  unchanged); (iii) the T3 designer's probes ran the DGX venv's Python directly, read built test files in the DGX
  repository's `runs/` (outside its scratch copy) and computed per-cell gene-count shares and a median library size on
  them, ran the T3 scorer on synthetic truths (not a try, on the synthetic-toy precedent), and wrote a probe with
  hand-typed gene lists; no hand-typed list or pathway set enters a model without an `EXTERNAL_DATA.md` row,
  `docs/PLAN.md` 4.4 provenance and A15; (iv) the T1 expression designer made one general web search, whose snippet
  (a published late-fetal trend, outside the banned range) was not used; (v) the T2 designer fitted the committed
  mixtures to each counted embryo target's best mixing fraction (so E7.25 and E7.5 cannot count for a change to the
  mixing fraction; not quoted anywhere), read the digest's lines on how cell limits are set, and displayed a
  validation board's published spec; (vi) two full test runs on the DGX that the audit could not attribute were the
  main session's harness initialisations at 01337d0 and 54fccb0 (finished 19:06 and 19:19 UTC on 2026-11-04; the
  audit read the DGX's local times, UTC+1, as UTC), and a uv lock four minutes later is unattributed (a probe's
  `uv run` fits); (vii) the AI assistant's task prompt gave every design agent the per-board figures of the screened
  leaderboard read, the misuse (g)(v) struck, and one design used them to say where the gap lies. It is struck: board
  priority comes only from the headroom tables, and prompts carry no per-board leaderboard figure (at most the overall
  figure, labelled as calibration); (viii) residues that predate the campaign, kept until the captain rules: the JAX
  source tables outside the repository (`~/v/vec-gapconc/ext/jax/meta/`), which still hold boundary embryos 50-53
  (agents read JAX only under `data/external/jax_qiu2024`); raw copies of challenge-site pages, leaderboard downloads
  and board-endpoint responses from 2026-10-18 to 2026-10-20 in the AI assistant's session scratchpad (no agent lists
  or searches the scratchpad recursively, and probes write into their own subfolder); and the DGX scratch folder
  `rebuild_2026-11-04/` (the rebuild design's probe scripts). The campaign left no scratch copy on the DGX.
- 2026-11-05 (c), a disclosure gap the campaign's audit found, traced by a screening agent (the main agent read only
  the screened answer and the digest): the T1 s2's placement (cardiomyocytes at 60% of heart cells) leaned on three
  E14.5 cardiomyocyte shares from the research of 2026-10-19 (digest, 2026-10-19, the composition-cost item): one
  whole-heart dataset (redacted, (e)) and Feng's chamber dataset in two strains (41.5% and 76.6%). All three came from
  label-only per-cell metadata, before any register row covered them: X002 covered Feng's CD1 strain only, and the
  whole-heart source (best match Feulner 2026 with Wünnemann 2020, CELLxGENE-coded Theiler stage 22; the digest does
  not name it) had no row. Under this policy's "Nearness counts" and the staging rule's working answer
  (`docs/ORGANIZER_QUESTIONS.md` Q7: an embryo either label puts inside a window stays out), a sample coded Theiler
  stage 22 (E13.5-E15) stays out of T1 and the T2 heart, even though its authors call it E14.5. So: X002 now discloses
  the C57BL/6 label use (E14.5, legal); X006 discloses the whole-heart source and withdraws its figure from the s2
  rationale; X005 lists the published code the method runs on. No rule changes: Feng's two strains alone span 42-77%,
  and the frozen s2's 60% lies inside that range (`docs/S2_RULES.md` does not cite the withdrawn figure; the s2_hedge
  README and `docs/PLAN.md` 4.1 are restated). The figure stays in the screened digest as a record and is not cited
  again (superseded by (e): its value is now redacted from the digest and every working file). No history purge is
  proposed: it is an authors-labelled E14.5 summary share, outside the organisers' window as written, kept out by our
  stricter working answer; the captain can rule otherwise.
- 2026-11-05 (d), process notes from the gate-estimator workflow (wf_ae139ac6-f35: a builder, a blind independent
  recompute and a code audit; local only, aggregation of committed results, no target read and no score computed): (i)
  reading the assigned decision-log entries and `docs/PLAN.md` section 8 showed the builder and the reviewer the
  screened field figures and the gate's source figure, and a grep for the outlook figures printed `progress.md` and
  `docs/PLAN.md` lines that carry them; none was used, written or quoted, and the estimator holds no leaderboard
  figure (a test pins every number in it); (ii) key dumps printed practice-task fields of committed results (the T1
  grid's re-sampled E10.5 truth sizes and prediction counts, prediction counts in `gata4_summary.csv`, and the
  board-spec counts of the released-data K2 pseudo-task); no test board is among them, and the estimator uses none of
  these fields; (iii) a capped read of `docs/CHALLENGE.md` showed its validation-anchor line (not test anchors); (iv)
  the builder opened one memory note of the AI assistant (the compliance rules; no values), which is not a session
  log, transcript or saved tool output, but sits under `~/.claude-account`: from now on agents leave memory files
  alone too; (v) probe folders were deleted, and the tests delete their temporary copies.
- 2026-11-05 (e), the withdrawn whole-heart figure ((c) above, X006), after a side agent's note to the captain: its
  value is now redacted from every working file a method choice could read: `EXTERNAL_DATA.md` X006, (c) above,
  `docs/PLAN.md` 4.1, the s2_hedge README and the screened digest (the composition-cost item). It had stood in the
  digest and `docs/PLAN.md` 4.1 since 97d0c47 (2026-10-19) and in the s2_hedge README, as the low end of a range, since
  the s2 rules (38faf73); the disclosure commit (1435acf) added it to X006 and (c). Why: the staging rule keeps the
  sample out of T1 and the T2 heart, yet the value still faced every reader of a T1 or heart choice, among them any
  prior for D3 (the cardiomyocyte share), which is open. Two lists that named the dataset as legal without that
  exception (`docs/PLAN.md` 4.1's legal external data, `ASSUMPTIONS.md` A9) now carry it. What stays: the source and its
  withdrawal stay disclosed (X006; Rules s.10 asks for sources, not values); git history keeps the value as the record
  of what informed the s2; the unscreened notes under `docs/research/` stay quarantined as before (only screening agents
  read them). No history purge, the AI assistant's recommendation under the delegation of 2026-11-02: the figure is an
  authors-labelled E14.5 share from a published paper, outside the organisers' window as written, and a purge would
  rewrite every commit since 97d0c47, changing each commit hash cited as pre-registration evidence and, as in the purges
  of 2026-10-18 and 2026-10-19, recreating the GitHub repo and losing its push record. The captain can order one (a
  force-push needs his yes; the script would be dry-run on a scratch clone first, and he runs it). The side agent's note
  repeated the value to the captain, so a D3 prior, if he writes one, rests on registered legal data (Feng's two
  strains, X002), not on it. The captain ruled on 2026-11-05: no purge (`docs/PLAN.md` decision log 2026-11-05 (c)).
- 2026-11-05 (f), the cleanup the captain ruled on 2026-11-05 (`docs/PLAN.md` decision log 2026-11-05 (c)): (i) DGX:
  `~/v/vec-gapconc/ext/jax/meta/df_cell.legalT1.csv` and `lpm_legalT1.csv`, which also held boundary embryos 50-53, are
  deleted; the reduction's output in `data/external/jax_qiu2024/` holds embryos 1-22 and 54-74 only, and a rerun needs a
  fresh download (X003). `jax_labels.csv` (embryos 1-22 and 54-74) and the gene table stay, and the per-embryo folders
  hold no file for embryos 23-53 (checked by name). (ii) DGX scratch: `~/v/scratch/rebuild_2026-11-04/` and, under the
  standing rule for scratch copies, `~/v/scratch/births/` were archived to
  `runs/archive/scratch_rebuild_2026-11-04.tar.gz` and `runs/archive/scratch_births.tar.gz` (every file's hash checked
  after extraction; sums in `runs/archive/SHA256SUMS.scratch_2026-11-05`) and removed; neither held a repository copy or
  a board file. (iii) The AI assistant's session scratchpad: a screening agent classified its top-level entries of
  2026-10-18 to 2026-10-20 by name and file type, relaying no values. 29 hold raw challenge-site copies (site pages,
  Nuxt and JS bundles, API and leaderboard captures, a Data-page copy, a copy of the board-anchor index, inventories of
  challenge files with their sizes), 3 could not be opened past a boards/, panels/ or leaderboard subfolder, and 130 are
  our own work. The auto-mode classifier blocked the AI assistant's deletion, so the captain ran the script it prepared
  (`cleanup_2026-11-05/delete_site_copies.sh` in that scratchpad: a dry run unless given --delete, with a date guard) on
  2026-11-05: 32 entries and 514,772,564 bytes deleted, nothing left behind. The screen also flagged the same kind of
  material from later dates (manifest tables with per-file sizes of 2026-10-27 and 2026-10-28, site copies of
  2026-10-28, the leaderboard reads and captures of 2026-11-04) and, among our own work, anchor-check and parity-hedge
  outputs, a raw Slack copy, a repository bundle of 2026-10-19 that may keep pre-purge history and the raw
  research-agent reports; a second script (`delete_later_copies.sh`, 17 entries) was recommended, and the captain ran it
  the same day: 32,163,793 bytes deleted, nothing left behind. (iv) Slips: the screening agent once read the first 1.2
  KB of an inventory file through `grep -o` to see its field names (it relayed no values); the AI assistant wrote two
  scripts of this task into an existing scratchpad folder (`r1/`) that the screen found to hold site downloads of
  2026-11-04, without listing or opening its other contents.
- 2026-11-05 (g), process notes from the search of decision log 2026-11-05 (d): (i) the three reviewers
  (wf_6b41d823-d5a) report: one parsed committed CSV and JSON files with the system `python3 -I -B` rather than `uv run`
  (no repository code ran that way); one opened the harness decision files of the T3 pseudo-tasks (practice scores only)
  and read section 2.5 of the screened digest, which it says holds a broad-phenotype paragraph on the held-out knockouts
  (not used); one grep over the decision log printed two lines that quote screened per-board leaderboard figures (not
  used or repeated). The last two go to a screening agent: per-board leaderboard figures do not belong in files every
  agent reads, and the digest paragraph is checked for gene-level knockout findings. (ii) The AI assistant cut the
  upload client off with `| head -2` during two format checks, which then never reached `results/submissions.csv`; both
  were rerun and logged, and the client showed no unresolved upload. (iii) A side agent's note told the captain that the
  T1 blood hedge's top-end gain rests on the grid's design prior; it was right, and the hedge was dropped.
- 2026-11-05 (h), a screening agent's cleanup of the files every agent reads, after (g): six per-board leaderboard items
  were removed and overall totals kept. Five are in `docs/PLAN.md`'s decision log (per-board score ranges of the Human
  top 3 and of the teams in the test top 10, a per-board validation-to-test rank correlation, the 6th Human place per
  board, and a claim of which board our gap to the field is largest on); the sixth is that same claim in this file's
  2026-11-04 (d) rationale. Line counts are kept, since `experiments/gate/estimate.py` cites `docs/PLAN.md` by line. In
  the screened digest's section 2.5 it removed a count of truth-DE genes on a T3 board whose source cannot be determined
  (it is not the validation board's, and the digest names no board); an exact-value search found the number nowhere else
  in the repository [both statements corrected in (i)]. Under the AI assistant's brief, which was stricter than this
  policy (the ruling of 2026-10-15 allows phenotype-level knockout knowledge), it had also replaced the digest's
  paragraph on the held-out knockouts' phenotypes; that edit was undone, so the paragraph stays, like the same kind of
  statement in `docs/PLAN.md` 4.4 and `METHOD_SUMMARY.md`. Git history keeps the removed text; no figure was relayed to
  the main agent.
- 2026-11-05 (i), corrections to (h) after a code review of 030ca89 (13 findings, 4 confirmed) and a second screening
  pass. (i) The count: the second pass searched every tracked file, `docs/research/` included (`boards/`, `panels/`,
  `data/`, `runs/` and index.json files excluded), and found the same quantity once more, in
  `docs/research/t3_methods.md`, now marked; the first search had skipped `docs/research/`. No T3 rule, script or result
  holds the value or anything derived from it (n_up and n_dn are counted from the truth at run time). Its likely source
  is an example in the organisers' scorer docstring, written before the validation answers were released and naming no
  board; the validation board's own truth-DE count under the scorer protocol varies by seed over a range that includes
  it, so it may describe either board, and the digest's mark now says "possible". From 97d0c47 (2026-10-19) it stood in
  the digest's section 2.5, which the AI assistant read when it built PLAN v1.2 and the T3 rules of 2026-10-18 and
  2026-10-19 were written; a reviewer of 2026-11-05 (wf_6b41d823-d5a) read that section before the cleanup. Nothing in
  the repository shows it was used. (ii) History: older revisions of `docs/PLAN.md`, `docs/DATA_POLICY.md`, the digest
  and `docs/research/t3_methods.md` hold the removed text, pushed to GitHub and in the DGX's .git; AGENTS.md forbids
  printing them. Rewriting history is the captain's call (`docs/PLAN.md` decision log 2026-11-05 (e)); the AI assistant
  recommends no rewrite, for the reasons of (e), and that any public release of the code be a fresh repository built
  from the final tree. (iii) The digest: the second pass removed 18 more per-board leaderboard items (leader scores and
  ranges per board, per-task prize lines, per-metric point pools at leader levels, per-board gaps) and its legend now
  keeps overall totals and ranks only. Other screened copies under `docs/research/` (t2_spatial.md, t3_methods.md) still
  hold such figures, behind the rule that the main agent reads only the digest. The edit (h) undid was the AI
  assistant's own `git checkout` of the whole digest, made before the count was removed and with no diff printed; a
  screening read in the review found section 2.5 phenotype-level only for the held-out knockouts. (iv) Fixed after the
  review: the `docs/PLAN.md` citations in `experiments/gate/` are entry labels (they had drifted); the one placeholder
  that named a board no longer does; the Scores section holds the standing screening rules. Each "[per-board leaderboard
  figure removed 2026-11-05]" mark in `docs/PLAN.md` and this file is explained by (h) and (i). (v) Slips: the review's
  agents had one oversized grep output saved in the session store by the harness (not opened) and first wrote working
  files into an existing scratchpad folder (names listed, none opened); the second pass ran the system `python3` once on
  its own scratch script, and, as its brief asked, searched `docs/research/`, which showed hit lines from unscreened
  reports to it only. The captain ruled the same day: no rewrite, and any public release of the code is a fresh
  repository built from the final tree (`docs/PLAN.md` decision log 2026-11-05 (e)).
- 2026-11-05 (j), a record update, no exposure.
  - A screened sweep (wf_cf28f5b9-9aa) read Slack, the organiser email and the Rules, Timeline, FAQ and Prizes pages,
    each digest audited (`progress.md`, Blockers). It found nothing new from the organisers.
  - The pages showed text this record had not summarised, among it s.10's midpoint nearness rule. A screening agent
    then read s.10's window sentences word for word. They are now under External data windows above: for the
    extrapolation targets the explicit window relaxes the midpoint rule, so legal post-E13.5 data stays legal.
  - The other new text is listed in `docs/CHALLENGE.md` (Sources).
  - The digests skipped five Yukon passages on the site pages, recording only place and date.
