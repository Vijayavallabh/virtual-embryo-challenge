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
