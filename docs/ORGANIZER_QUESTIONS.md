# Questions for the organisers

Send privately (email virtual.embryo.moonshot@gmail.com, or a Slack DM to Weize Xu), not in the public channel: some
questions reveal method choices. Rules s.10: "an answer given in advance is not a finding of misconduct afterwards".
Record each answer with its date in `docs/DATA_POLICY.md`.

Status: **sent** by email to virtual.embryo.moonshot@gmail.com from be23b041@smail.iitm.ac.in on 2026-10-18. Answers:
_pending_ (record them in `docs/DATA_POLICY.md`). The leak report went separately by Slack DM the same day. As of
2026-10-20 22:43 IST neither email nor the leak report has a reply (the sender's mail tracker reported both emails
unopened 24 hours after sending), and no organiser has posted on Slack since 2026-10-16. Question 6 (the deadline's
time zone) was answered for everyone by the organisers' email of 2026-10-20 22:03 IST: "through 2 December 2026
(UTC)", no hour. A third email was saved as a draft in the sender's mailbox on 2026-10-20 for the user to review and
send. After the 22:03 IST email and the new partner rules page it was cut to two questions (the test phase's start
time; validation answers as training data and model inputs) plus a reminder of the two open emails: the Timeline
page already says the validation answers come with the test phase, Rules s.4 and the partner page settle what
ordinary teams do about Yukon's work, the deadline is UTC, and whether the validation boards stay open this week
changes nothing we do (the upload path was tested on them on 2026-10-20). The Generality question (question 4 of the
18 October email, asked again) was added the same day (prize strategy, `docs/PLAN.md` decision log 2026-10-20 (c)),
and the user sent the email with three questions on 2026-10-21 at 01:24 IST (2026-10-20 19:54 UTC); its text is
under "Third email" below. The user sent a fourth email (the HTTP 403 on the panel files, final entries, confirmation
of the A20 reading) on 2026-10-28 at 15:41 IST; its text is under "Fourth email" below. As of then none of the four
emails had a reply, and the last screened Slack sweep (2026-10-23) found no organiser post or reply.

## Screened reads of 2026-11-19: two organiser replies of 16 November, found inside Slack threads

Two reads, each by a screening agent whose released block a second agent audited before the main session read it
(workflows wf_56e44340-f54 and wf_a539cc31-54a); the first read was made on 2026-11-19 at about 04:38 UTC (its own statement) and the second between 04:41 and 04:44 UTC
(the main session's clock readings before and after it).
The first is the read the captain approved for 19 November (`docs/SCREENED_READ_SEVENTH_EMAIL.md`), with one item its
reminder added (whether the public Slack question of 2026-11-12 has an organiser answer). The second repeats the Slack
part of the check of 2026-11-18 with the channel's threads included, on the open questions that check listed; the AI
assistant started it on its own when the first read showed that the earlier check had read the channel's top-level
messages only. Recorded by the AI assistant.

- **Seventh email (processing batches of the two held-out samples): no answer found** in the organiser mail thread
  (no organiser reply), in the Slack channel or on the Data Protocol page (HTTP 200; whether its text changed since
  2026-11-10 cannot be told without an earlier copy). By the reading fixed beforehand, nothing changes.
- **The public Slack question on a floored Task 3 severity: answered** by an organiser (Weize Xu) on 2026-11-16 at
  19:53 UTC, in the question's thread, as released through the screen: "Yes, the current P3 final-test scorer uses the
  same behavior. No-response predictions, or cases with beta ≤ 0 or R² ≤ 0.01, receive the finite severity_slope
  fallback of log(1e-3), matching wt_identity and mapping to skill 0.5 for that component. There is no additional
  penalty specifically for triggering this fallback. Other T3 components are scored independently, so this does not
  imply an overall T3 score of [n]." (The first read's auditor had also replaced the two gate values and the phrase
  on wt_identity; the second read's auditor let them stand as scoring conventions, which they are in
  `docs/CHALLENGE.md`.) By the reading fixed on 2026-11-11 ("Public Slack question" below), a "yes" keeps T3's last
  upload as it is; the switch to the barred s1 that an answer of zero would have made does not apply. No follow-up is
  drafted.
- **External datasets that span a protected window, and where sources are disclosed: an organiser's reply to another
  participant** (Weize Xu, 2026-11-16, 19:52 UTC; the participant had asked whether a public dataset spanning the
  protected window may be used once every cell inside the window is removed, and where to disclose external
  sources), as released: "Yes, provided all cells within the protected windows are removed before training, and the
  retained data do not carry information derived from the excluded cells. Using raw expression counts and stage
  labels from the permitted cells is acceptable. Please exclude annotations, clusters, embeddings, or fitted
  transformations learned using protected cells, and recompute any such processing on permitted data only. Describe
  the filtering in your method report. On the submission page, use the P3 method report → Data used field to list
  external datasets, links/accessions, stages used, and filtering details. Disclose any pre-trained models or
  external code in their corresponding fields, and attach the report to your P3 submission."
  - **What it touches** (the AI assistant's reading; a rule question for the captain, `AGENTS.md`, "Escalation"):
    register row X002. T1's direction is computed per cell group from Feng's permitted cells (E9.5 and E14.5), and
    the groups are the authors' own cell-type labels (`cell_type_spec` in the authors' metadata table;
    `experiments/external/feng2022_geo.py`, `experiments/t1_m1/feng_dirs.py`). If the authors made those labels in
    an analysis of all their stages, the in-window ones among them, they are "annotations ... learned using protected
    cells" in the organiser's words. Both T1 files rest on them (the uploaded first file and the built second), and
    the second file's cardiomyocyte share of 0.60 was placed inside a range read from the same labels (X002). No
    in-window cell, count or expression value entered anything: the removal the organiser asks for was done. What is
    in question is the labels of permitted cells.
  - It also answers the open question of where sources are stated: in the submission page's method-report fields,
    with the report attached to the submission. The runbook had planned an email.
- **No organiser statement since 2026-11-13** on: the validation files and the Gata4 replicates as model inputs; the
  scorer build for occupancy_dice or a recomputation of shown scores; how final entries are nominated; the use of
  returned scores, deadlines, the number of uploads or cell-count limits. No announcement to all teams since
  2026-11-10. The third organiser reply of 16 November says only that questions belong in Slack.
- Coverage, as the screen states it: the channel's four top-level messages since 2026-11-10 and all three threads
  with replies; the only organiser posting in that window is Weize Xu, three thread replies on 2026-11-16. Withheld by
  the screens, by class only: one passage on test-board limits, two on scores or leaderboard figures, three of other
  participants' text, three Data Protocol passages not of the asked kind, three mailbox items that are not
  organiser replies; one figure inside the severity quote reads [n].
- **A correction to the check of 2026-11-18 (next section):** its line "The Slack channel has no organiser post
  since 2026-11-13" holds for top-level posts only. Three organiser replies of 2026-11-16 sat inside threads, which
  that check did not open; its "no answer found" on how a floored severity scores was therefore wrong on the day it
  was written. A screened Slack check reads every thread with replies from now on.

## Disclosure of the X002 labels (drafted 2026-11-19, rewritten the same day; in the captain's drafts since 2026-11-19; the captain sends)

Put into the drafts of be23b041@smail.iitm.ac.in, to virtual.embryo.moonshot@gmail.com, on 2026-11-19 at 11:24 UTC
(the draft's own time stamp), at the captain's word of that day ("Yes, both into drafts (Recommended)", `docs/PLAN.md`
decision log 2026-11-19 (f)), as the text below for the no-outside-data case, and read back from the drafts folder
against it. Whether and when it was sent, and any edit before sending, are recorded here when the captain says.

First drafted and reviewed in wf_63571f3d-9a2 (the text stood here until the rewrite; it promised recomputed groups
and left a bracket on the 60% share). Rewritten in wf_038f7126-b04 for the case the AI assistant now advises, a
last T1 upload with no outside data (`docs/PLAN.md` decision log 2026-11-19 (d) and (e)): a writer, three reviewers
(accuracy, strategy, tone: eight points, applied by a fourth agent). The captain has ruled on none of it. Send
conditions, unchanged: (1) it goes only after the captain has ruled which file is T1's remaining upload, with its
paragraph on that upload matching his ruling; (2) the sources note of 18 November (below) goes in the same sitting
in its amended form, and the copy in his drafts, which still lists `T1_s2.h5ad`, is replaced or deleted first. Order
advised by the reviews: this email first, the sources note within minutes, both as new emails. 295 words from the
greeting to the signature. Answers: _none; not sent_.

**Subject:** Human Team (Mithril): disclosure on the cell-type labels in our Task 1 entry, and one question

To: virtual.embryo.moonshot@gmail.com

Hello organisers,

One disclosure about our Task 1 entry, and one question; please keep it private.

Our first Task 1 upload used one external expression dataset: Feng et al., Nat Commun 13:7960 (2022), GEO GSE193346,
raw counts of wild-type CD1 hearts at E9.5 and from E14.5 on. Cells of E10.5 to E13.5 were removed on the per-cell
stage column before any model used the data. One step of our method used the permitted cells: it took each cell type's
direction of expression change from E9.5 to E14.5.

Those cell types were the authors' labels. After your public reply of 16 November, we checked the paper's Methods.
They say the authors merged every sample of every stage into one object, batch-corrected and clustered it, so the
labels were learned using protected cells.

Our first Task 1 upload used these labels: T1:test, 18 November 2026, 14:37 UTC, file T1_s1.h5ad, SHA-256
a294062f5f7bbb2dae19ef695979cc21f03f38df40d22be9ccdd30af36e3bdef. Apart from these labels, no in-window cell, count or
expression value entered any file.

We also built a second Task 1 file. Its cardiomyocyte share was placed within a range read from the authors' labels.
That file has not been uploaded.

Our remaining Task 1 upload uses no external dataset: it is our first file's method with that one step taken out. Two
parts it keeps were made from released data only but chosen with evidence that used the authors' labels.

Does our first Task 1 upload stand as an eligible entry, or should only our remaining one count?

As your public reply asks, the P3 method report's Data used field will list this dataset, its accession, the stages
used and the filtering, and the attached report will describe this history and every choice resting on the authors'
processing.

Best regards,
Vijayavallabh J
Team captain, Mithril

**If the captain rules for the recomputed-labels file instead** (Feng's direction kept, groups recomputed from the
permitted cells), two sentences change, by the same reviews:

If the captain rules for the recomputed-labels file (Feng's direction kept, groups recomputed from the permitted
cells), make two changes to email A. The email is then 298 words (counted by wc), so it is still under 300.

1. Replace "Our remaining Task 1 upload uses no external dataset: it is our first file's method with that one step
taken out." WITH: "Our remaining Task 1 upload keeps that step, with groups recomputed from the permitted cells' raw
counts using markers from released data and published statements."

2. Replace "Two parts it keeps were made from released data only but chosen with evidence that used the authors'
labels." WITH: "Two other parts were made from released data only but chosen with evidence that used the authors'
labels." (That file keeps g_mom and the cell base too.)

Everything else stays, the share sentence included. The recomputed-labels file is the first file's design and carries
no 0.60 (decision log 2026-11-19 (d), the cold challenge).

The words "raw counts and stage labels only" are not used. The groups also rest on marker genes, from the released
E9.5 stage and from published statements on adult mouse tissue in eight papers (decision log 2026-11-19 (d); proposed
register row X010, corrections still owed). Those papers must be listed in the method report.

Do not add "by a procedure fixed before it runs" today. The labeller's final review says not ready, with three
blockers, and its settings are a draft (decision log 2026-11-19 (d)). Once the settings are signed off, the clause ",
by a procedure we fixed and recorded before it ran" may be added. To make room, cut "batch-corrected and" and "only"
from "released data only", or the email goes over 300 words. [TO COMPLETE AT THE RULING]

**What the email leaves to the method report** (the reviewers' note): which candidates and priors the share range
bounded; the authors' quality control and their zone labels, both made on the full dataset; the two kept parts by
name (g_mom's adoption evidence, and more weakly the pricing behind the cell base).

## Optional public Slack question on the first upload's eligibility (drafted 2026-11-19; not posted; the AI assistant advises against posting)

From the same workflow, as an option the cold challenge of the T1 plan raised: the organisers answered two public
Slack questions within three days in November, while the team's emails and notes have had one reply (17 November).
The AI assistant advises against posting it: it would put a detail of the team's P3 method into a public channel,
against the rule that P3 methods stay private until the announcement, which only the captain can set aside for one
question, as he did on 2026-11-11 for a scorer question; and other teams could read from it that a team has this
problem. If the captain wants it at all, it goes only after the email has had no answer for some days. Text:

> Hi organisers, a question on your 16 November reply about external datasets that span a protected window. Suppose a team made an official upload with an external dataset's cell-type labels, and later learned that the authors had made those labels with protected-window cells included. If the team has disclosed this to you and its remaining upload does not use those labels, does the first upload stay eligible? A yes or no is plenty. Thanks!

## Screened check of 2026-11-18 (after the first uploads): no answer to any open question

Run by a screening agent and audited by a second agent before the main session read it (wf_bfc2939d-aa0); the reads
were made on 2026-11-18 from about 15:40 to 15:46 UTC. The check approved for 19 November was brought forward by the
captain's word of 18 November ("go ahead with whatever is best for us and immediately as weell instead of waiting";
`docs/PLAN.md` decision log 2026-11-18 (k)). Sources: organiser messages in the captain's mailbox since 2026-11-13,
organiser posts in the public Slack channel since 2026-11-13 and the captain's direct messages with the organiser,
the Rules page, the Data Protocol page, and the leaderboard page for overall totals only.

- **No answer found** to: the validation files, the Gata4 replicates included, as model inputs (fourth email,
  question 3); how a floored severity scores (fifth email); the scorer build and parity (October); how final entries
  are nominated; where a method summary goes and how sources are stated; the processing batches of the two held-out
  samples (seventh email; the Data Protocol page makes no such statement); the further reply the organisers
  announced on 17 November. The Slack channel has no organiser post since 2026-11-13.
  (Corrected on 2026-11-19: this check read the channel's top-level messages and no thread. Three organiser replies
  of 2026-11-16 sat inside threads, one of them the answer on a floored severity; see "Screened reads of 2026-11-19"
  above.)
- **The Rules page on score use (Section 10), as released** (two quotes are shortened by the reading tool):
  "Returned scores may be used to choose among predictions, not to compute them." "Using score returns to solve for a
  property of a held-out target, by inverting the published skill map, ... is recovering held-out data by another
  route and is prohibited." "The test is whether the truth enters your prediction's content or only its selection."
  "Trying a parameter value, seeing how it scores, and trying another is ordinary leaderboard work and remains
  permitted however many rounds it takes; ... what is prohibited is solving for a property of the target and
  building a submission out of the solution."
- **Section 9, as released:** "Choosing which completed run to submit is a person's decision and always has been."
  "A final entry nominated for prizes must carry its evidence." No nomination procedure was found on the page.
- **The leaderboard page** showed no figures (its table did not load), so no overall total or rank was read.
- Withheld by the screen: two passages (one on test-board limits and hidden-target properties, quoted from the
  team's own earlier text; one Data Protocol sentence on how the E12.5 samples were collected, marked unsure).
- **What the team does:** by the rules fixed beforehand nothing switches: T3's remaining upload stays a file that
  takes no number from the Gata4 replicates, and the working answers stand as working answers.

## Organiser messages received (2026-11-17): one answer, one notice

Both were pasted into the main session by the captain on 2026-11-17; `docs/PLAN.md` decision log 2026-11-17 (a) says
what each changes, and `docs/DATA_POLICY.md` 2026-11-17 (a) how they were handled. A real answer replaces a working
one the day it arrives.

1. **The answer on the exposure notice of 4 November** (an email to the captain; the first reply the organisers have
   given the team). In full, as pasted: "Thanks for flagging the metadata issue, and sorry for the delayed reply. We
   have removed the hidden-target file sizes and cell-limit derivation from the website. The scoring rules and
   submission limits remain unchanged. Thank you for documenting the exposure and explaining that your predictions did
   not depend on it. Please include that explanation and your precommitted procedure in your method report so we have
   a clear record. This reply addresses the metadata issue; the other questions in this thread still need a separate
   response."
   - It replaces the working answer in the table of 2026-11-10 ("Exposure notice of 4 November") and agrees with it.
   - What the team does: the block of 2026-11-04 stays; `METHOD_SUMMARY.md` carries the explanation and the
     precommitted procedure in a section of its own; the correction that the working answer set aside for the
     pre-upload email (T3's count default was fixed on 2026-11-01, after the release and before the exposure) goes
     into that section as well.
   - The separate response it announces to the thread's other questions has not come. When it does, a screening
     agent reads it first and it is read as the sections of those questions below fix.
2. **The deadline notice to every team** (through the challenge's message system; thread
   b6a676c673a047d2ac013675e532a302). As pasted: "We are moving the final submission deadline one week earlier, from
   December 2 to November 25, 2026, at 23:59 AoE. Final-entry nominations and Agent Team evidence are due at the same
   time. [...] Official evaluation will begin on November 26, and winners will be notified by December 1. Nothing else
   changes: the tasks, data, metrics, two official submissions per board, and eligibility rules remain the same.
   Final rankings will still be based on the hidden test set."
   - It replaces "through 2 December 2026 (UTC)" (the organisers' email of 2026-10-20). 23:59 AoE on 25 November is
     11:59 UTC and 17:29 IST on 26 November.
   - It bears on the working answer for final entries ("editable until 2 December UTC"): nominations are due at the
     new deadline. Whether an unmarked board counts its best upload is still not stated.
   - What the team does: every upload date and the upload gate's exception date are the captain's to move; the
     AI assistant's proposal and what stays unchanged until he rules are in the decision log.

## Working answers (2026-11-14): the four questions that had none

The captain, 2026-11-14, after being told that two questions for the organisers waited on a decision and that the
leading total is out of reach: "so far organizers have not answered any of my single mail too", and "Assume logical
replies from the organizers for relevant mails". Four questions had no working answer that anything rested on: the two
of the seventh email, the Task 3 question cut from it, and the drafted eighth email's. These are working assumptions,
never organiser rulings. A real answer replaces one the day it arrives, and where only an answer given in advance
counts (Rules s.10) an assumed reply gives no protection.

How they were made (`docs/PLAN.md` decision log 2026-11-14 (a)). Three independent predictors (the lab that made the
samples; an organiser answering a short email; a statistician on base rates) worked from one digest of released-file
facts and the screened Data Protocol text, and two reviewers attacked their answers, one for where they suit the team
too well and one for the opposite (workflow wf_48f6fcf6-157, read only). The figures are the five agents' range and
its middle; they are judgements.

| Question | Assumed reply (the agents' range) | What the team does |
|---|---|---|
| Final entries (eighth email, drafted, not sent) | Marked per board; the marked upload counts and, where nothing is marked, the best official one (0.60 to 0.80, middle 0.75). One pick over a set of boards: 0.05 to 0.13 | No change. Every board is marked explicitly after its last upload, which removes the dependence on what an unmarked board counts. The draft stays unsent unless the captain decides otherwise |
| Task 1: E12.5's libraries and E10.5's library `_3` (seventh email) | Prepared separately, with the same protocol (0.45 to 0.54); prepared together and sequenced in one run (0.20 to 0.40); another protocol (0.10 to 0.30) | No change. Either of the first two confirms the first upload as built (row R6, drawn from library 3). Request R-B is not reopened by an assumed reply: its fixed reading below names an organiser's answer |
| Embryo: E7.75's sections (seventh email) | No single reply: processed with E8.0 (0.30 to 0.55, middle 0.50) or in a run of its own (0.33 to 0.62, middle 0.42); with E7.5 0.08 to 0.12 | No change. Without an answer the fixed reading below keeps the pair as built, and the pair rule that chose it holds from 0.30 to 0.60. The first upload built for a target like E7.5 stays declined: its trigger names the reply the agents rate least likely |
| Task 3: beta-catenin and the Gata4 replicates (cut from the seventh email) | On a slide and in a run of its own (0.65 to 0.80, middle 0.73): each released knockout line has its own slide (M125, M164) and both littermates of one line share it | No change: both uploads are fixed by frozen rules. The question need not be sent |

What the assumption settles. Nothing in the plan waits for an organiser reply. The screened reads approved for 19 and
23 November stay as checks, since a real answer would still replace a working one. The eighth email and the Task 3
question need not be sent; sending either stays the captain's choice.

What it does not do. It moves no file, reopens no request and lifts no closed lever: the readings fixed below for each
possible answer are readings of real answers, and nothing changes without the captain's ruling. It does not move the
committed estimate either. The forecast's judgement weights are not defined by these replies;
`experiments/gate/forecast_weights.py` gives the forecast's sensitivity to them for the files as built: 194.4 as
committed, 182.5 to 203.7 with every weight at one end, so no set of replies reaches the leading total.

## Working answers (2026-11-10): every question sent, predicted, not ruled

The captain, 2026-11-10: "so far they have not replied to a single email. So assume their logical position with the
contrast of the leaderboard scores to judge", and "i meant for all the emails sent to them so far". These are working
assumptions, never organiser rulings. A real answer replaces one the day it arrives, and acting on a prediction gives no
protection under Rules s.10, where only an answer given in advance counts. Confidences are the AI assistant's, after an
independent prediction from the rules alone and two adversarial reviews, one for where the team is too strict and one
for where it is too loose (`docs/PLAN.md` decision log 2026-11-10 (x)). They replace the table of 31 October below where
they differ.

The overall standings bear on only two rows: the stakes of the Generality request, and weakly the board specs, since 91
ranked entries show that the specs reach teams. Their main use is scale, in the decision log.

| Question | Logical position (confidence) | What the team does |
|---|---|---|
| Released data across tasks (18 Oct Q1; A15) | Yes: released and unlocked official files may train any task; the windows bind external data (0.8) | No new cross-modal input. What it would unlock is measured as about nothing (heart concordance null; T3 programme +0.02 to +0.05). Feng on the heart is legal external data rather than an A15 matter, and its null holds anyway |
| Textbook background in a protected window (18 Oct Q2) | Yes if qualitative: it may choose populations, with no measured in-window values (0.75; their 15 October answer on broad background) | No change: births are estimated at about -0.2 (range -2 to +0.6), and today's T1 pair is within about 1.5 of the best under every stated prior |
| truth_cells for the count (18 Oct Q3) | Yes where published; the test boards publish only limits (0.7) | Moot: max_cells by the committed count rule |
| Generality scope (18 Oct Q4; 21 Oct Q3) | Must underlie every task entered; heads and training may differ; judged from the summary (0.85) | See the next row |
| Generality as corrected (sixth email) | No ruling in advance. As described, more likely not one architecture under every task: about 0.35 eligible if our request is reached | The request stays, at no cost; no further work and no method change for it. The standings make reaching it unlikely: every team above us has to not ask or not qualify |
| Scorer (18 Oct Q5; 19 Oct Q1) | The live board's metric set and build through the final evaluation; PR #6 not adopted; parity pairs not endorsed (0.8) | No change: no parity pair; the row-order seed fixed in code |
| Leak report (18 Oct; 19 Oct Q8) | Thanks; no penalty for reporting; leaked values stay barred; targets and scoring unchanged (0.8) | No change: never used |
| Wrong-sign severity (19 Oct Q2; fifth email; Slack DM; the public Slack question of 2026-11-11) | Skill 0.5: the scorer's comment calls the case a finite floor, "not credited and not dropped", and P3 scores are published on upload, so a later change would rescore published results. The published wt_identity 50 fixes only the no-response branch, and the "undefined" label of our questions is our own paraphrase (decision log 2026-11-11 (e)) (0.85) | Today's T3 s2 stands (decision log 2026-11-06 (b)). If they ever say zero, that ruling's switch applies |
| Gata4 replicates in T3 (19 Oct Q3; 28 Oct Q3) | Yes: the unlocked files may inform the beta-catenin prediction; external perturbation data of either gene stays barred (0.88) | Allowed branch for s1; s2 holds no Gata4 number. The pre-upload email says so |
| Duplicated and interpolated records (19 Oct Q4) | Pipeline artefacts; truth sets processed alike, not cleaned (0.6) | No change |
| z coordinate and test stacks (19 Oct Q5) | Points to the Data Protocol page; stored z as given; nothing more on the test stacks (0.6) | No change |
| truth_cells for tissue size (19 Oct Q6) | No: counts are format limits, not a route to a hidden size (0.75) | No change: sizes come from the pre-registered procedure |
| Stage conventions for external data (19 Oct Q7) | Standard somite and Theiler conversion; conflicting labels stay out (0.6) | No change: ambiguous embryos stay dropped |
| Validation answers as training data and inputs (21 Oct Q2; 28 Oct Q3; A20) | Yes, E10.5 as the starting stage included (0.92) | No change |
| Board specs, HTTP 403 (28 Oct Q1) | An unintended Cloudflare block; a browser copy of the linked file is fine (0.8) | Done: saved file, read by code |
| Final entries (28 Oct Q2) | Marks per board, editable until 2 December UTC; an unmarked board counts its best official upload (0.7) | Every board is marked explicitly after its last upload |
| Notice of 31 October | No objection to the browser save, validation inputs or per-board marks (0.8) | No change; silence is not an answer under s.10 |
| Exposure notice of 4 November | Seeing what every signed-in team sees is no breach and carries no penalty; keep it out of the methods (0.8) | The block of 2026-11-04 stays. One correction goes in the pre-upload email: T3's count default was fixed on 2026-11-01, after the release and before the exposure. **Answered on 2026-11-17** (above): the exposed metadata was removed from the site, the rules are unchanged, and the explanation and the precommitted procedure go into the method report |
| Embryo-series stage labels (Slack DM, optional) | Timed matings refined by morphological staging; probably unanswered (0.55) | No change |
| External sources "with each submission" (not asked; `EXTERNAL_DATA.md`) | The form has no disclosure field, so a dated private statement before the upload meets it (0.8) | Private email before the first upload (`docs/P3_RUNBOOK.md`); `--model sources-disclosed-by-email` on every upload |
| Where the method summary goes (not asked) | The organisers' address, since no route is published for a claimant outside the places (0.55) | Sent privately in the plan's 24-30 Nov slot unless the screened check finds a route |
| Human Team status (not asked) | A team steering AI tools, with a person deciding every upload, is a Human Team; the Agent track is locked, unsteered runs (about 0.97, low confidence: the Rules' definition is in no file the review could read) | The screened check before the first upload quotes the Rules' definition |
| Processing batches of two held-out samples (seventh email, 2026-11-13) | Probably unanswered. If answered: a one-line procedural reply, perhaps posted for every team (a useful reply before the freeze at about 0.05 to 0.1, the strategy reviewer's estimate of 2026-11-13; not part of the review of 2026-11-10) | R-A and R-B as they stand. An answer is read as "Seventh email" below fixes it, through a screening agent |

## Working answers (2026-10-31): predicted, not ruled

Notice sent: the captain replied on the fourth-email thread on 2026-10-31 at 18:59 UTC (drafted by the AI assistant):
unless told otherwise, the team saves the board index with a browser not before 2026-11-04 (UTC), uses the unlocked
validation files (the Gata4 replicates included) as inputs, and marks one final entry per board before 2 December.

Follow-up sent: the captain replied on the same thread on 2026-11-04 (drafted by the AI assistant). The board index
still answered 403 to our tool; instead of saving it unopened, the captain displayed it in a browser and pasted it,
and then the Data page, into the AI assistant's session; its test boards carry no anchors or truth_cells; read with
the Data page's rule for the upper cell limit, some test-board limits give away the size of a hidden truth half, and
the page lists hidden targets' file sizes; nothing in our predictions depends on what was seen (the cell counts and
sizes follow rules committed before the release); the file was then saved for our intake. The same day the captain
checked the Gmail Spam folder, which the screening connector cannot see: nothing from the organisers since 2026-10-28.

No organiser had replied to any of the four emails or to the private leak note by 2026-10-31, and on Slack they have
answered no team's rule question since 2026-10-16 (screened sweep, 2026-10-31). The captain then told the AI assistant:
"Now assume that they don't reply to the emails and you then come up with their rational responses and act on it and proceed". Seven AI agents (a Slack screen; three independent predictors with the lenses "challenge
operator", "rules text" and "precedent"; a judge; two red-team adversaries; workflow wf_65dab548-d7d) predicted the
organisers' most likely answer to every open question, and a screening agent read the organisers' Data Protocol page
(posted 2026-10-29). These are the team's working assumptions, never organiser rulings: a real answer replaces a
prediction the day it arrives, and every disclosure says which is which. Confidence is the judge's.

| Question | Working answer (confidence) | What the team does |
|---|---|---|
| Board specs (4th email Q1) | Not a decision to withhold the specs; the 403 is Cloudflare's browser-signature block (error code 1010, seen 2026-10-31); saving the Data page's linked file in a browser for the team's own tool is fine (0.8) | 2026-11-01: the captain lifted his two rules for this route; dated notice on the email thread first, save not before 2026-11-04 (UTC), "Save link as", the intake strips the anchors; test-board counts never quoted |
| Final entries (4th Q2) | An unmarked board counts its best official submission; marks are per board and changeable until 2 Dec UTC (0.55) | Mark each board after its last upload: of the entries every ruling allows, the higher published score (a tie or an unscored s2: s1); recheck on 1 Dec IST |
| Validation files as inputs (A20; 3rd Q2, 4th Q3) | Yes (0.88) | No change |
| Gata4 replicates in T3 (follow-up Q3) | Yes (0.8) | Allowed branch for s1; s2 = the refused branch's s1 (no Gata4-derived number), uploaded first (the red team's hedge, approved by the captain on 2026-11-01) |
| Released data across tasks (A15; 18 Oct Q1) | Yes (0.75) | No new cross-modal input in any upload unless the organisers answer; the committed links stay |
| Textbook background (18 Oct Q2) | Fine if qualitative (0.6) | No change: new populations only from the legal route; the T1 scope rule rests on the published E12.5 dissection note |
| truth_cells for the count (18 Oct Q3) | Yes (0.75) | No change |
| truth_cells for size (follow-up Q6) | No (0.65) | No change; the size outputs are committed (f5be247) |
| Generality scope (18 Oct Q4; 3rd Q3) | Every task entered; no eligibility ruling in advance; fits if the summary shows the shared part doing the main work (0.55) | No change |
| Scorer build and parity (follow-up Q1; 18 Oct Q5) | Current metric set, PR #6 not adopted, no endorsement of a parity pair (0.7 for that core; the build itself is near a coin) | No parity hedge; the row-order seed fixed in code |
| Wrong-sign severity (follow-up Q2) | Skill 0.5 (0.65) | No change |
| Duplicated records (follow-up Q4) | Reconstruction artefacts; truth sets as distributed (0.45) | No change |
| z coordinate (follow-up Q5) | As far as the Data Protocol page goes: use the stored z as given; its unit is not stated | No change |
| External staging (follow-up Q7) | Standard references; an embryo either label puts inside a window stays out (0.45) | No change |
| Leak report (follow-up Q8) | Thanks and removal; nothing more from the team (0.55) | No change |
| Generality as corrected (sixth email, 2026-11-10) | No ruling in advance and probably no reply; judged from the method summary after final standing, best claimant first. "Not an unrelated model per task" favours the request; the shared code being the generic part, with each board's main step its own, counts against it (about 0.4 eligible if our request is reached; the AI assistant alone, 2026-11-10) | The request stays, at no cost; no further work on it and no method change for it |
| Stating external sources "with the submission" (not asked, 2026-11-10) | The upload form has no disclosure field, so a private statement to the organisers on the day of the upload, with the method report, meets it (0.6; the AI assistant alone) | Private email on the day of the first official upload: sources and stages per board. The form's free-text fields carry a neutral label only |
| Where the method summary goes (not asked, 2026-11-10) | The report is due 14 days after notification for winners; no route is published for a claimant outside the places, so the organisers' address is the natural one (0.55; the AI assistant alone) | Sent privately in the plan's 24-30 Nov slot unless the screened check of 26-27 Nov finds a published route |

Dropped because they are already answered (see `docs/DATA_POLICY.md`):
- T3 literature: gene-level knockout findings are prohibited, broad phenotype background is allowed (2026-10-15).
- T2 heart extrapolation window equals T1's (2026-10-15).
- External data is allowed for T3 within the windows (GSE208162 ruling, 2026-08-18).
- Scoring conventions: no renormalisation; cell limits inclusive and enforced (2026-10-16).

Deferred to the release (2026-10-27, moved from 2026-10-20; ask only if `panels/index.json` leaves them unclear): each
test board's reference stage (`ref_cells`) and cell bounds. 2026-10-28: the panel files answer HTTP 403 to our client
since 2026-10-27, so question 1 of the fourth email asks how to get the test boards' submission specs.

Not asked because our plan already takes the strict reading: somite-to-day boundary cases (we drop ambiguous embryos),
Tier-B mutant datasets (unused), cross-species foundation models (unused), selection of the second submission among
pre-registered candidates (Rules s.10 already allows it).

---

**Subject:** Human Team: four rules questions before the final phase (cross-task data, developmental background, board metadata, Generality scope)

Hello organisers,

We are a Human Team (captain: Vijayavallabh J). Before the final phase we would like to confirm four points, so that
we follow the rules exactly.

1. **Released competition data across tasks in the final phase.** Once the validation answers are released, may a
   released file from one task be used as training input for another? Specifically: (a) Task 1 scRNA (E8.5, E9.5,
   E10.5) for the Task 2 heart E12.5 board and for Task 3, for example to learn wild-type gene programmes and to calibrate
   scRNA against MERFISH; (b) Task 2/3 MERFISH (heart E8.25 to E10.5, wild type E8.75 and E9.5) for the Task 1 E12.5
   board. Our reading is that released competition data is training data for every task in the final phase, and that
   the "after E9.5 up to and including E13.5" exclusion applies to external data only. Is that correct?

2. **General developmental-biology background inside a protected window.** For Task 1 and the Task 2 heart board we
   would like to use textbook-level qualitative knowledge of heart development between E10.5 and E12.5 (for example
   that the epicardium covers the heart, the coronary plexus begins to form, valve primordia appear and the outflow tract
   septates) to decide which cell populations should appear in our predictions. We would not use any measured values
   from that window (no proportions, sizes or expression data). By analogy with your 2026-10-15 answer about broad
   phenotype background for Task 3, we read this as allowed. Is it? (This overlaps with a question from another team on
   16 October in the public channel.)

3. **Published board metadata.** May we use the `truth_cells` value that `panels/index.json` publishes for a board to
   choose how many cells we submit? Neighbourhood statistics depend on spatial density, so the count matters for the
   Task 2 boards. We will not use the anchor values (floor and ceiling metric values) for anything other than reading
   scores.

4. **Generality Award scope.** If a team enters all three tasks, must the one shared architecture underlie all three, or
   is it enough that it answers at least two of them? And does a shared population-transport model count (unbalanced
   optimal-transport coupling, a conditional velocity/growth network, decoding by retrieval of real cells and one
   calibration protocol, with task-specific encoders, a coordinate head for Task 2 and a genotype head for Task 3)?

Also: (5) will final evaluation use veckit `main` at commit cc3859a (severity_slope, scale_log_ratio), or the metric set
of veckit PR #6? And (6) which time zone applies to the 2 December deadline?

Thank you for running the challenge.

Vijayavallabh J

---

## Follow-up (sent 2026-10-19)

User's decision (2026-10-19): send the full follow-up. Sent by the user on 2026-10-19 by email from
be23b041@smail.iitm.ac.in. Answers: _pending_, no reply as of 2026-10-28 15:41 IST (record them in
`docs/DATA_POLICY.md`). Until question 3 is answered the T3 build makes both branches (allowed and barred); the
fourth email asks it again and names question 1 (the scorer build) as still open.
Duplicate counts come from `experiments/qc/` (update them if the QC run changes them). Points whose answer would help
every team ask for a public answer, so we never hold private information about held-out data.

**Subject:** Human Team: follow-up on the scorer build, Task 3 data use and two data issues

Hello organisers,

Thank you for your earlier answers. A few follow-up points from checking the released files and the scorer. Where an
answer would help every team (points 1, 2, 5 and 6), please feel free to answer publicly instead of to us.

1. **Scorer build and canonicalisation.** Which veckit build will the test boards and the final evaluation use? The
   public validation scores look consistent with the build before the right-handed-frame fix that was merged on
   veckit `main` on 2 October. Under that build a mirror copy scores the same as the original, and `occupancy_dice`
   depends on the sign pattern of the SVD used for canonicalisation, which changes with row order. So two otherwise
   identical submissions can get different Dice values, and a team could deliberately submit opposite "parities". We
   have not done this and will not unless you tell us it is acceptable; we are reporting it so you can decide. If the
   build changes, will the published shape anchors be regenerated?
2. **Severity with the wrong sign.** When the predicted severity slope has the wrong sign (or R² ≤ 0.01), the current
   scorer returns log(1e-3), which equals the floor (skill 0.5). The Evaluation page and the proposal call this case
   "undefined". Which applies in the final evaluation: skill 0.5 or skill 0?
3. **Released Gata4 data in Task 3.** Once the Gata4 validation answers are released, may we use those released
   knockout files as training input for the beta-catenin test prediction? Rules s.10 still lists the Gata4 knockout at
   E8.75 among the held-out genotypes, while the Data page says validation answers become training material in the
   final phase. (We will keep excluding external beta-catenin and Gata4 perturbation data either way.)
4. **Duplicated and interpolated records.** Some files contain the same cell more than once: identical integer count
   rows in adjacent z sections, either at the same xy position (443 in E8.0, 1,164 in E8.25_late, 317 in E8.75) or
   shifted in xy (about 2,460, 560 and 520 more). In the Mab21l2 knockout file, sample_id 35 appears to be stored twice
   and sample_ids 18 and 19 look like one section placed at two z levels (2,811 records together). Separately, E8.0,
   E8.25_late and E8.75 contain 200, 2,558 and 614 records with half-integer counts, each equal to the average of two
   cells in the neighbouring sections. We found none of these in E6.75, E7.25 or E9.5. Are they artefacts of
   registration or section interpolation, and do the validation and test truth sets keep or remove them?
5. **z coordinate.** Is `spatial_3D` z the section position times the section thickness, in the same unit as xy, for
   every file? Were sections registered rigidly, or with expression-aware or non-rigid alignment (some knockout
   sections are not rigid copies of their slide coordinates)? If you can say so publicly: do the test stacks follow
   the same acquisition protocol (section thickness, similar number of sections) as the training files?
6. **truth_cells and tissue size.** Further to our earlier question about using `truth_cells` to choose the submitted
   cell count: may it also inform the predicted tissue size (for example a cell-count-to-volume relation fitted on the
   released files)? We will not use it for that unless you confirm.
7. **Stage conventions for external data.** Could you share your convention for converting somite counts, Theiler
   stages and morphometric stages to E-days, or the somite counts of your own E8.25_late, E8.75 and E9.5 embryos? Two
   concrete cases, which we currently exclude: atlas embryos staged E9.5 by somite count but harvested from a litter
   labelled "E9.75", and embryos staged E13.75-E14.0 morphometrically but harvested from gestational E12.5 or E13.5
   litters. Are they inside the excluded window?
8. **Leak report, addendum.** Further to our private Slack note of 18 October: the same held-out values also appear
   in the veckit 0.1.2 wheel on PyPI and in the description of veckit PR #6. We have not used them.

Thank you,

Vijayavallabh J

---

## Third email (sent 2026-10-21)

Drafted by the AI assistant on 2026-10-20 and saved in the sender's drafts; sent by the user on 2026-10-21 at 01:24
IST (2026-10-20 19:54 UTC) by email from be23b041@smail.iitm.ac.in to virtual.embryo.moonshot@gmail.com. Answers:
no reply as of 2026-10-28 (record them in `docs/DATA_POLICY.md`; question 2 decides assumption A20, question 3 A21).
Question 1 is answered by the final-phase Rules s.11 (P3 opened 2026-10-27 at 09:00 America/Los_Angeles, 16:00 UTC).
Question 2: the user read the final-phase Rules s.10 and Data page as allowing it (A20 passed, 2026-10-28); the
fourth email asks the organisers to confirm. Text as sent, rewrapped:

**Subject:** Human Team: three questions before 27 October, and two open emails

Hello organisers,

Thank you for the two updates of 20 October. Three short questions remain before the test phase:

1. Start time. At what time, and in which time zone, does the test phase open on 27 October?

2. Validation answers as training data. Once released, may the validation answers be used as training data and as
   model inputs for the test predictions (for example E10.5 when predicting E12.5)? The Data page says they become
   training material, while Rules s.10 still lists those stages as held out. Another participant asked the same in
   #virtual-embryo-challenge on 19 October.

3. Generality Award. Does this count as one shared architecture: a single operator that builds every prediction from
   real cells (composition re-weighting, then count-space expression edits, with coupling and interpolation between
   stages where a task needs it), with task-specific heads only (a 32,285-gene or 500-gene encoder, a coordinate head
   for Task 2, a genotype head for Task 3)? This repeats question 4 of our 18 October email.

Our emails of 18 October ("four rules questions before the final phase") and 19 October ("follow-up on the scorer
build, Task 3 data use and two data issues") are still open, as is our private Slack note of 18 October to Weize Xu
and Zehua Zeng. Two points from them affect what we do on 27 October: which scorer build the test boards use, and
whether the released Gata4 files may be used in the beta-catenin prediction.

Best regards,
Vijayavallabh J
Team captain, Human Team track

## Fourth email (sent 2026-10-28)

Drafted by the AI assistant on 2026-10-28 at the user's request and saved in the drafts of be23b041@smail.iitm.ac.in;
sent by the user on 2026-10-28 at 15:41 IST (10:11 UTC) to virtual.embryo.moonshot@gmail.com, with the signature
changed to the team name ("Team captain, Mithril"). Answers: _pending_ (question 3 confirms A20 and selects the T3
branch, `docs/P3_RUNBOOK.md` section 5). Text as sent, rewrapped:

**Subject:** Human Team: board specs blocked (HTTP 403), final-entry nomination, use of the unlocked validation files

Hello organisers,

Thank you for opening the final phase. We unlocked the P2 validation release on 28 October and have the six files.
Three questions:

1. Board specs. Since 27 October, https://virtualembryo.ai/challenge/panels/index.json answers HTTP 403 (Cloudflare)
   to our scripted download, while the Data page still links the panel files and mentions no outage. We have not tried
   to get around it. Is the block intended, and how should a team get the test boards' submission specs (panel genes,
   cell-count limits, required obsm keys)? A spec-only file served through the API with the team key, like the data
   files, would work for us.

2. Final entries. The account page says "In the final phase the captain picks which submission counts", while Rules
   s.11 says a board ranks a team by its best official submission. If the captain marks nothing on a board, does the
   best official submission count? Is the pick made per board, and can it be changed until 2 December?

3. Validation files as model inputs. Rules s.10 and the Data page allow the unlocked files for P3 preparation. We read
   this as also allowing them as inputs when we produce the test predictions: for example the T1 and T2 heart E10.5
   files as the starting stage for E12.5, and the two Gata4 knockout replicates in the beta-catenin prediction. Please
   tell us if that reading is wrong. This would also answer question 2 of our 21 October email and question 3 of our
   19 October follow-up.

The start time we asked about on 21 October is now in Rules s.11; thank you. Still open from our 19 October
follow-up: which veckit build scores the test boards.

Best regards,
Vijayavallabh J
Team captain, Mithril

## Fifth email (sent 2026-11-05)

Drafted by the AI assistant on 2026-11-05 at the captain's request ("yes draft the organiser email") in the drafts of
be23b041@smail.iitm.ac.in. The captain sent it on 2026-11-05, by 18:18 UTC. The sent copy is addressed to
virtualembryomoonshot@gmail.com, the same inbox. It
asks again question 2 of the 19 October follow-up, which decides T3's second upload (`ASSUMPTIONS.md` A24;
`experiments/research_2026_11/REPORT.md` section 3). Answers: _pending_. Text as drafted:

**Subject:** Human Team: how the final evaluation scores an undefined severity_slope (Task 3)

Hello organisers,

One short scoring question about Task 3. Its answer affects every team the same way.

In Task 3, the severity metric is undefined when the predicted slope has the wrong sign, beta <= 0, or R^2 <= 0.01.
The current scorer then returns log(1e-3), which equals the floor, so the live leaderboard gives that metric skill
0.5. The Evaluation page and the competition proposal call this case "undefined", which would map to skill 0.

Which will the final evaluation of the test boards use: skill 0.5 or skill 0?

This was question 2 of our 19 October follow-up. A public answer in #virtual-embryo-challenge would be welcome, since
it applies to all teams.

Best regards,
Vijayavallabh J
Team captain, Mithril

## Sixth email (sent 2026-11-10)

Drafted by the AI assistant on 2026-11-10 on the captain's "go ahead with your recommendation" (`docs/PLAN.md` decision
log 2026-11-10 (u)), in the drafts of be23b041@smail.iitm.ac.in, to virtual.embryo.moonshot@gmail.com. The captain sent
it on 2026-11-10 (the captain's words in the main session: "i sent the email"; the draft had left the drafts folder when
checked, and the sent copy was not opened, so any edit before sending is not recorded). Answers: _pending_ (record them
in `docs/DATA_POLICY.md`); the captain, 2026-11-10: "so far they have not replied to a single email", so a working answer stands in
the table above (decision log 2026-11-10 (w)). It corrects the architecture described in
question 4 of the 18 October email and question 3 of the third email (the D3 review, decision log 2026-11-10 (s)). Text
as drafted:

**Subject:** Human Team: correction to our Generality Award question (team Mithril)

Hello organisers,

This corrects how we described our entry when we asked about the Generality Award (question 4 of our 18 October email
and question 3 of our 21 October email). Neither description matches what we built, so please do not judge our request
on them.

The 18 October text mentioned unbalanced optimal transport, a conditional velocity/growth network and one calibration
protocol. Our entries use none of these.

The 21 October text said a single operator builds every prediction, with task-specific heads only. That overstates what
is shared. What our tasks have in common is the design: every submitted cell is a real released cell, and its expression
is edited only in count space or not at all. On Task 1 and both Task 2 boards, shared library functions draw the cells
and make the count edits, and each task's code runs the steps it needs; no single function runs on more than two of our
four boards. Task 3 is built by a script we pre-registered before that library existed. It calls none of the library's
functions, though five of its own are identical to the library's copies, which a test enforces.

Our method summary will set this out and ask to be considered for the award. Does this structure count as one shared
architecture? If you judge that it does not, we will withdraw the request, and our entries would not change.

Please reply privately, as this describes our method.

Best regards,
Vijayavallabh J
Team captain, Mithril

## Seventh email (sent 2026-11-13)

Request R-C of `docs/PLAN.md` decision log 2026-11-12 (c). On 2026-11-13 the captain forwarded a side agent's note to the
main session ("If you want it, ask the main agent to draft R-C now"); the AI assistant read the forwarded note as that
request and drafted the email in the drafts of be23b041@smail.iitm.ac.in, to virtual.embryo.moonshot@gmail.com. The
captain sent it on 2026-11-13 (the captain's words in the main session: "i also sent the draft"; the draft had left
the drafts folder when checked, and the sent copy was not opened, so any edit before sending is not recorded).
Answers: _pending_ (record them here and in `docs/DATA_POLICY.md`). A reply is not pasted into the main session: a
screening agent reads it first and passes on the procedural answer only, as the classification below requires.
Screened read of 2026-11-18, about 15:40 to 15:46 UTC (the section at the top of this file): no answer found in the
mailbox, the Slack channel or the Data Protocol page.
Screened read of 2026-11-19, about 04:38 UTC ("Screened reads of 2026-11-19" at the top of this file): no answer
found in the organiser mail thread, the Slack channel (threads included) or the Data Protocol page.

**Classification, stated before any answer.** What it asks is procedure: which released sample a held-out sample was
processed with, and whether a protocol was the same. That is the class of the collection text released through the
screen on 2026-11-10, not a measured property. A reply that carries a figure about a held-out sample (a count, a depth,
a size, a share) is not read into any file or decision; using one would need the captain's named exception, fixed
before the read, as on 2026-11-11. The captain confirmed this classification by sending the email.

**Where an answer can appear, and who reads it** (added 2026-11-13 after a side agent's note to the captain). The
email lets the organisers answer by reply or with a line on their Data Protocol page, and they have also answered in
public before. The captain approved the two reads on 2026-11-13 ("yes to the screened reads on 19 and 23 Nov"); their
brief, written before any read, is `docs/SCREENED_READ_SEVENTH_EMAIL.md`. None has run yet.
- A reply by email reaches the captain's inbox. The captain tells the session that one has arrived and does not paste
  it; a screening agent reads it and passes on the procedural answer only.
- A line on the Data Protocol page. Only a screening agent reads that page, at the captain's say-so, and the screened
  check of 26-27 Nov did not name it. That check is also too late for a Task 1 build change, which has to be
  rehearsed by about 24 Nov. The AI assistant therefore asks the captain for two screened reads of the page, about
  19 and 23 Nov.
- A post in #virtual-embryo-challenge, read through the same screened check as the thread of 2026-11-11.
- Scope of every such read, fixed before any of them: the screen passes on only a statement of which released sample
  the E7.75 sections or the E12.5 Multiome libraries were processed with, or of the same or another protocol, with
  every number replaced by [n]. Everything else stays withheld: any measured or relative property of a held-out
  sample, anything about other held-out samples, and any other change to the page beyond the fact that it changed.

**The review.** Three reviewers (accuracy, strategy, tone; single agents given only the draft, the recorded quotes and
the stakes) each returned "send with edits". What changed:
- A Task 3 question (was the beta-catenin knockout processed with the Gata4 replicates) was cut: both T3 uploads are
  fixed by frozen rules and the better one counts, so no answer changes what the team does. The AI assistant had added
  it to R-C's two questions; the captain can restore it.
- "Run" alone would get a useless answer, since a slide is usually its own imaging run. The embryo question now names
  the slide and the hybridisation and imaging batch, and asks which stage, without offering only E8.0 and E7.5.
- "The dissected hearts" as a name for E10.5's third library was our own inference, beside the page's "E10.5 five hearts
  per replicate". The library is named by its cell-name suffix alone (`_3`, checked on the released file's names).
- The Task 1 question separates library preparation from the sequencing run and asks about the protocol.
- The guard against figures is explicit, a date gives a reason to answer now, and a public answer is offered as the
  organisers' option.
The strategy reviewer put the chance of a useful reply before the freeze at 5 to 10%.

**How an answer is read, fixed before any answer.** Nothing changes without the captain's ruling.
- Embryo, "processed with E8.0, E8.25_late or E8.5" (the released stages that share one profile, diagnostic N23): the
  advice on R-A becomes option 1, `fgw_aligned_series`. On the three practice tasks built that way it scores +6.69,
  +3.56 and +3.42 against today's pair (tries E28).
- Embryo, "processed with E7.5": the advice becomes to keep today's fgw_shape0. On the practice task of that kind every
  second upload scores 0.00 against today's pair. A specialist for E7.5's profile (the later bracket's cells rendered to
  it) is not built; it would be registered as a try before any use.
- Embryo, "a batch of its own", or with any other stage: the advice becomes to keep today's fgw_shape0. On the two
  practice tasks whose target has a profile of its own, option 1 scores +0.86 and -1.42 and option 2 +0.77 and -1.14.
- Embryo, "not recorded" or no answer: R-A stays as `docs/PLAN.md` decision log 2026-11-13 (a) puts it.
- (Since 2026-11-13, row R6, T1's first upload is already drawn from library 3; the two Task 1 lines below are kept as
  they were fixed. An answer "with library `_3`" or "the same protocol" now confirms the built file; an answer
  "separately and with another protocol" goes to the captain with try T30's rows, and nothing changes by itself.)
- Task 1, "prepared with library `_3` and sequenced with it", or "the same protocol": R-B goes back to the captain with
  the advice changed to yes, read on try T30's rows at library 3's depth alone. There the switched pair against today's
  pair is +1.60, 0.00 and +8.24 at a true cardiomyocyte share of 0.45, 0.60 and library 3's own. This is the organiser's
  procedural answer that decision log 2026-11-12 (c) names as a clean way to reopen R-B.
- Task 1, "separately and with another protocol", "not recorded" or no answer: R-B stays declined.

Text as drafted:

**Subject:** Human Team: processing batches of two held-out samples (procedure only)

Hello organisers,

A short procedural question about two held-out samples. Please leave any figures about held-out samples (counts,
depths, sizes) out of the answer: we cannot use them, and a word or two per item is plenty.

The Data Protocol page says that knockout and control embryos were processed in separate runs. As with any assay,
samples processed at different times can differ for reasons unrelated to stage, and we would like to account for that
correctly.

1. E7.75 (Task 2 embryo): were its sections processed together with those of any released MERFISH stage, meaning on the
   same slide or hybridised and imaged in the same batch? If so, which stage?

2. E12.5 (Task 1): were its Multiome libraries prepared together with the released E10.5 library whose cell names end in
   "_3", and sequenced in the same run? If not, were they made and sequenced with the same protocol as that library?

An answer before 27 November would reach our final entries. If you would rather answer every team, a line on the Data
Protocol page would suit us equally.

Best regards,
Vijayavallabh J
Team captain, Mithril

## Planned eighth email: how final entries are marked (drafted 2026-11-13; not sent; the captain decides)

Why now. Question 2 of the fourth email (28 October) has no answer on record, and the sixth campaign round showed
that it decides more than how the marks are made: the value of every second upload as "the better of two" holds only
if the pick is per board (`docs/S2_RULES.md`, "What an s2 is worth"), and a cold check on 2026-11-13 stopped an option
for the embryo's first upload on exactly that point (`docs/PLAN.md` decision log 2026-11-13 (hh)). The AI assistant
drafted this on its own initiative after that check; the captain has not asked for it and decides whether and when it
goes. Two figures in the team's files differ and are left as they are: `docs/S2_RULES.md` (hash-pinned) rates
"per board" at 0.55, the table of 31 October; the table of 2026-11-10 above rates the whole final-entries row at 0.7.

Review. A three-lens panel (wf_11847608-0af: accuracy, strategy, tone; read-only, the repository's own files only)
returned "send with changes" three times. What changed: the opening no longer says the question "is still open",
which implied the email's other two questions had been answered; the reason given is that the answer applies to every
team the same way, not that it decides how second uploads are planned; "while Rules s.11 says" became "and"; question 1
names what a pick over several boards could be; question 2 is open ("what counts there?") and no longer offers the
answer the team assumes. On timing and channel the strategy lens advised: send it on its own, soon (an answer helps
only before the first uploads on 27 November), as a new short email and not on the old thread or in the Slack thread
on Task 3 severity; if nothing has come and the Final column has not appeared by about 24 November, the captain may
post the same three questions as one neutral top-level message in #virtual-embryo-challenge, which is his own call and
brings the runbook's public ask of 28 November forward. Expect no reply (none of seven emails has one); nothing in the
plan waits on it. The captain, 2026-11-14, on being reminded that this draft and the Task 3 batch question wait: "so
far organizers have not answered any of my single mail too". The draft stays unsent unless the captain decides
otherwise; the public message of about 24 November stays the captain's own call.

Readings fixed beforehand (they stay in this file; none is sent):
- "Per board", alone or with the other points: no change. A point left unanswered keeps its working answer.
- "One pick covers all of a team's boards" (a set or bundle): the per-set clause of `docs/S2_RULES.md` applies as
  written (the second uploads still stand; no second set is kept), and the marks follow the runbook's rule for that
  case ("the bundle with the higher overall standing", `docs/P3_RUNBOOK.md` section 5). The AI assistant's own
  addition, which binds nobody until the captain says so: it then proposes no file for a first-upload slot that is
  far below the entry on most practice tasks.
- "One pick per task" (the two Task 2 boards together, Task 1 and Task 3 each alone): the per-set clause is applied
  to the two Task 2 boards only, and the captain is told the same day.
- "The team may choose, per board or as a bundle": treated as per board; every board is marked on its own.
- An unmarked board "counts its best": the working answer; the captain still marks every board as soon as its last
  official upload is scored. "Counts nothing", "the first" or "the latest": the same marking, with the recheck on
  1 December (IST) no longer optional. "Marks are ignored, the best always counts": no change.
- "Picks cannot be changed", or close before 2 December: each mark is made only after that board's last upload is
  scored and before the stated close; T3's case in `docs/S2_RULES.md` (s1 unscored by 30 Nov IST: mark s2) goes to the
  captain that day, since such a mark could not be undone.
- A one-liner that does not say which ("yes"; "the captain picks one"): counts as no answer; one clarifying reply at
  most, at the captain's choice.
- A reply that cites the notice of 31 October as agreed ("as you said, one per board"): recorded as a per-board
  answer.
- An answer after 27 November: the first uploads are designed by then; only the readings about marking apply.
- An answer that changes the rule for every team: the captain is told the same day; the per-set clause is applied
  without new tries.
- No answer by 27 November: the working answers stand.

Where it is. A draft with the text below is in the drafts of be23b041@smail.iitm.ac.in, created by the AI
assistant on 2026-11-13 without being asked, with the two signature lines left as placeholders so that it cannot
go out unread. Sending it, changing it or deleting it is the captain's alone. Answers: none can exist yet.

Text as drafted (to virtual.embryo.moonshot@gmail.com, as a new email):

**Subject:** Human Team: how final entries are marked (one question from 28 October)

Hello organisers,

We are asking one question from our email of 28 October again, on its own; the answer applies to every team in the
same way.

Final entries. The account page says "In the final phase the captain picks which submission counts", and Rules s.11
says a board ranks a team by its best official submission.

1. Is the pick made separately for each board, or does one pick cover several boards together (both boards of a task,
   or all of a team's boards)?

2. If the captain marks nothing on a board, what counts there?

3. Can a pick be changed until 2 December (UTC)?

A one-line answer would be enough. If you would rather answer every team at once, a line on the Rules page or in
#virtual-embryo-challenge would serve as well.

Best regards,
Vijayavallabh J
Team captain, Mithril

## Disclosure email for the uploads of 18 November (drafted 2026-11-18; in the captain's drafts; the captain sends)

The runbook's email of sources and hashes (`docs/P3_RUNBOOK.md`, "Disclosure ids for each official upload"), which was
due before the first official upload and follows it: the captain uploaded by hand on 2026-11-18 (`docs/PLAN.md`
decision log 2026-11-18 (j)). Drafted by the AI assistant; reviewed for accuracy, strategy and tone by three
independent agents, rewritten, verified by a fourth, and its later additions (the two second uploads it lists)
checked by a fifth. Put into the drafts of be23b041@smail.iitm.ac.in, to virtual.embryo.moonshot@gmail.com, at the
captain's word of 2026-11-18 (decision log (k)). The runbook also has it copied by Slack DM to Weize Xu; that is the
captain's to do. Whether and when it was sent, and any edit before sending, are recorded here when the captain says.
Answers: none can exist yet.

**Amended on 2026-11-19 (wf_038f7126-b04: a writer, three reviewers, a fixer; not yet in the captain's drafts, whose
copy still has the text below and is not to be sent as it stands).** Eight changes and nothing else: (1) the subject
ends at "our official uploads of 18 November"; (2) the opening sentence no longer says "and for the two second
uploads that follow it"; (3) to (5) the paragraph "Second uploads that follow this note" and its items 5
(`T1_s2.h5ad`) and 6 (`T2embryo_s2.h5ad`) are removed: the first is not uploaded, the second is not decided; (6) in
the X002 paragraph, after "no absolute expression level.", one sentence is added: "Those cell types are the authors'
labels; a separate note sent with this one is about them."; (7) the sentence on 19 October ends at "of its E14.5
hearts.", without "the share in our second Task 1 upload lies within that range"; (8) "What comes next" begins
"Before each of our second uploads we will send its board, file hash and sources in a short note." Left for the
captain: the validation paragraph still says "before our heart and Task 3 second uploads"; the last T1 file also
uses the released E10.5 file, so he may widen it to "before our second uploads".

**Three more changes, and the amended note put into the captain's drafts (2026-11-19, 11:25 UTC by the draft's own
time stamp; at his word of that day, decision log 2026-11-19 (f)).** Before replacing the draft the AI assistant
found one sentence the eight changes had left false: "no measured data from the held-out stages or genotypes is
used in any upload by any route" cannot stand beside the labels email. A further agent reviewed the proposed edits
in three passes (accuracy, strategy, tone) and reworded the first: (9) that sentence now reads "Apart from the
official validation files above, used under Rules s.10 as our email of 28 October described, and apart from the
[X002] authors' cell-type labels and processing that our separate note describes, no measured data from the
held-out stages or genotypes is used in any upload by any route, and no gene-level knockout finding from outside
the released data."; (10) the validation paragraph ends "please tell us before our second uploads.", since every
board's second upload uses an official validation file; (11) the register sentence begins "Our Task 1 upload above
uses X001, X002, ...", since the remaining T1 upload uses no external dataset. The reviewer asked whether the X006
sentence "That figure no longer supports any choice" still holds: it does (the atlas share that supported keeping
E10.5's own mix was X002's CD1 share, audit item 7, not X006's figure). The draft in his folder now holds the text
below with changes (1) to (11), read back against it; the earlier copy that named `T1_s2.h5ad` is gone. One thing
for him to glance at before sending: how the veckit address in the X001 line is shown (the mail tool's read-back
wraps it in a redirect link). Text as drafted on 18 November:

Subject: Human Team (Mithril): sources and file hashes for our official uploads of 18 November and two second
uploads

Hello organisers,

This note lists, for each of our four official uploads of 18 November 2026 and for the two second uploads that
follow it, the file's SHA-256 and every external source used, as Rules s.10 asks. It should have reached you before
the first uploads; it follows them. As it describes our method, please keep it private.

The uploads (18 November 2026, times in UTC as the site lists them)

1. T1:test, 14:37. File T1_s1.h5ad, SHA-256 a294062f5f7bbb2dae19ef695979cc21f03f38df40d22be9ccdd30af36e3bdef.
2. T2:embryo:test_interp, 14:38. File T2embryo_s1.h5ad, SHA-256
  dc910d4b3fb1d62b643bf96e5765741f8f0ad011f7b2d928a7955928ddee1303.
3. T3:bcatenin, 14:38. File allowed_s1.h5ad, SHA-256
  3d6c94a97b02b8f378b275070dad2767e54005e4194d093e4e01bdf3881b4203.
4. T2:heart:test_extrap, 14:42. File T2heart_s1_half.h5ad, SHA-256
  2ca632d399502943da3f50d6609ab630927692073521735769fdacba0ea3a5a5. At 14:38 the site refused the file our procedure
  had built (T2heart_s1.h5ad, SHA-256 3d6548382b24e01445148d238b5b9ec832f09c4ef6ccddbf819ccfb6ade7c46c) for its cell
  count; for your information, the same file passed the site's format check on 13 November. The uploaded file is a
  random half of that file's cells (drawn without replacement with a fixed seed, 1000001, cell order kept) and is
  otherwise identical. We halved it, rather than cutting it to the board's limit, so that the number of cells
  follows from the built file alone. This is the one change to a cell number since our note of 4 November, and the
  method report will list it with the other later decisions that bear on the quantities that note covers.

Second uploads that follow this note. Both files were built on 13 November, before any score, by the rules we had
fixed in advance for our second uploads.

5. T1:test, second upload. File T1_s2.h5ad, SHA-256
  f169ea492e8dbac7e8b298f7dfa13fcceba71b5ff35e12ccfdc3a9504346304f. Sources as for upload 1.
6. T2:embryo:test_interp, second upload. File T2embryo_s2.h5ad, SHA-256
  41982d8afd784ecb07a4021077f701c3806c3bdb06eda92c2cb99fcc41b19b73. Sources as for upload 2.

External sources

Codes in brackets are the rows of our internal source register. Task 1 uses X001, X002, X003, X005, X006, X007, X008
and X009; both Task 2 boards and Task 3 use X001, X004 and X005.

On every board:
- [X001] veckit from its public repository (github.com/aristoteleo/veckit, commit cc3859a, MIT licence): local
  scoring only.
- [X005] Open-source Python libraries, pinned in our lock file (as locked on 5 November 2026): POT 0.9.7 and rdata
  (MIT); numpy 2.5, scipy 1.18, scikit-learn 1.9, pandas 2.3, anndata 0.13, scanpy 1.12 and h5py 3.16
  (BSD-3-Clause). Library calls only; no data or pretrained weights come with them.

On both Task 2 boards and on Task 3, also:
- [X004] One function copied from veckit 0.1.1 (commit 46d41e6, common/shape_metrics.py, _canonicalise; MIT licence)
  into our local scorer. Never part of a prediction.

On Task 1 only:
- [X002] Feng et al., Nat Commun 13:7960 (2022), GEO GSE193346: wild-type CD1 hearts at E9.5 and from E14.5 to P9.
  Cells of E10.5 to E13.5 were removed on the per-cell stage column. It supplies a per-cell-type direction of change
  (E9.5 to E14.5) and no absolute expression level. For completeness: GEO's full object (all stages) and the full
  metadata tables of both strains in the study, CD1 and C57BL/6, were on our compute machine from 18 to 19 October
  2026 and were then deleted. While they were there we printed once the median library size and mitochondrial
  fraction of cells from the removed stages; nothing in any upload was set from that print. On 19 October, while
  those tables were there, we computed each strain's cardiomyocyte share of heart cells from the cell-type and stage
  labels alone (no expression) of its E14.5 hearts; the share in our second Task 1 upload lies within that range. On
  1 November the reduction was rerun from the raw GEO files by code that read the CD1 metadata table as a stream,
  parsed the object in memory, kept only the counts of the E9.5 and E14.5-and-later cells, and deleted the raw
  files.
- [X003] Qiu et al., Nature 626:1084 (2024): embryos 1 to 22 (labelled E8.5 to E9.5) and 54 to 74 (E14.25 to P0).
  Diagnostics only, not a model input. Expression of embryos 23 to 49 was never downloaded; embryos 23-24 and 50-53
  were dropped as boundary or in-window litters. A table of labels and metadata for all 74 embryos (no expression)
  was on our compute machine from 18 to 19 October 2026 and was then deleted; nothing was set from it.
- [X006] Feulner et al., Biology Open (2026), doi:10.1242/bio.062493, whose E14.5 hearts come from Wünnemann et al.,
  Nat Genet 52:40 (2020): on 19 October our research took a cardiomyocyte share of heart cells from the cell-type
  labels of its E14.5 whole hearts (the authors' E14.5; CELLxGENE codes the sample as Theiler stage 22). Our notes
  match the figure to this source without settling it; the alternative is the Qiu et al. labels above. That figure
  no longer supports any choice, and no cell, count or expression from this source enters any upload.
- [X007] van den Brink et al., Nat Methods 14:935-936 (2017): two published gene lists (Supplementary Tables 1 and
  5), from adult mouse muscle satellite cells. Gene symbols only; no cell, count or expression.
- [X008] Bakken et al., PLoS One 13(12):e0209648 (2018): one published per-gene table (S6 Table), from adult mouse
  visual cortex. Gene properties only; no cell, count or expression.
- [X009] Young and Behjati, GigaScience 9:giaa151 (2020): a published statement on haemoglobin genes in ambient RNA
  (a method paper, no stage), which informed one step.

Official validation files. Each of the uploads above also uses official validation files that we unlocked on 28
October: Task 1 the E10.5 file, the heart the E10.5 and E8.5 files, the embryo the E7.5 file, and Task 3 the two
Gata4 knockout replicates. This is the use we asked about in question 3 of our email of 28 October and announced in
our note of 31 October. That question has not been answered yet; if our reading of Rules s.10 is wrong, please tell
us before our heart and Task 3 second uploads.

No pretrained model is used in any upload. Apart from the official validation files above, used under Rules s.10 as
our email of 28 October described, no measured data from the held-out stages or genotypes is used in any upload by
any route, and no gene-level knockout finding from outside the released data.

One correction to our note of 4 November. That note said our cell counts and sizes follow rules committed before the
data release. For Task 3, the rule for the number of cells (use the board's maximum) was fixed on 1 November 2026:
after the release, though still before what we saw on 4 November. The method report will say the same.

What comes next. Before the heart's and Task 3's second uploads we will send their board, file hash and sources in a
short note. Our method report will set out the procedure we fixed in advance for the second uploads and any change
we make to it after the first scores.

Best regards,
Vijayavallabh J
Team captain, Mithril

## Scorer-build note (drafted 2026-11-18; in the captain's drafts; the captain decides)

A repeat of follow-up question 1 of 19 October (the scorer build and the handedness of the canonical frame), raised
again after the first Task 2 scores. Drafted by the AI assistant and reviewed for accuracy, strategy and tone by
three independent agents, whose changes are in the text: it asks only which version applies and asks for no change
to any score, because under the right-handed version every file is compared directly, which could score our own
embryo file lower as easily as higher; it names no board, hash or score of ours; and it states beforehand that no
file or row order of ours is chosen by its handedness. Put into the same drafts folder at the captain's word of
2026-11-18. Sending it, changing it or deleting it is the captain's alone. Answers: none can exist yet. Text as
drafted:

Subject: Human Team (Mithril): which canonicalisation will score occupancy_dice in the final evaluation (Task 2)?

Hello organisers,

One scoring question about occupancy_dice in Task 2. The answer matters to every team with a Task 2 entry.

In veckit before the right-handed-frame fix (for example at commit 46d41e6), common/shape_metrics.py::_canonicalise
rotates each cloud into its own PCA basis and leaves the handedness of that basis to the SVD. occupancy_dice then
compares the prediction with the truth either directly or as mirror images. Which case applies is set by the SVD
sign patterns of the two clouds as the scorer subsamples them, so for a given truth it changes with the order and
the number of rows in a submitted file. It does not depend on the handedness of the predicted shape: under that
version a mirrored copy of a file scores exactly like the file. On a released heart stage (E9.5), comparing two
random halves of the stage with each other, the mirror case scored on average about 0.10 to 0.14 lower in raw Dice
than the direct case; the size of the difference depends on the stage's shape. veckit main at commit cc3859a (2
October), the version we pin for local scoring, makes the basis right-handed, which removes the effect.

In our 19 October email we wrote that the public validation scores looked consistent with the earlier version. When
the validation answers were released to us, we checked the published floor and ceiling Dice values of the validation
boards: all of them come out exactly the same under both versions, because on those boards the SVD happened to
return matching handedness. So those values do not show which version the server runs.

Two questions:

1. Which version of _canonicalise returns the test-board scores now, and which will the final evaluation use?
2. If it is the earlier one, do you plan any change before the final evaluation? We are not asking for any change to
  our own scores; we only need to know which version applies.

For transparency: our first Task 2 uploads of 18 November are what made us raise this again; we do not know whether
the handedness affected any of their scores. We have not submitted the same prediction in two row orders, and no
file or row order of ours is chosen by its handedness. Our heart upload is a random half of the file we built, for
the reason given in our note on sources of the same day, and our second Task 2 uploads are different predictions
chosen on other grounds. We will not submit the same prediction in two row orders unless you say that is acceptable.

We asked which build the final evaluation uses on 18 October and raised this handedness question on 19 October
(question 1 of that follow-up). A public answer in #virtual-embryo-challenge would be welcome, since it applies to
all teams.

Best regards,
Vijayavallabh J
Team captain, Mithril

## Notes before Task 1's and the embryo's second uploads (drafted 2026-11-19; not in the captain's drafts; the captain decides and sends)

The amended sources note promises "Before each of our second uploads we will send its board, file hash and sources
in a short note." The heart's and Task 3's notes are below, in the captain's drafts since 18 November. These two were
drafted by an agent of the AI assistant, reviewed by a second in three passes (accuracy, strategy, tone: four points)
and fixed by a third (wf_e7336078-e77). They go into his drafts at his word. Task 1's hash is the rehearsed build's
(`experiments/builds/results/t1none_385ff58.txt`); the embryo's file name is the set embryo2's (decision log
2026-11-19 (j)), and its hash is filled when that set is built.
Answers: _none; not sent_.

**Subject:** Human Team (Mithril): second Task 1 upload, file hash and sources

To: virtual.embryo.moonshot@gmail.com

Hello organisers,

As our email on sources and file hashes of 18 November promised, here are the hash and sources of our second Task 1 upload, before we make it. Please keep this private.

T1:test, second upload. File T1_s2none.h5ad, SHA-256 55eaee547ff18d7477da5b730c34e7de9ab74f85cc12fd01a9eff6450367611c. It is our first Task 1 upload's method without its one step that used external data. It reads no external dataset.

Sources: X001 (local scoring only) and X005 (pinned libraries). The file reads nothing of X002 (Feng et al.), nor of X007 to X009, which acted only on the removed step. We still list X002: two parts kept here were chosen with evidence that used its labels, as our labels note explains. X003 and X006 never entered a model. Official validation file: E10.5, as before.

It replaces the second file fixed in advance; the method report will say so.

Best regards,
Vijayavallabh J
Team captain, Mithril

**Subject:** Human Team (Mithril): second embryo upload, file hash and sources

To: virtual.embryo.moonshot@gmail.com

Hello organisers,

As our email on sources and file hashes for the uploads of 18 November promised, here are the hash and sources of our second embryo upload, before we make it. Please keep this private.

T2:embryo:test_interp, second upload. File T2embryo_s2x.h5ad, SHA-256 a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9. It has the cells and expression of our first embryo upload, with each cell placed at a position measured at E7.5 in place of a moved position. Sources as for our first embryo upload (X001, X004 and X005), and the same official validation file (the embryo's E7.5 file).

This replaces the file we had fixed in advance. We built and chose it after our first scores, which the site returned to us as to every team; no number in the file was computed from a score. The method report will set out this change.

Best regards,
Vijayavallabh J
Team captain, Mithril

## Note before the heart's second upload (drafted 2026-11-18; in the captain's drafts; the captain sends)

The short note the disclosure email promises before the heart's second upload: the board, the file's hash and its
sources, and what the file changes against the procedure fixed in advance (`docs/S2_RULES.md`, amendment of
2026-11-18). Drafted by the AI assistant and reviewed for accuracy, strategy and tone by three independent agents,
whose changes are in the text: it names every departure (the size mode and its figure, the dropped expression step
and with it the read counts, the scope and the cell number against the file fixed for this upload), says that the
choice followed the first scores and the published metric definitions, and claims only what the record backs (both
size figures were recorded on 28 October; no number was computed from a score). It goes only after the disclosure
email has gone. Sending it, changing it or deleting it is the captain's alone. Answers: none can exist yet. Text as
drafted:

Subject: Human Team (Mithril): second heart upload, file hash and sources

Hello organisers,

As our email on sources and file hashes for the uploads of 18 November said we would, here are the file hash and
sources of our second heart upload, before we make it. As it describes our method, please keep it private.

T2:heart:test_extrap, second upload. File T2heart_s2iso_half.h5ad, SHA-256
dc436d9225bb28e7997d2ec51e95254cd125e569d44c059596298d0a7b6af5c2. The external sources are the same as for our first
heart upload (X001, X004 and X005 in that email), and so are the official validation files used (the heart's E10.5
and E8.5 files).

This is not the file we had fixed in advance for this upload, which trimmed the tissue to a heart-focused half and
sized it for that scope. We decided on this one after our first scores, which the site returned to us as it does to
every team, and from the published metric definitions. No number in the file was computed from a score: the size
figure and the halving rule were fixed before the scores, and the scores decided only which of our existing steps to
keep.

Compared with our first heart upload, it differs in two steps. First, the size: the cells are scaled evenly on all
three axes to the figure our method had computed for that way of scaling, where the procedure we had fixed in
advance enlarged the two in-plane axes and held the third, at a slightly different figure. Both figures were
computed and recorded on 28 October. Second, the expression step of our first upload is dropped: the cells keep the
E10.5 file's released expression and read counts. Its scope and its number of cells (a random half, by the rule of
our first heart upload: fixed seed 1000001, cell order kept) are the first upload's.

The method report will list each of these changes against the procedure we fixed in advance.

Best regards,
Vijayavallabh J
Team captain, Mithril

## Note before Task 3's second upload (drafted 2026-11-18; in the captain's drafts; the captain sends)

The short note the disclosure email promises before T3's second upload. Drafted by the AI assistant and reviewed for
accuracy, strategy and tone by three independent agents, whose changes are in the text: it says "computed from"
where the draft said "comes from", and names the Gata4 replicates as practice tasks for the two choices made on
them (the knockout factor left out by the rules of 19 October, and the wild-type correction the captain added on
13 November by his own decision); it says that the order was reversed before any Task 3 score and that the file is
unchanged; it names the wild-type stage files the correction is fitted on and that three of them are official
validation files; it reports the practice try of 18 November (rows K16 and K17), which changed nothing; and it
explains the file's name. Both notes now refer to the disclosure email by its subject and not by a sending date.
Before sending: the file chosen in the upload form must have the hash the note states (the archived copy, not a
local rebuild). It goes only after the disclosure email has gone. Sending it, changing it or deleting it is the
captain's alone. Answers: none can exist yet. Text as drafted:

Subject: Human Team (Mithril): second Task 3 upload, file hash and sources

Hello organisers,

As our email on sources and file hashes for the uploads of 18 November said we would, here are the file hash and
sources of our second Task 3 upload, before we make it. As it describes our method, please keep it private.

T3:bcatenin, second upload. File refused_s1.h5ad, SHA-256
0394575a11fcd4f915e0daf4972e5a783b6b89efa065ed4f85fe2445cbbb1545. ("refused" in the file name is our internal label
for the version that takes no number from the Gata4 knockout files; the site has not refused this file.) The
external sources are the same as for our first Task 3 upload (X001, X004 and X005 in that email).

We built this file on 13 November, before any test score, as one of our two Task 3 files, and our procedure had it
go first and the other file last. On 18 November, before any Task 3 score, we uploaded the other file first, so this
one now goes second. The file itself is unchanged.

It consists of the released wild-type E8.75 cells, with two changes. First, the cell-type proportions are changed by
rules we fixed on 19 October: the cell states that change follow the broad phenotype background your answer of 15
October allows (cell states only, no genes), and the size of the change is computed from the released Mab21l2
knockout against its wild type. Second, each cell group's expression is adjusted by a correction estimated only from
released wild-type stage files: E8.75's own deviation from a trend through E7.5, E8.0, E8.25_late, E8.5, E9.5 and
E10.5. Three of these (Task 2's E7.5, E8.5 and E10.5) are official validation files we unlocked on 28 October, so
this upload uses official validation files beyond the Gata4 replicates that our earlier email listed for Task 3. Its
number of cells follows the same rule as our first Task 3 upload.

No number in the file is computed from the Gata4 knockout replicates. They entered only as practice tasks, scored
locally, together with the released Mab21l2 knockout. On those practice scores our rules of 19 October left a
knockout-derived per-gene factor out of this file, and on 13 November I, as team captain, decided to add the wild-
type correction described above; that was my decision and not the outcome of a rule fixed beforehand. On 18
November, after our first scores, we ran one more practice comparison of that correction's strength on the same
practice tasks. By the readings we had fixed before reading them, it changed nothing: the file is the one built on
13 November.

The method report will list the change of order.

Best regards,
Vijayavallabh J
Team captain, Mithril

## Data Protocol question on the E12.5 heart stack (drafted 2026-11-19; NOT reviewed, not in the captain's drafts; the captain decides)

Why: the heart's last upload keeps E10.5's full scope, and the released E10.5 stack is 23% heart core. Whether the
held-out E12.5 stack was dissected more closely is a test-stage premise no released task carries (`experiments/harness/TRIES.md`
row H26; `docs/PLAN.md` decision log 2026-11-19 (m) to (r)). Scope is a lever closed on 4 November, and the block's own
route to reopen one is "released data or a screened organiser answer". This is a first draft by the AI assistant,
written before the captain has chosen a route. It goes to three independent reviewers (accuracy, strategy, tone)
only if he chooses to ask, and into his drafts only at his word. It asks about the collection procedure and nothing
else; it states one fact about the released E10.5 file and nothing about our method.

**The reading of each answer, fixed before it is sent** (any answer is read through a screen and audited first;
numbers are withheld):
- The E12.5 hearts were dissected with most or all of the tissue outside the heart removed, or hold less surrounding
  tissue than the released E10.5 stack: the condition of row H26. Its conditional outcome applies once the captain
  has signed it: by the read of 2026-11-19, the full trim, which then needs his ruling on its size; or the half trim
  if he rules so.
- The two stacks were collected alike: the file in force stays.
- No answer, or an answer that does not say: the file in force stays. Nothing is assumed from silence.
- An answer that gives a number of the held-out file: the number is not read by the main session, and only the
  qualitative statement, if any, is used.
Answers: _none; not sent_.

**Email (draft).** Subject: Human Team (Mithril): a Data Protocol question on the E12.5 heart stack (Task 2)

To: virtual.embryo.moonshot@gmail.com

Hello organisers,

One short question about the Data Protocol page, for the Task 2 heart board.

The page describes the released E10.5 MERFISH hearts and the held-out E12.5 hearts in similar words: dissected, and sectioned whole. In the released E10.5 file, a large share of the cells belongs to tissue around the heart, not to the heart itself.

Could you say whether the E12.5 hearts were collected the same way, with neighbouring tissue left attached, or dissected more closely, as the page says of the E12.5 Multiome hearts ("most of the surrounding tissue removed")?

We are asking about the collection procedure only. If you would rather answer for every team at once, the Slack channel is fine for us.

Best regards,
Vijayavallabh J
Team captain, Mithril

**Slack thread (draft), if the captain chooses to ask there as well.** "A Data Protocol question for Task 2 (heart):
were the E12.5 MERFISH hearts dissected the same way as the released E10.5 ones, with neighbouring tissue left
attached, or more closely, as the page says of the E12.5 Multiome hearts? We are asking about the collection procedure only."

## Planned Slack DM (drafted 2026-11-06; sent by the captain the same day)

The state-of-the-art research recommends it (`experiments/research_2026_11_sota/REPORT.md` section 8; `docs/PLAN.md`
decision log 2026-11-06 (c)), and the captain approved it on 2026-11-06 (decision log 2026-11-06 (e)). He sent it
himself the same day, by Slack DM to Weize Xu, the private route at the top of this file. His words in the main session
were "i sent the slack dm". The main session had suggested the first paragraph only; whether the optional line went is
not recorded. No answer yet. He posts nothing in public. Any answer is recorded here and in `docs/DATA_POLICY.md` with its
date. An answer that arrives after T3's s2 is up changes nothing (decision log 2026-11-06 (b)). The second paragraph is
optional.

> Hi Weize, one short scoring question for Task 3, the same as in our email of 5 November; it affects every team the
> same way. When the severity metric is undefined (predicted slope <= 0, or R^2 <= 0.01), the current scorer returns
> log(1e-3), which equals the floor, so the live leaderboard gives skill 0.5, while the Evaluation page and the proposal
> call this case "undefined", which would map to skill 0. Which will the final evaluation of the test boards use: skill
> 0.5 or skill 0? A public answer in #virtual-embryo-challenge would be welcome, since it applies to all teams.
>
> (Optional) Also, are the embryo-series stage labels gestational (plug-dated) or morphological?
>
> Thanks, Vijayavallabh (captain, Mithril)

## Public Slack question (drafted 2026-11-11, rewritten the same day; the captain posts it)

The captain, 2026-11-11: "yes, draft the Slack question", on the AI assistant's offer, then "rewrite the draft after
thorough introspection". It asks in public what the 19 October follow-up (Q2), the fifth email and the Slack DM asked;
none has an answer. The 2026-11-06 choice to post nothing in public is superseded for this question only (`docs/PLAN.md`
decision log 2026-11-11 (e)). By the ruling of 2026-11-06 (b), an organiser's answer of zero switches T3's s2 to the
barred s1. Skill 0.5, a "yes" to the confirmation, or no answer before T3's s2 goes up keeps it. Any other answer (for
example "dropped") goes to the captain before T3's s2. An answer after T3's s2 is up changes nothing.

A three-lens panel (wf_b4e17573-4cb: accuracy, strategy, public tone) found the first draft worse than not posting:
- It put "undefined" in the mouths of the Evaluation page and the proposal. No record we hold quotes either document
  saying so; the label is our own paraphrase, first used in the 19 October follow-up. Our condensed Evaluation text gives
  this case "worst = log(1e-3)". The scorer's comment calls it "a finite floor, not credited and not dropped", and the
  scorer returns NaN ("undefined") only for board-level failures.
- Citing "undefined, which would suggest skill 0" and offering 0 as an equal option argued in public for the outcome
  that costs T3 9 to 12.5, and invited competitors to argue for it.
- "Wrong sign" could be read as a negative log-slope, which is scored normally; the gate is on beta.
- "We also asked by email" and the team signature pointed at unanswered mail and tied Mithril to the gated case.

The rewrite asks for a yes or no on the code's present behaviour. It names the scorer's own conditions, with the
no-response case beside them, and gives the value as the code's, not as a live observation. It says nothing about our
entry. Expected value against not posting: about +0.03 to +0.1 T3, with an answer by the screened check at about 0.2.
The first draft's was about -0.1 to -0.2.

Handling: post it once, as one top-level message with no @-mention, soon and in US Pacific working hours. The captain
does not argue in the thread and does not answer competitors, and replies only to clarify an organiser's unclear
answer. The screened Slack check reads the replies. A reply quoting test-board numbers, veckit 0.1.2 or PR #6 is
skipped. If nothing has come by about 19 November, one short reply may follow in the same thread (an @-mention is fine
there). Posted on 2026-11-12 by the captain (the captain's words: "i posted it"). The posted copy was not opened, so any
edit before posting is not recorded. Answer: **yes** (an organiser, 2026-11-16 19:53 UTC, in the thread; read through
the screen on 2026-11-19; quoted under "Screened reads of 2026-11-19" at the top of this file and in
`docs/DATA_POLICY.md`). T3's last upload stays as it is; no follow-up is posted. Text as drafted:

> Hi organisers, a quick Task 3 scoring check. Currently, veckit's `severity_slope` returns `log(1e-3)` when a
> prediction shows no response, or when `beta <= 0` or `r2 <= 0.01` (wrong sign, or almost no fit). That is the value
> `wt_identity` gets, so these cases score at the floor (skill 0.5). Can you confirm the final evaluation will score
> them the same way? A yes or no is plenty. An answer here before the final uploads keeps every team on the same rule.
> Thanks!

## Notes before the second uploads of the heart, the embryo and Task 3 (drafted 2026-11-19 after the captain's rulings of entry (cc); NOT reviewed by the three-lens review, not in the captain's drafts, nothing sent)

Drafted by an agent of the AI assistant from the captain's rulings of `docs/PLAN.md` decision log 2026-11-19 (cc),
checked by one reviewer whose corrections are applied here, and placed by a further agent (wf_1c3e47d4-cad and the
placing run). They are written for the files now named: the heart's full trim (row R10), the embryo's series twin
(row R11) and T3's barred s1 (row R12). The sections above (the heart and Task 3 notes drafted on 2026-11-18, the
embryo note drafted on 2026-11-19 for `T2embryo_s2x.h5ad`) are left as drafted: they are the record of what was
drafted then, and they name files that are no longer the last uploads. These three have not had the accuracy,
strategy and tone review the earlier notes had. They go into the captain's drafts only at his word, after that review,
and he sends them after the disclosure email. [TO CONFIRM: whether any earlier note or the disclosure email was sent;
each note below carries a marker for it.] The two hashes below were filled on 2026-11-19 from the build records, after each file's two builds agreed (decision log 2026-11-19 (ff)). Answers: none can exist yet.

### Note before the heart's second upload (the full trim)

Held: not to be finalised or sent before row H27's read is with the captain and he has said which heart file goes (docs/PLAN.md decision log 2026-11-19 (ee)).

**Subject:** Human Team (Mithril): second heart upload, file hash and sources

To: virtual.embryo.moonshot@gmail.com

Hello organisers,

As our email on sources and file hashes for the uploads of 18 November said we would, here are the file hash and sources of our second heart upload, before we make it. As it describes our method, please keep it private.

T2:heart:test_extrap, second upload. File T2heart_s2core_half.h5ad, SHA-256 04cde6361c04840ee47fff32eb9a20f4781bd649387cc540ed9fdea71e65fa35. It keeps only the heart core of the released E10.5 file: the cardiomyocytes and endocardium, and any other cells inside their in-plane spread. Every cell keeps its released expression and read counts. The cells are scaled evenly on all three axes to our forecast of the heart core's size. The uploaded file is a random half of the built file's cells, by the rule of our first heart upload (fixed seed 1000001, cell order kept). The external sources are the same as for our first heart upload (X001, X004 and X005 in that email), and so are the official validation files used (the heart's E10.5 and E8.5 files). [TO CONFIRM before sending: that this file uses the heart's E8.5 validation file, through the size forecast.]

This is not the file we had fixed in advance for this upload. That file trimmed the tissue around the heart to half its in-plane extent, sized it for that scope and changed its expression. [TO CONFIRM: whether the earlier heart note was sent. If it was: "It also replaces the file our earlier note described (T2heart_s2iso_half.h5ad, SHA-256 dc436d9225bb28e7997d2ec51e95254cd125e569d44c059596298d0a7b6af5c2), which kept the full section." If it was not, drop this sentence.]

Compared with the file fixed in advance, four things change. The scope is the core only, not half the section. The pool of cells the rows are drawn from is therefore smaller: the number of rows follows the same count rule as our earlier heart files, and the file is then halved by the rule of our first heart upload. The count rule does not depend on how many cells the trim keeps, so rows repeat only if the count exceeds the pool, each repeat with a small positional jitter. The size is an even scaling to the core forecast, where our procedure enlarged the two in-plane axes and held the third, at a figure for the half-trimmed scope. The expression step is dropped.

Scope, cell number and size are the quantities that the metadata exposure we reported on 4 November bears on. Since then we had kept all three fixed unless I decided otherwise. I, as team captain, made this choice after our first scores, which the site returned to us as to every team, and after that exposure. I had seen the figures that exposure showed, including the size of the held-out E12.5 heart file, so this choice cannot be shown to be free of them. The practice comparisons I had before me used released heart stages only. On them the trimmed file did worse than the file it replaces when the held-out heart is sampled like E10.5, and better when most of the surrounding tissue is removed. Those comparisons scaled every file to the practice heart's true size. No released stage shows how the held-out hearts were dissected.

No number in the file was computed from a score. Every number in it comes from a rule or a figure we had before this choice. The trim is the operation our fixed procedure applied at half strength, here at full strength. Our procedure had priced this full trim on 1 November and not chosen it. The count rule and the halving are our first heart upload's. The size is the figure our method computed and recorded on 28 October. The choice among them is mine. I made it on 19 November, between 19:25 and 19:31 UTC, in answer to a question put by the AI assistant we work with. It stated no preference and named no heart file. Before the question it had told me that an aim for the top three weighs each board's last file by its best case, with the first score kept, and the options were put with their practice gains and losses. [TO CONFIRM at the three-lens review: whether the next sentence goes to the organisers.] It had stepped back from every choice on the heart's scope and size after a lapse of its own in reasoning about the heart's size; nothing from that lapse was recorded or used. It is not the outcome of a rule fixed beforehand.

The method report will list each of these changes against the procedure we fixed in advance.

Best regards,
Vijayavallabh J
Team captain, Mithril

### Note before the embryo's second upload (the series twin)

**Subject:** Human Team (Mithril): second embryo upload, file hash and sources

To: virtual.embryo.moonshot@gmail.com

Hello organisers,

As our email on sources and file hashes for the uploads of 18 November promised, here are the hash and sources of our second embryo upload, before we make it. As it describes our method, please keep it private.

T2:embryo:test_interp, second upload. File T2embryo_s2xs.h5ad, SHA-256 bc50d012cd2c9481cd18c1f5bcf5b01093e22fd16ac65582a5c9f34321f13115. It has the cells of our first embryo upload, each placed at a position measured at E7.5. The expression of its E7.5 cells is re-rendered to the E8.0 file's per-gene detection and read depth; its E8.0 cells keep their expression as released. Its number of cells and its size are our first embryo upload's. Sources as for our first embryo upload (X001, X004 and X005), and the same official validation file (the embryo's E7.5 file).

This replaces the file we had fixed in advance for this upload, built on 13 November before any score, which used the same rendering on different positions. [TO CONFIRM: whether the earlier embryo note was sent. If it was: "It also replaces the file our earlier note named (T2embryo_s2x.h5ad, SHA-256 a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9). It has the same cells and positions and differs only in this rendering." If it was not, drop this sentence.]

I, as team captain, chose this file after our first scores, which the site returned to us as to every team. The practice comparisons behind the choice used released stages only. On them this rendering did better where the target was measured like the later neighbouring stage and worse where it was not; no released file shows how E7.75 was measured. Neither a rule we had fixed before our practice comparison nor a later check chose this file: the rule kept a file on our earlier positions, and the check kept the same cells and positions without this rendering. The choice is mine. It is not the outcome of a rule fixed beforehand. No number in the file was computed from a score. The method report will set out this change.

Best regards,
Vijayavallabh J
Team captain, Mithril

### Note before Task 3's second upload (the barred s1)

**Subject:** Human Team (Mithril): second Task 3 upload, file hash and sources

To: virtual.embryo.moonshot@gmail.com

Hello organisers,

As our email on sources and file hashes for the uploads of 18 November said we would, here are the file hash and sources of our second Task 3 upload, before we make it. As it describes our method, please keep it private.

T3:bcatenin, second upload. File barred_s1.h5ad, SHA-256 3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2. ("barred" in the file name is our internal label for the version our procedure kept ready for a set of organiser rulings listed below. It takes no number from the Gata4 knockout files. It is built from the decision file of our frozen rules' Gata4 run, whose rules fix this version as the Mab21l2 factor at full strength with no depth step. No ruling of yours has barred anything.) The external sources are the same as for our first Task 3 upload (X001, X004 and X005 in that email). [TO CONFIRM before sending: which official validation files, if any, this file reads.]

It consists of the released wild-type E8.75 cells, with two changes. First, the cell-type proportions are changed by rules we fixed on 19 October. The cell states that change follow the broad phenotype background your answer of 15 October allows (cell states only, no genes). The size of the change is computed from the released Mab21l2 knockout against its wild type. Second, each cell's expression is moved by a per-gene factor computed from the released Mab21l2 knockout: its E9.5 knockout cells against the wild-type E8.75 cells. Its proportions are those of the file it replaces, and its number of cells follows the same rule as our first Task 3 upload. It reads no Gata4 knockout file.

We built this file on 13 November, before any test score. It is on the list of file hashes we fixed that day for our uploads. Our procedure sent it to this slot only if, before this upload, a ruling of yours restricted or barred the Gata4 files, or barred our use of the validation answers, or you said that a floored severity scores zero; none is on our record. It replaces the file our procedure had for this upload (refused_s1.h5ad, SHA-256 0394575a11fcd4f915e0daf4972e5a783b6b89efa065ed4f85fe2445cbbb1545), which carried a wild-type-only correction in place of the Mab21l2 factor. [TO CONFIRM: whether the earlier Task 3 note, which named refused_s1.h5ad, was sent. If it was: "Our earlier note named that file; this one replaces it." If it was not, drop this sentence.]

I, as team captain, chose this file after our first scores, which the site returned to us as to every team. The Gata4 replicates entered only as practice tasks, scored locally. On 18 November we compared the two files on them. By the reading we had fixed before that comparison, and by a later check, the other file stayed. The choice is mine. It is not the outcome of a rule fixed beforehand. No number in the file was computed from a score.

The method report will list this change.

Best regards,
Vijayavallabh J
Team captain, Mithril

### Lines of the disclosure email draft that the rulings of 2026-11-19 make inaccurate

[TO CONFIRM: whether the disclosure email was sent; if it was, these corrections go into the note that precedes each upload instead]

The email is quoted as it stands in the section "Disclosure email for the uploads of 18 November" above (the text of
18 November with changes (1) to (11) of 2026-11-19). That section is not edited.

1. The paragraph of upload 4. Quote: "This is the one change to a cell number since our note of 4 November, and the
   method report will list it with the other later decisions that bear on the quantities that note covers."
   Why it now misleads: the heart's last file keeps the count rule and the halving, so the number of rows is not
   changed again; but it changes the scope and the size, and with the scope the pool its cells are drawn from, and
   the sentence reads as if nothing else that note covers had changed.
   Replacement: "This is the one change to a cell number since our note of 4 November. Our second heart upload uses
   the same count rule and halving, but changes its scope and size, and with the scope the pool its cells are drawn
   from, as the short note before that upload sets out. The method report will list each of these with the other
   later decisions that bear on the quantities that note covers."
2. The correction paragraph. Quote: "One correction to our note of 4 November. That note said our cell counts and
   sizes follow rules committed before the data release. For Task 3, the rule for the number of cells (use the
   board's maximum) was fixed on 1 November 2026: after the release, though still before what we saw on 4 November.
   The method report will say the same."
   Why it now misleads: it is not false sentence by sentence, but read beside the heart's full trim it suggests that
   T3's count rule is the only departure from the note of 4 November. The heart's last file has a scope and a size
   chosen after the first scores and after the exposure.
   Replacement: keep the paragraph and add, before "The method report will say the same.": "Our second heart upload
   also departs from that note: its scope and size, and with the scope the pool its cells are drawn from, were chosen
   by me after our first scores, as the short note before that upload sets out."
3. No other sentence of the email is made inaccurate by these three rulings:
   - "Before each of our second uploads we will send its board, file hash and sources in a short note." still holds.
   - The validation paragraph ("please tell us before our second uploads.") still holds. The reason the record gave
     for change (10), "every board's second upload uses an official validation file", may no longer hold for T3's
     barred file (which validation files it reads is still to confirm; see the marker in the Task 3 note above), but
     that reason is not in the email.
   - The register sentence ("both Task 2 boards and Task 3 use X001, X004 and X005") holds for all three new files.
   - "no gene-level knockout finding from outside the released data" holds: the barred file's factor comes from the
     released Mab21l2 knockout.

## Note before Task 1's second upload (drafted 2026-11-19 after the captain's ruling of entry (ii); NOT reviewed by the three-lens review, not in the captain's drafts, nothing sent)

Drafted by an agent of the AI assistant from the captain's ruling of `docs/PLAN.md` decision log 2026-11-19 (ii), row
R13 of `experiments/harness/TRIES.md` and the amendment of `docs/S2_RULES.md` for the recomputed-groups file. It is
written for that file (`T1_s2rc.h5ad`). The section "Notes before Task 1's and the embryo's second uploads" above is
left as drafted: its Task 1 note names `T1_s2none.h5ad`, which is no longer the last upload. This note has not had the
accuracy, strategy and tone review the earlier notes had. It goes into the captain's drafts only at his word, after
that review, and he sends it after the disclosure email. [HASH] is filled from the t1rc build record ([MAIN SESSION:
the build record]); the file is built only after row T63's read ([MAIN SESSION: row T63's outcome]). [TO CONFIRM:
whether the captain sent the disclosure of the X002 labels or the sources note.] Answers: none can exist yet.

**Subject:** Human Team (Mithril): second Task 1 upload, file hash and sources

To: virtual.embryo.moonshot@gmail.com

Hello organisers,

As our email on sources and file hashes for the uploads of 18 November said we would, here are the file hash and sources of our second Task 1 upload, before we make it. As it describes our method, please keep it private.

T1:test, second upload. File T1_s2rc.h5ad, SHA-256 [HASH]. It is our first Task 1 upload's method with one change: the step that uses the external heart atlas (Feng et al. 2022, GEO GSE193346) is computed on cell groups we recomputed ourselves from the atlas's permitted cells, in place of the authors' cell-type labels. The cells, the cell number and the composition are those of our first Task 1 upload.

[TO CONFIRM: whether the earlier Task 1 note, which named T1_s2none.h5ad, was sent. If it was: "It replaces the file our earlier note described (T1_s2none.h5ad, SHA-256 55eaee547ff18d7477da5b730c34e7de9ab74f85cc12fd01a9eff6450367611c), which used no external data." If it was not, drop this sentence.]

What entered the recomputation: the raw counts of the atlas's wild-type CD1 cells at E9.5 and E14.5, after every cell from E10.5 to E13.5 was removed on the per-cell stage column; the gene names; and the cell name, stage, zone and strain columns. Zone is used only to split each group's cells into two halves for checks (the direction itself pools every zone), and strain only to confirm that each cell is CD1; neither is a cell-type label. Each cell was scored on fixed sets of marker genes against control genes of similar expression, and assigned to a group by a fixed threshold, margin and endothelial split. The marker genes of the four groups we use (cardiomyocytes, endocardium, endothelium, epicardium) come from the released Task 1 E9.5 file's cell-type means. Those of the other cell types, which keep non-chamber cells out, come from eight published papers on adult mouse tissue. No clustering, embedding or classifier was used, and the authors' cell-type labels and every other annotation column were not read. The settings were fixed, and a dry run on the released E9.5 file set one threshold by a rule fixed beforehand, before the atlas's cells were read.

One thing cannot be undone. The authors' quality filter, which they applied to all stages together, decided which permitted cells exist in the files we have. The zone labels also come from their sample processing. Earlier choices in our method (this dataset as the direction source, its E9.5 to E14.5 bracket, and a list of genes held flat) were made with analyses that used the authors' labels; our labels note and the method report list them.

The practice comparisons behind this file were scored locally against released Task 1 stages only. On them the step on the recomputed groups did better than the same method without external data on all five tasks we ran. Those tasks look one day ahead with half the step the test needs, so they do not show what the file will score.

Sources: X001 (local scoring only), X005 (pinned libraries), X002 (Feng et al., on the recomputed groups), X010 (the published marker genes), and X007 and X009, the gene lists held flat in the atlas's direction as in our first Task 1 upload. X008 informed that list's premise. X003 and X006 never entered a model. Official validation file: E10.5, as before. [TO CONFIRM before sending: the captain confirms this list of register ids.]

I, as team captain, chose this file after our first scores, which the site returned to us as to every team. I followed the advice the AI assistant we work with had fixed before its practice comparison; a later check let that advice stand and named the risk that you may not accept a second file that uses this atlas. No number in the file was computed from a score. The method report will describe the recomputation, its inputs and the filtering, as your reply of 16 November asks.

Best regards,
Vijayavallabh J
Team captain, Mithril
