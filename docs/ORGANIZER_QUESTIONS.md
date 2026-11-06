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
