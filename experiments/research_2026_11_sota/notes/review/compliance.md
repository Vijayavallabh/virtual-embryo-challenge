# Compliance audit: experiments/research_2026_11_sota/REPORT.md and its notes/ folder

Reviewer: compliance auditor (AI agent), 2026-11-06. Read-only; this file is the only write.
Standards: AGENTS.md; docs/DATA_POLICY.md (windows, Scores, do-not-use list, the post-exposure block 2026-11-04
(a)-(g), change notes to 2026-11-05 (k)); notes/brief.md section 6 (classes A, B, W, L, Y, Q); docs/PLAN.md section 8
and decision log 2026-11-05 (c)-(g), 2026-11-06 (a)-(b); docs/S2_RULES.md; docs/ORGANIZER_QUESTIONS.md.

## Verdict

Not ready to commit as is. REPORT.md itself carries no hidden-target value, test-board limit, per-board leaderboard
figure of any team, gene-level beta-catenin/Gata4 knockout finding, in-window result, Yukon material or quarantined
artefact, and no sentence in it ranks, chooses or gates a method by 209. But its compliance front matter points to a
change note that does not exist and misstates the class-B record; three agents opened memory files and the report does
not say so; task-e.md keeps a residue of a class-A removal; four design files use 209 as a planning threshold (their
own critics flagged it); and four recommendations drop conditions that the rules or the critics attached (organiser
channel, per-upload go-ahead and T3 screened check, the provenance of the captain's embryo prior, the composition and
A15 conditions). Fix list below; each fix is text-level and none needs a run.

## 1. The four checks

**Check 1, barred material in REPORT.md and notes/.**
- REPORT.md: no issue for classes A, B, L, Y, Q and W values. The organisers' published baseline rows (48.58, 47.76,
  55.22, 42.48, 42.50) are the allowed class. The Mesp1-Cre phenotype paragraph (section 4) is phenotype-level, which
  the 2026-10-15 ruling allows. One wording can be read as a test-stack property (issue 11).
- notes/task-e.md: one residue of a class-A removal (issue 4); its compliance log describes barred items beyond "class
  and source" (issue 4). The class-L removal in Deep Read [25] is clean.
- notes/task-b.md: no barred content found. The L1/L2 E12.5 illustrations are derived from legal anchors and released
  E10.5 (the auditor's call holds). The legality conditions its audit attached are dropped by the report (issue 10).
- notes/task-a.md, task-a2.md, task-c.md, task-d.md, task-f.md: no barred content found. task-c's phenotype items keep
  only phenotype-level statements; its class-B meetings are in the compliance log only (issue 2 is about the report).
- notes/design/*.md: no per-board leaderboard figure, hidden value or bytes-per-cell figure. 209 misuse (issue 5).
- notes/registry.md: entry [134] (dxyang42 blog) is "Approved" although its only role in task e is as the place a
  possible class-A item was met (issue 12).
- notes/brief.md: the overall-total line lacks the calibration label (issue 14).

**Check 2, uses of 209.** REPORT.md: no issue. Every use (lines 11, 27, 67, 111, 144-158, 538) reports it; section 2
labels it "for information only ... It chooses nothing here", and the verdicts in section 8 rest on corrected expected
gains, not on 209. The design notes do not pass: t1.md section 4, t3.md section 5 (and line 11), heart.md section 6
(and line 42), embryo.md section 6 (and line 30), cross.md (title, lines 30-36, the P(>= 209) column of section 2, and
section 4) measure candidates against 209; task-a2.md line 144 does so in passing (issue 5).

**Check 3, recommendations against the rules.** Issues 6 (organiser channel), 7 (embryo prior provenance), 8 (A15 and
composition conditions, adoption-rule override not stated), 9 (go-ahead and T3 screened check), 10 (Leonard legality),
17 (T3 C1's full set of reversals). No issue for: "scores choose, never compute" (no recommendation reads a score);
"do not front-load" (consistent with S2_RULES); recommendation 5 (T3 C1 not recommended); recommendation 6 (watch
through a screening agent).

**Check 4, process notes for docs/DATA_POLICY.md.** The note the report cites, 2026-11-06 (a), does not exist (issue
1). Section 3 below lists every item the notes record.

## 2. Issues (most severe first)

### Issue 1 (high). The cited change note does not exist
- Location: REPORT.md front matter, line 19: "Process notes are in `docs/DATA_POLICY.md` 2026-11-06 (a)."
- What is wrong: docs/DATA_POLICY.md ends at 2026-11-05 (k); there is no 2026-11-06 entry. DATA_POLICY's header
  requires a dated note for every change, and the notes record slips that the policy has recorded every time
  before (memory files opened, harness-saved outputs, git without GIT_OPTIONAL_LOCKS=0, practice-truth displays).
- Fix: write docs/DATA_POLICY.md 2026-11-06 (a) in the same commit as the report, from section 3 below, headed
  "process notes from the state-of-the-art research workflow of 2026-11-06 (notes/ in
  experiments/research_2026_11_sota/); no test-board limit, anchor or truth count was seen". Do not commit the report
  while line 19 points to a missing entry.

### Issue 2 (high). The report says class B was not met; task c met it twice
- Location: REPORT.md front matter, line 19: "No gene-level knockout finding (class B) was met." Also line 13: "No
  agent saw a test target or a test-board limit".
- What is wrong: notes/task-c.md compliance log, lines 274-275: "Class B met in one WebSearch result summary (query on
  canonical Wnt signalling, second heart field and RV/OFT phenotype; ...)" and "Class B met at
  https://pubmed.ncbi.nlm.nih.gov/18000065/ (abstract text quoted in a search summary)". Nothing was recorded, but
  class B was met, on the held-out gene. The screen counts (B 0) count removals, not meetings. Class L (task-e.md
  276-277) and class Y (task-e.md 278-280) were also met in listings without being recorded, and task e's research
  agent met organiser statements that may describe hidden targets (the six class-A removals).
- Fix: replace line 19's first sentence with "Class B was met twice, in search-result summaries read by task c's
  research agent while searching the Mesp1-Cre phenotype; nothing was recorded or used, and the screen found none in
  the note. Class L and Y items were met in search listings by task e (not opened, not recorded)." Replace line 13's
  first clause with "No agent opened a test file, a board file or the Data page, and none saw a test-board limit;
  task e's agent met organiser statements that may describe hidden targets (class A, removed below)".

### Issue 3 (high). Three agents opened memory files; the report is silent and the rule is not holding
- Location: REPORT.md front matter (lines 13-19) lists removals only; the slips are in notes/task-a.md line 232
  (screener: one memory file), notes/task-e.md lines 304-306 (screener: one memory file, "named in its loaded memory
  index"), notes/task-b.md line 242 (auditor: two memory files, "the project's challenge note and its
  compliance-handling note").
- What is wrong: AGENTS.md (Quarantined material) and DATA_POLICY 2026-11-05 (d)(iv) and (k) bar agents from memory
  files. This is the second day running (2026-11-05 (k) recorded the same slip). The harness injects the memory index,
  which names the files, into every agent (reported in at least 13 agents' logs, section 3), and the task-e screener
  says it found the file through that index. task-b's auditor also says the files "held overall leaderboard totals",
  while 2026-11-05 (k) says the totals were taken out of the project memory file.
- Fix: add to line 19: "Three agents (task a's and task e's screeners, task b's auditor) opened memory files under
  ~/.claude-account against AGENTS.md; nothing from them was copied (DATA_POLICY 2026-11-06 (a))." In the change note,
  record all three and a systemic fix for the captain's choice: a first line in every workflow brief and prompt, before
  anything else, "do not open any file the injected memory index names", or turning off auto-memory injection for
  workflow subagents (a settings change, his call). The main session reconciles the "overall totals" remark against
  (k) from its own memory, without asking an agent to open it.

### Issue 4 (medium). task-e.md keeps a residue of a class-A removal

[Body removed in the main session, 2026-11-06: it described the topic of a removed class-A statement. The fix was applied: notes/task-e.md lines 71, 273 and 281-284 now carry neutral markers only (docs/DATA_POLICY.md 2026-11-06 (a)).]

### Issue 5 (medium). Design notes use 209 as a planning threshold and are about to be committed
- Location: notes/design/t1.md section 4 (lines 429-453: "the gains T1 would need", "about a fifth of the 17.5 points
  between reading A's mean and 209"); t3.md line 11 and section 5 (lines 305-319: "The candidates here move that
  probability by well under 0.01"); heart.md line 42 and section 6 (lines 486-508, which ends in a try-budget
  recommendation); embryo.md line 30 and section 6 (lines 355-371); cross.md title, lines 30-36, the P(>= 209) column
  of the ranking table (lines 401-412) and section 4 (lines 446-475: "The remaining +16 or so ... has to come from");
  task-a2.md line 144 ("moves the estimate materially toward 209"). Source of the per-board figures:
  experiments/gate/path_to_target.py, `whatif_p_target_after_sure_shift` and `d_for_p_half` (used at t1.md 434-439 and
  cross.md 56).
- What is wrong: DATA_POLICY 2026-11-05 and docs/PLAN.md section 8: 209 "gates uploads only; it is never a threshold or
  evidence for a try, an adoption or a method choice". The T1 critic (t1_critique.md 196-200: "Compliance. It uses the
  gate figure as a planning threshold") and the T3 critic (t3_critique.md 51-53) flagged it; the cross critic asked to
  take the column out (cross_critique.md 233-234). Precedent 2026-11-04 (g)(v): misused design passages were struck and
  the design texts not committed. REPORT.md cites notes/design/, so they will be committed.
- Fix, before commit: (a) put one line at the head of each passage listed: "> Struck 2026-11-06: this measures
  candidates against the captain's 209, which gates uploads only and is never a threshold or evidence for a try or
  method choice (docs/DATA_POLICY.md 2026-11-05; docs/PLAN.md section 8). Information only; nothing here chooses.";
  (b) delete cross.md's P(>= 209) column and report P(>= 209) once in its section 4, labelled; (c) drop the per-board
  209 what-ifs (`whatif_p_target_after_sure_shift`, `d_for_p_half`) from path_to_target.py and its JSON before the
  commit (REPORT.md uses neither; keep `p_at_least.target_209` and the overall-gain table, labelled); (d) record the
  strike in 2026-11-06 (a). If all five designs answered a prompt that asked for a "Path to 209" section, record that
  too (cf. 2026-11-05 (b)(vii)).

### Issue 6 (medium). Recommendation 2 adds an organiser recipient outside the documented channels
- Location: REPORT.md section 8, recommendation 2, line 495: "to the organiser contact `docs/ORGANIZER_QUESTIONS.md`
  names for private questions; optionally also the proposal's lead for metrics [108]." Same in the executive summary,
  item 2, line 52 ("by Slack DM to an organiser").
- What is wrong: the documented routes are the organisers' email or "a Slack DM to Weize Xu"
  (docs/ORGANIZER_QUESTIONS.md line 3) and, for rule questions, the email or #virtual-embryo-challenge (AGENTS.md,
  Escalation). The proposal's metrics lead (identified only in task-e's deep read of [2]) is not a documented contact.
  The report also drops the T3 critic's conditions (t3_critique.md 153-158): ask for a public answer, post nothing in
  public unless the captain overrides the ask-privately rule, record any answer with its date in DATA_POLICY.
- Fix: line 495 becomes "to Weize Xu, the private route `docs/ORGANIZER_QUESTIONS.md` names, word for word as the fifth
  email, asking, as it did, for a public answer in #virtual-embryo-challenge. A second recipient only if the captain
  first adds that contact to `docs/ORGANIZER_QUESTIONS.md`. Any answer is recorded with its date in
  `docs/DATA_POLICY.md`; one that arrives after T3's s2 is up changes nothing (ruling 2026-11-06 (b))." Line 52: "by
  Slack DM to Weize Xu, the private route ...".

### Issue 7 (medium). The embryo s2 ruling rests on the captain's prior about a target whose cell limits he has seen
- Location: REPORT.md recommendation 4, line 508: "Embryo s2 at fraction 1/4: worth it only if he himself puts P(E7.75
  is E7.5-like) at about 0.2 or more." Executive summary lines 57-58.
- What is wrong: the deciding input is the captain's own prior on a property of the hidden E7.75 target. On
  2026-11-04 he displayed the board index and the Data page, which show every test board's cell limits, the embryo
  test board's included, and those limits "can encode a held-out stage's cell count" (DATA_POLICY 2026-11-04 (a), (b)).
  The cross critic flagged the provenance of the "E7.5-like" risk, logged on the day of that exposure, and asked that
  its basis be the released-stage diagnostic of 2026-11-03 alone (cross_critique.md 107-110). The report drops the depth
  reason but not the provenance issue, and does not name S2_RULES's change procedure.
- Fix: append to line 508's sub-bullets: "His prior rests only on committed released-stage evidence (the committed
  labels' staging spread, which the embryo critic frames as sigma >= 0.15 day; the 2026-11-03 depth diagnostic; the
  E7.5 and E7.25 backtests), never on anything seen on the Data page or in the board index on 2026-11-04; the ruling
  says so, and goes in its own commit with his words and the new sha256 (`docs/S2_RULES.md`, Changes after the
  freeze)." Add "(a prior from released-stage evidence only)" to line 58.

### Issue 8 (medium). Recommendation 4's A15 and cycling-quota bullets drop required conditions
- Location: REPORT.md recommendation 4, lines 511-514 ("A15 read as not covering T1-only gene sets ... It needs a
  register row per source and the gene classes fixed before any build." and "A cycling-state quota in T1 categories:
  about +0.05, needs A15's reading or the released labels."). Executive summary lines 59-60.
- What is wrong: (a) the cycling quota is anchored on Leonard's E14.5 WT nuclei (t1.md 303-306), so GSE272154 needs
  the Q7 Theiler-stage check and an EXTERNAL_DATA.md row first (task-b.md line 80; AGENTS.md: every dataset logged
  before it enters a model); (b) it is a composition change: rule and every constant in its TRIES row before any run
  and before anything is computed on E10.5's within-category shares, with the library-size report (DATA_POLICY
  2026-11-04 (d), (f)(1), (f)(3); t1_critique.md 157-164); (c) it cannot pass rule 5, so it enters only by ruling, as R1
  did (t1_critique.md 149-152), which overrides AGENTS.md rule 6 and should be said; (d) for the A15 bullet, the cross
  critic's condition "no class chosen by reading E10.5's DE genes (a counted truth)" and PLAN 4.4-style provenance
  (cross_critique.md 179-181) are missing.
- Fix: line 511-513 becomes "A15 read as not covering T1-only gene sets, for a Feng-direction hygiene step: about
  +0.10. Conditions: an EXTERNAL_DATA.md row per gene-set source with docs/PLAN.md 4.4-style provenance; gene classes
  fixed in the ruling row before any build, none chosen by reading E10.5's DE genes; R1 untouched; T1's s2 back to
  him. No counted evidence exists, so it enters by his ruling, not the adoption rule." Line 514 becomes "A cycling-state
  quota in T1 categories, anchored on Leonard's E14.5 WT nuclei: about +0.05. GSE272154 first passes the Q7 staging
  check and gets its EXTERNAL_DATA.md row; the rule and every constant go in a TRIES row before any run and before
  anything is computed on E10.5's within-category shares, with the (f)(3) library-size report; cycling calls come from
  the released labels or a gene list under his A15 reading; below rule 5's bar, so by his ruling only; T1's s2 goes
  back to him."

### Issue 9 (medium). The dated runbook omits the per-upload go-ahead and the screened check before T3's first upload
- Location: REPORT.md recommendation 3, line 500: "27 Nov 05:30 IST: rebuild at the freeze commit, checks and the first
  round;" (executive summary item 3, lines 54-55).
- What is wrong: the first round includes T3's s2, which the 2026-11-06 (b) ruling chooses "at the screened check of
  mail, Slack and the Rules page that precedes T3's first upload", run no earlier than 2026-11-26; S2_RULES (T3, Order)
  says every T3 upload follows such a check; and each official upload is run by the captain himself on his go-ahead at
  the time (docs/PLAN.md section 8; S2_RULES, Go-ahead). The cross designer's version had "with T3's s2 per the screened
  answer check" (cross.md 432) and the critic's "pre-filled go-ahead sheets" (cross_critique.md 87 and 240); the report
  lost both.
- Fix: line 500 becomes "27 Nov 05:30 IST: rebuild at the freeze commit and checks; then the first round, each upload
  run by the captain on his go-ahead at the time, with T3's s2 chosen at the screened check of mail, Slack and the Rules
  page that the ruling of 2026-11-06 (b) requires (run no earlier than 26 Nov);". Add after line 503: "Every T3 upload
  follows a screened check (`docs/S2_RULES.md`, T3)."

### Issue 10 (medium). "New legal resources" drops the audit's legality conditions on Leonard and mdCLA
- Location: REPORT.md section 3, lines 205-210 ("New legal resources ... Leonard et al. (2025) ... about 29%
  cardiomyocytes and 26% fibroblasts; GSE272154 [50]"; "mdCLA ... and 18 points at or before E9.5"); key finding 2,
  line 532 ("New legal anchors exist (Leonard 2025, mdCLA 2026)").
- What is wrong: task-b's audit left one window question on Leonard: "check GSE272154's sample metadata for a
  Theiler-stage code before registration, since an E14.5 sample coded Theiler stage 22 stays out under Q7 (X006
  precedent)" (task-b.md line 80), and said the pooled 29% "may motivate the check but never set the prior" (line
  108). D3 (the CM-share prior) is open, and 2026-11-05 (e) redacted an unvetted CM share from every file a D3 reader
  could see. mdCLA's early points include E7.5 and E7.75 (embryo window) and E8.5 (heart interpolation window), usable
  for T1 alone or after removal by code for any shared-core step (task-b.md line 55).
- Fix: line 206-207 ends "...about 29% cardiomyocytes and 26% fibroblasts, pooled with three mutant (PMIS) libraries:
  motivation only, never an input to a share prior. Legal by its E-day label, pending the Q7 Theiler-stage check of its
  sample metadata and an EXTERNAL_DATA.md row [50]." Line 210 ends "...18 points at or before E9.5 (E7.5, E7.75 and
  E8.5 among them sit in the embryo and heart-interpolation windows: T1-only use, or removed by our code for any
  shared-core step) [33][34]." Line 532: "New legal anchors exist, each pending its register row and staging checks
  (Leonard 2025, mdCLA 2026)".

### Issue 11 (low). A released-stage fact reads as a property of the test stack
- Location: REPORT.md section 5, line 320: "But our scored heart cloud is not the organ. At E9.5, 63% of cells are
  non-cardiac shell; every stage has 30-40 sections;".
- What is wrong: the source says "Every file has 30-40 sections" of the released files (experiments/coords/REPORT.md
  line 49) and lists the test stacks' section counts as an open question (line 170); task d calls the test stack "class
  A territory". After "our scored heart cloud", "every stage" can be read as including E12.5.
- Fix: "every released file has 30-40 sections (the test stacks' counts are unknown and not sought)".

### Issue 12 (low). Sources with barred content are not marked; one citation points to the wrong entry
- Location: REPORT.md section 12 (lines 560-654) and line 327 ("[13 in task d, cited for method only]"); notes/
  registry.md entry [134].
- What is wrong: the notes' compliance logs record class W on the pages of [6], [8], [10], [38], [39], [40], [79],
  [80], [82], [85], [98], [105], class B in the text of [78] (the registry says "read with a no-genes prompt"), and
  possible class A on [108], [109], [111]. The 2026-11-05 (k) precedent marked its method-only source in the report
  ("[87] TIMO (method point only; in-window results not used)"). Global [13] in the list is Chreode, which the body
  never cites; the line-327 citation is task d's Raiola et al., which the registry skipped. Registry [134] is approved
  although task e met a possible class-A item there and cites it nowhere else.
- Fix: add under the section 12 header: "Marks: (W) the page also reports in-window results, used for method or
  legal-stage points only; (B) the text holds gene-level knockout findings, phenotype only, read with a no-genes
  prompt; (A) the page holds statements that may describe hidden targets, read only through a screening prompt." Tag
  [6][8][10][38][39][40][79][80][82][85][98][105] (W), [78] (B, W), [108][109][111] (A). Replace line 327's bracket
  with "(Raiola et al., eLife 2026, task d's source 13; method only)" and drop [13] from the list. In the registry,
  move [134] out of "Approved" or tag it "(A) do not open".

### Issue 13 (low). Notes record full session-store paths and file names
- Location: notes/task-f.md lines 174-175 (full absolute paths with the session id and file names), task-a.md line
  222 (session path), task-b.md line 158 (file names bet24qofk.txt, toolu_01NQYGHtZf688ygAAotKaDx8.txt).
- What is wrong: committed notes then map class-Q files that agents must never open; the other notes write
  "~/.claude-account/.../tool-results/".
- Fix: replace each with "~/.claude-account/.../tool-results/ (a harness-saved output of this agent's own call; not
  opened)".

### Issue 14 (low). The brief carries the overall figure without the calibration label
- Location: notes/brief.md section 3, line 58: "A screened read of overall totals puts a Human top-3 place near
  202-210. Our gate: 209."
- What is wrong: DATA_POLICY 2026-11-05 (b)(vii): prompts carry "at most the overall figure, labelled as calibration".
  The brief went to every research and design agent.
- Fix: append "(overall totals, calibration only: never a threshold or evidence for a try, adoption or method choice;
  the gate figure gates uploads only)". Use the same label in future briefs.

### Issue 15 (low, accuracy). The organisers' rows are not all at or below the floor
- Location: REPORT.md executive summary, line 44: "**The organisers' own models** sit at or below the floor on the
  practice boards."
- What is wrong: task-e's spot-check 1 corrected this: pseudobulk_shift scores 55.22 on embryo interpolation, above 50
  (task-e.md line 44); section 6 of the report lists 55.22 itself.
- Fix: "sit at or below the floor on the T1, heart-extrapolation and Gata4 practice boards (pseudobulk_shift reaches
  55.22 on embryo interpolation)."

### Issue 16 (low, operational). Recommendation 1 covers the local tree only
- Location: REPORT.md executive summary item 1 (lines 50-51) and recommendation 1 (lines 491-493); section 2, line 114
  ("`results/path_to_target_2026-11-06.json`").
- What is wrong: cross.md lines 506-511 found path_to_target.py untracked on the DGX too, plus two older harmonise CSVs
  there. Uploads go `--from-dgx` (S2_RULES, Builds) and need the DGX tree clean (guard "dgx provenance",
  src/vec/submit.py); a pull of the commit that adds the file may be refused while the untracked copy exists. The JSON
  sits at experiments/gate/results/, not the top-level results/.
- Fix: add to recommendation 1: "On the DGX, compare the untracked copy with the committed file, remove it, pull, and
  confirm the tree passes the clean-tree rule (the two harmonise CSVs included)." Line 114: use
  `experiments/gate/results/path_to_target_2026-11-06.json`.

### Issue 17 (low). T3 C1's reversals are listed incompletely
- Location: REPORT.md section 8, T3 bullet, lines 446-447 ("C1 is close to the 2026-11-04 exposures ... It would also
  reverse the 2026-11-01 (b) decision.").
- What is wrong: the T3 critic adds that C1 "reopens the closed T3 depth line (TRIES rule 6; the QC20 entry) and departs
  from D4 (the frozen T3 rules)" (t3_critique.md 58). The captain should see all three reversals.
- Fix: append "It would also reopen the closed T3 depth lever (TRIES rule 6, the QC20 entry) and depart from the frozen
  T3 rules (D4)."

### Issue 18 (low). A tracked repository file holds a possible class-A/W item outside the quarantine list
- Location: notes/task-e.md Sources [2] (line 9) and compliance log line 283: docs/reference/neurips-2026-competition-
  proposal.pdf, p.7 (an organiser evaluation on held-out knockout data) and its Data section, skipped on purpose.
- What is wrong: the file is tracked (`git ls-files docs/reference`), and AGENTS.md's quarantine paragraph does not
  name it, so the next agent that reads the proposal meets the item unwarned.
- Fix: record it in 2026-11-06 (a) and add to AGENTS.md (Quarantined material): "Read
  docs/reference/neurips-2026-competition-proposal.pdf skipping its Data section and its page-7 evaluation on
  held-out knockout data." (An AGENTS.md change is the main session's, under the delegation.)

## 3. Process notes for docs/DATA_POLICY.md 2026-11-06 (a)

Memory and harness files (class Q):
- Memory files opened, against AGENTS.md and 2026-11-05 (d)(iv), (k): task a's screener (one, the compliance-handling
  note; task-a.md 232); task e's screener (one, found through the injected index; task-e.md 304-306); task b's auditor
  (two, the project's challenge note and the compliance-handling note, before reading the brief; it reports overall
  leaderboard totals in them; task-b.md 242). Nothing copied.
- Memory index injected by the harness, not opened by a tool call: task-a auditor (240), task-a2 screener (241) and
  auditor (315), task-b research agent (158) and auditor (242), task-f research agent (176) and auditor (223), t1.md
  (476), t1_critique.md (239), t3.md (345), t3_critique.md (259), cross.md (496), cross_critique.md (278). It also
  reached this audit agent; not opened.
- Harness-saved tool outputs under ~/.claude-account, none opened: task-a (222), task-a2 (228, two), task-b (158, two),
  task-c (280, two), task-d research agent (234, two) and auditor (279, one), task-e (286, two), task-f (174-175, two),
  embryo_critique.md (293, one), t3.md (343, one).

Git:
- t3_critique.md (254): one `git log --oneline -5` without GIT_OPTIONAL_LOCKS=0.
- Read-only git looks without stating GIT_OPTIONAL_LOCKS=0, to confirm: t1.md (479-482), t1_critique.md (243-248),
  embryo.md (397-399), embryo_critique.md (299-302).
- cross.md (506-511) and cross_critique.md (287-293) used GIT_OPTIONAL_LOCKS=0; cross.md also read the DGX repository's
  git state over SSH (read-only).
- The design prompts relayed a user request to sync the DGX and the local machine and check for unpushed commits; every
  design agent left it to the main session.

Displays of practice-truth and validation files (no test board):
- task-a auditor: a grep displayed validation-board anchor lines of
  experiments/calibration/results/calibration_summary.txt (task-a.md 240).
- embryo designer: a committed practice-task file displayed the E7.5 validation board's spec and anchors and that task's
  truth DE counts (embryo.md 393-394).
- cross designer: payoff.txt's footer printed the re-sampled sizes of the T1 grid's practice truths (cross.md 492-494);
  cross critic: a grep printed the first 60 characters of payoff.txt line 20, a practice truth's mix and the start of
  its size (cross_critique.md 274-276).
- heart designer and critic: the first lines of experiments/labels/results/E10.5__vs__E9.5.csv (released and unlocked
  label counts; heart.md 523-524, heart_critique.md 341-342).
- embryo critic: the header and E7.5 row of experiments/coords/results/p3/file_summary.csv (embryo_critique.md 289).
- task-c auditor: the D10 and W1 blocks of the Gata4 and Mab21l2 decision files (task-c.md 364).
- t1 designer and critic: in-memory re-weighting of the committed per-seed T1 grid (score and skill means printed;
  t1.md 470, t1_critique.md 26 and 235).
- Kept class shown: one overall total in the digest (task-e auditor, 376-377); screened overall totals in the decision
  log (t3_critique.md 260-261).

Other:
- System `python3 -I -B`, stdlib only, no repository code, instead of `uv run` (precedent 2026-11-05 (g)(i)):
  heart_critique.md (353), t1.md (470), t1_critique.md (26), t3_critique.md (254-255).
- Writes beyond the one output file the brief allows: t1_critique's helper script in the scratchpad, deleted within a
  minute (236-237); task-d auditor's scratch copy of TRIES rows, deleted (280).
- Skill-map inversions on committed practice results of released files (allowed; Scores bars them for hidden targets):
  task-c (D11b's severity skill, 211 and 347); cross critic (fgw's mmd_u skill on the E7.5 backtest, 282-285).
- Barred classes met and not recorded: class B, two search summaries (task-c 274-275); class L, search summaries of two
  community repositories (task-e 276-277); class Y, two repository listings and the site's partner passages (task-e
  278-280); possible class A, the Task 1 and Task 3 pages, the dxyang42 blog and the proposal PDF p.7 (task-e 281-284);
  class W at many pages (task-a 214-221, task-a2 221-227, task-b 157 and 243-247, task-c 276-279, task-d 227-233,
  task-e 285).
- task-f's research agent cited docs/research/t3_methods.md (unscreened; still holds per-board figures, 2026-11-05 (i))
  for VCC's design; its log does not say how it knew; the auditor repointed it to the digest (task-f 9, 190, 226).
  Confirm it was not opened.
- The designers read and cited another agent's untracked work (experiments/gate/path_to_target.py and its JSON) as
  evidence before any review (t1.md 434-439; cross.md 49-50); the T1 critic objected (t1_critique.md 202-203).
- The design texts' 209 passages are struck (issue 5), on the precedent of (g)(v).
- A workflow restart created a second task-a agent; both notes were screened and audited (REPORT.md section 11). Not a
  slip.

## 4. This audit's own log

- Read: AGENTS.md; docs/DATA_POLICY.md; REPORT.md; notes/brief.md; notes/task-a, a2, b, c, d, e, f (b, c, d, e, f and
  a2 in full; task-a in full); notes/design/* (compliance sections, the 209 sections, candidates D, C1, C2, C4a);
  notes/registry.md (header, map, stats, targeted greps); docs/PLAN.md sections 7-8 and decision log lines 1150-1294;
  docs/ORGANIZER_QUESTIONS.md lines 1-45, 118-132, 245-277; docs/S2_RULES.md lines 1-82 and a grep;
  experiments/gate/path_to_target.py; its JSON (keys parsed with the system `python3 -I -B`, stdlib json only);
  src/vec/submit.py (grep of the clean-tree lines); experiments/coords/REPORT.md (grep; it printed released-stage
  geometry and "[redacted 2026-10-19: A]" marks, no values); experiments/research_2026_11/REPORT.md (grep).
- Git: `git status --short`, `git log --oneline -12`, `git ls-files` on three paths, all with GIT_OPTIONAL_LOCKS=0.
- Displayed, kept class: docs/PLAN.md section 8's gate figure (209.3); practice results in the decision log.
- Not opened: boards/, data/, runs/, docs/research/ (the digest included), any memory file, session log, transcript or
  harness-saved output. No web. The harness injected the auto-memory index into this agent's context; not opened.
