> Record for ledger row K27 (T3 candidate CHIPQ-S): the points a score-blind rule writer settled before the read's
> code existed, the corrections after a cold attack (section D) and a rule agent's addendum on the keyword filter
> (section E). Copied from the session's scratch at the row's registration commit; see `docs/PLAN.md`, decision log
> 2026-11-23 (i) to (o). The row in `experiments/harness/TRIES.md` governs.

# K27 (CHIPQ-S): points settled before the read's code exists

Authority, highest first: S = the challenger's final specification (drastic_spec_t3_1.txt); N = the challenger's notes;
D = the designer's candidate (its practice_kill_test is the base text of the rule); E = the database engineer's
hand-back (chipq_data_build.json); K22/K23/K24/K26 = those ledger rows. Tie-break used where no text fixes a point: what
K24 fixed for the same point; where K24 fixes none, the reading that reads less and refuses more. No constant is new.

## A. Points the texts leave open or state two ways

1. **Placebo draw (S: "one antigen drawn ... with seed 20261124 + s" against "ten seeds give ten different placebo
   genes").** Both statements are kept, and the draw method is what changes. Take the pool sorted by name (Python's
   default string order). Go through the seeds s = 0 to 9 (vec.adoption.SEEDS) in that order. For seed s, the
   remaining pool is the sorted pool minus the antigens already drawn for earlier seeds, and the antigen is
   remaining[numpy.random.default_rng(20261124 + s).integers(len(remaining))]. Source: N (c), "a placebo drawn fresh
   per seed gives ten different placebo genes, so the arm is tested against a null over antigens"; that sentence
   carries the inference. The seed formula is S's. The pool needs at least 10 antigens, which is D's P2 bar ("at least
   10 placebo antigens qualify"). Fewer than 10 stops the run at P2.
2. **The placebo pool's bars.** An antigen is in the pool when it has at least 2 experiments that the keyword filter
   keeps and |B_p| lies in [0.75 K_q, 1.25 K_q], both ends included. K_q = min(40, |B_q|) is the QUERY's K (Gata4 on
   this fold). |B_p| is the placebo's whole bound set, counted before any top-K cut. Excluded: q, Ctnnb1, Gata4 and
   Gata6. The band uses K_q because S step 4 defines K as the query's. "Both ends" comes from D's wording "a bound set
   of 0.75K to 1.25K genes". Each placebo is then built by steps 3 to 7 with its own K_p = min(40, |B_p|), at the same
   size S. The least pool size is 10 (point 1). The pool and its size are computed from the database only, uncounted
   (P2).
3. **Ties in the ranking by s_g.** Sort by s_g, highest first, with a stable sort. Genes with equal s_g keep the
   barred factor's file order (transfer_mab21l2_all.csv), both inside B and at the K-th place. A tied gene's weight
   comes from its position. Source: S step 3 fixes that file order for every table, and K24 computes everything "in
   file order". E's reading stands.
4. **The finish.** NO reference profile: finish's den is a zero profile. That means exp, then clip to [0.25, 4], then
   exactly 1 where w = 0, with no rescale, and nothing is rescaled or clipped after any product. The two texts:
   - "finished exactly as every outside table is finished": the earlier tables were divided by a den-weighted mean,
     whether pooled groups or whole_den.
   - S's own parenthetical lists only "set to exactly 1 at genes with w = 0, clipped to [0.25, 4]".

   The parenthetical yields nothing, and the general phrase yields its rescale, for two reasons:
   - S step 6 fixes the sign as "a factor's bound genes go down". A den rescale divides every value by a mean below 1,
     so a bound gene with w near 0.5 can end above 1 when the den is concentrated on bound genes. That breaks S's
     stated convention.
   - K24 fixes no finish for a new table (it used tables "as loaded"), so the tie-break applies, and no profile reads
     less.

   So FINISH_DEN = None stands. The target-free exponent then sets the table's size whatever the finish does.
5. **Tables folder, names, listing.**
   - Folder: ~/v/vec-ext/t3/tables_k27 on the DGX, outside every checkout, named to the script only by the
     environment variable VEC_T3_EXT_TABLES_K27 and checked by t3x_common.outside_repo. It is a folder of its own, so
     no K22 or K24 table sits in it. This follows K24's pattern of one variable per folder; K24 fixed no name for a
     new table.
   - Names: chipq_Gata4 (the query) and chipq_rand_seed0 to chipq_rand_seed9 (one placebo per seed; no antigen name
     in a file name, a manifest or a listing; the manifest records the target file's sha256 and the pool index only).
   - Listing: experiments/t3_ext/TABLES_K27.sha256, named and not yet filled at the row's registration commit, and
     filled by the row's second commit (after P2, before P3) with exactly those 11 lines ("<sha256>  <name>",
     sha256sum form, sorted by name), with K27_LISTING_SHA256 and the 11 table pins filled in k27_chip_query.py.
     Source: S steps 6 and 9; the two commits are D6.
   - `two` is not in this listing. It enters only as S fixes: K22's listing 884b294bb8ea... (experiments/t3_ext/
     TABLES.sha256) and table 9490f524e40e..., through K23's code path, from the folder VEC_T3_EXT_TABLES_TWO names
     (~/v/vec-ext/t3/tables_two_study), as K24's guard 2 read it.
6. **The chip_size twin.** One twin per size, of chip_q's own product P = x_e2 * C_q ** e. The twin is
   chip_size = k * phi_M ** a, with:
   - a = SD(ln P) / SD(ln phi_M);
   - k = den_total(P) / den_total(phi_M ** a);
   - population SD over the 500 genes in file order; den = WT E8.75's mean within-cell proportions (K23's formula, as
     K24 took it).

   a and k come from the pre-gate (P3 below), are recorded in pregate.json to full precision and are applied
   unchanged in both worlds. The twin holds no gene of C. The twin's product is the whole arm's factor, x_e2 times
   the table, because S's arm is x_e2 * C ** e and K24 took P as the whole arm's product. Source: S ("K24's size and
   depth twin, k * phi_M ** a with a and k from the pre-gate").
7. **The comparator of every pass clause.** Every paired difference d is taken per cell and paired by seed:
   - (a), (b), (c) and (e) take d against x_e2 (K23's named form, the base);
   - (d) takes chip_q against chip_rand and against chip_perm of the same size;
   - (d') takes chip_q against chip_size of the same size;
   - the size choice takes chip_q(4 S*) against chip_q(S*);
   - (g) takes each control (chip_rand, chip_perm, chip_size) of the size against x_e2 in every cell where chip_q is
     AHEAD (D1).

   barred_form is reported only. Source: N, "base x_e2", "a pass would not have shown it beats the form in force,
   which K24's clause (f) made the standard". On the size twin: S lists the twin as one of three controls, and N says
   they "separate size, depth and gene identity", while D's (d) names two. The twin is added at K23's and K24's bar
   (their (e)), which refuses more.
8. **The (a) clause.** D's (a) is "K22's QUALIFIED, unchanged": AHEAD in at least one plain cell and at least one dm
   cell, and on g1 in some world and on g2 in some world. K23's and K24's added gate sub-clause (an AHEAD dm cell whose
   gated count equals the base's) is also taken, because it refuses more and K24 fixed it for the same comparison.
9. **Either size qualifying.** Each size is decided on its own, with its own three controls and the same clauses
   ((a) to (g), D1). The
   row passes when at least one size qualifies. When both qualify, the larger is named only by K23's fixed pivot
   (point 7's last comparator). Source: K23's NAME rule, the one prior rule written for two exponents of one table.
10. **Gate clauses, counted or not.** By K26's ruling (a reading that compares a predictor with a fold target, or
    with the log of the barred factor, counts):
    - P1 and P2 read the database only and count nothing.
    - P3 computes corr(ln C_q ** e, ln phi_M) and the twin's a and k, which use phi_M, so P3 counts the eight tries
      from the start of its step, whatever follows. P3 also computes corr(ln C_q, ln two ** 2) for the second drop
      condition (D2); it changes nothing of the counting.
    - The exponent's own check (SD of ln C ** e equals S within a relative 1e-12) is P2's, because it uses the
      database only. It refuses, and counts nothing, as K24's check does.
    - K24's own table bounds (max abs ln 0.5, SD(ln product) at most 1.05 SD(ln phi_M), den 5%, 45 genes at the
      clip, e in [1/8, 8]) are NOT added. D's P2 fixes the gate's clauses, and K24 set those bounds for "no size above
      S*". Applied here they would refuse S's committed 4 S* on every table (max abs ln at 4 S* is about 0.95). They
      are reported, never used.
11. **Void versus fail.**
    - VOID: a guard fails, the tree is dirty, the commits or machines are mixed, a pin or listing is off, or the
      files disagree with pregate.json (K24's list). A void counts the eight tries, as K26 rules, because P3 runs
      before the guard part. Each guard arm that fails adds one (K24), so at most ten.
      A rerun after a guard failure finishes at the same guard exit with no scoring part started, and the VOID reading
      is made once, only when the reading folder holds no reading.json (K24's clause; D4).
    - A refusal before P3 (pins, listing, kept list, manifest, scaled SD, a duplicated placebo table or listing sha256)
      counts nothing.
    - An INCOMPLETE read is not a reading, and the eight tries stand counted from P3; a rerun past it happens only under
      rows K21 to K24's read-crash and code-fix rules and makes no new reading and no new count (D3).
    - FAIL: the clauses are read and nothing qualifies. It counts the eight tries.
12. **The size and e at the test.** A pass names a SIZE (S* or 4 S*), not Gata4's exponent. A ruled build computes
    Ctnnb1's e by the same formula, S / SD(ln C_Ctnnb1), from Ctnnb1's own table, and records it to full precision.
    Gata4's e never enters a Ctnnb1 file. S's "e ... the committed pregate.json records" cannot hold for a table that
    the pre-gate never opens (E: the tables step "never opens Ctnnb1").
13. **The keyword filter's scope for placebo antigens.** It runs on every non-query mm10 antigen's experiment rows,
    over the seven metadata fields E names, with classes (i) to (iii) and class (ii) included ("the same drops"). It
    keeps an antigen with at least 2 kept experiments. No class is added. Word variants and abbreviations of the
    challenger's own words are matched. A named agent without one of the challenger's words is not dropped, since the
    filter adds no words.

## B. The database engineer's ten literal readings

1. POOL_BAND (0.75, 1.25) x the query's K, both ends included, each placebo with its own K: **stands** (A2; S step 4,
   D's "0.75K to 1.25K").
2. Placebo draw by independent default_rng(20261124+s).integers(len(pool)): **changes to** the sequential draw without
   replacement of A1 (same seeds, same sorted pool, already drawn antigens removed before each draw; the pool must hold
   at least 10; source N (c) and D's P2).
3. FINISH_DEN = None: **stands** (A4).
4. Ties keep file order: **stands** (A3).
5. Bound means 2 x (count > 0) >= n and at least 1: **stands** (S step 4).
6. Table names chipq_{q}, chipq_rand_seed{s} into a free --out: **changes to** the fixed folder
   ~/v/vec-ext/t3/tables_k27, given by VEC_T3_EXT_TABLES_K27, with the names kept and the committed listing
   experiments/t3_ext/TABLES_K27.sha256 of those 11 tables pinned in the script, both filled by the row's second
   commit (A5; S steps 6 and 9; D6).
7. P1 bars 3 and 15, "0 0" and exit 3 on an empty kept list: **stands** (D, N).
8. Keyword filter as built: **stands** (A13). The class-(ii) words apply to placebo rows, and the abbreviations and
   spellings are of the challenger's words. The case-sensitive sh/si prefix stands: a case-insensitive prefix would
   match ordinary words, not the word quoted. No word is added by this row.
9. Gata6 filtered like any antigen and excluded from the pool only: **stands** (S).
10. chipq_kept.yaml form, committed and unchanged against HEAD, with its pin in CHIPQ_KEPT.sha256: **stands** (E's
    departure; the pin's place changes nothing of the rule).

## C. Open, for the main session (not settled here)

- Which agents read the Ctnnb1 and Gata4 rows is fixed by the main session: two independent screening agents, neither
  of them the main session (a declared departure from the specification, D11). The row only requires that
  chipq_kept.yaml is committed before P1.
- K24's closure ("no further size, weight or combination of these tables") is read as not reaching a new database
  table on top of x_e2. That is S (13)'s reading, recorded by the main session when it places the row, as K26's was.

## D. Corrections after the cold attack

Source: the cold challenger's attack (k27_rule_attack.json, verdict "place after the listed corrections"). Each of its
thirteen faults is applied with its exact_correction words at the place it names; the reviser added no constant,
loosened nothing and changed no clause the challenger did not fault. Cell numbers are the ledger's column order, 1 to 11.

1. **A control that is itself AHEAD (cell 8, QUALIFIED).** Added clause (g) after (f): in every AHEAD cell no control
   of that size (chip_rand, chip_perm, chip_size) is itself AHEAD of x_e2; a size at which any control is AHEAD of x_e2
   in a cell where chip_q is AHEAD does not qualify. Can only make naming harder.
2. **A table that duplicates `two ** 2` (cell 8, P3).** The drop clause now also drops the arm and its controls when
   abs corr(ln C_Gata4, ln two ** 2), over genes != 1 in either, is 0.5 or more or undefined (the same bar, no new
   constant; neither correlation depends on the size).
3. **INCOMPLETE vs the counting sentence (cells 1 and 8).** After the exit-1 INCOMPLETE clause: an INCOMPLETE is not a
   reading, the eight tries stand counted from P3, and a rerun past it happens only under rows K21 to K24's read-crash
   and code-fix rules and makes no new reading and no new count. A11 above carries the same words.
4. **VOID once-only (cell 8).** Added K24's clause after the guard-failure sentence: a rerun after a guard failure
   finishes at the same guard exit with no scoring part started, and the VOID reading is made once, only when the
   reading folder holds no reading.json.
5. **x2 asserted before the fetch (cells 3 and 9).** Cell 3's "x2 NOT TESTABLE and not run (the analysis list holds no
   ... antigen among its 733 mm10 antigens)" is replaced by: x2 is not read by this row; P2 records whether the
   analysis list holds an antigen for the barred factor's gene, and either way no x2 arm is built, no x2 try is counted
   and no x2 cell is read. Cell 9 point (4) now reads "the database is not read for the barred factor's gene by this
   row". The row no longer rests on the database engineer's format inspection (733 mouse antigens with target-gene
   files, Ctnnb1 and Gata4 held, none for the barred factor's gene); P2 records that point itself at the run.
6. **Hard-coded try numbers (cell 1).** "(74), so T3 75 to 82" is replaced by: from the running count in the counts
   table at placement, the eight tries take the next eight numbers after it, which the outcome cell records.
7. **Placement (cell 10).** The challenger's first form is taken: two commits. The registration commit, before P1,
   carries the script, chipq.py, chipq_fetch.py and their tests on synthetic data and hand-made tables, with
   TABLES_K27.sha256 and the 11 table pins named and not yet filled, and X018 placed in EXTERNAL_DATA.md before any
   fetch-to-model step. The second commit, after P2 and before P3, fills the listing's 11 lines, its sha256 and the 11
   pins in the script; every counted step runs at it. Nothing built from released data and nothing scored at either.
   (The alternative of moving the table build to the X018 step was not taken.)
8. **J5's own label (cell 9).** Added after J5's arithmetic: no counting task runs the entry on T3 and the four cells
   are one truth, so the read is a diagnostic under J5; the adoption rule counts nothing of it, and a pass is a
   labelled reading carried to the captain, never a kept change.
9. **What the controls do not separate (cell 9).** Added before Rule 9: chip_rand matches size and the bound-set band,
   not the number of weighted genes; chip_perm is one fixed shuffle and below the 0.5 corr bar does not separate
   residual alignment with phi_M from gene identity; chip_size matches size and depth, not the gene pattern; at 4 S*
   the per-gene move and den shift exceed every bound K23 and K24 applied, so a 4 S* pass is read as a larger move than
   any table before it, not as a like-for-like step.
10. **Reported figures never re-choose (cell 8).** Added after "K24's recorded-not-applied bounds": the recorded bounds
    and every reported figure never change the named size, the pass, the file named or a later row's or ruling's
    choice between the two sizes; the size a pass names is the one its fixed pivot gives and nothing else.
11. **The screen's readers (cell 3).** The row now declares the departure in its own words: two independent screening
    agents, neither of them the main session, against the specification's "the main session and an independent agent",
    so that nobody in the session reads the metadata rows.
12. **Duplicated placebo tables (cell 8, P2).** Added: a placebo whose finished table equals the query's byte for byte,
    or any sha256 appearing twice among the 11 listing lines, refuses the run (exit 2, nothing counted, before P3).
13. **Targets unchecked (cell 5).** The parenthesis now records the main session's check: checked by the main session
    on 2026-11-23 against this ledger's T3 headroom block at commit 56a357c (response identity 15.18 points left on the
    cross-knockout stand-in and 11.89 on stand-in A, as the cell has them); no clause of the reading turns on the figure.

### The two additions the attack found undeclared, now declared

- **The readers (attack item 10).** The screen's readers are two independent screening agents, neither of them the main
  session: a departure from the specification's "two readers, the main session and an independent agent". It reads
  less inside the session (nobody in the session reads the metadata rows); either reader may still drop, and a doubt
  still drops.
- **Clause (f) (attack item 11).** The de_score plus de_direction skill check against x_e2 in every AHEAD cell is
  carried over from the designer's text (its (e)) and from row K26's (e). It is an addition over the specification's
  QUALIFIED that can only make naming harder, as A8's additions are. (Attack item 12, the x2 assertion, is removed by
  D5 rather than declared.)

### Consistency edits the corrections force (each with its source)

- D2: P3's written figures now hold corr(ln C_Gata4, ln two ** 2) beside corr(ln C_Gata4, ln phi_M); "It counts
  because that corr" becomes "the first corr"; the old "(the corr does not depend on the size)" is removed as the
  corrected clause carries "neither correlation depends on the size"; the drop sentence now reads "not independent of
  the barred factor or of the named form's cross-knockout table", true for either correlation; A10 notes the second
  corr.
- D5: P2's printed counts include whether the analysis list holds an antigen for the barred factor's gene (1 or 0),
  which it records; the reported list carries it; cell 8's "x2 cannot be read" and cell 9's J4 "x2 not testable" now
  say "not read by this row".
- D3 and D6: cell 1 says the eight tries take the next eight numbers and count also for an INCOMPLETE after P3 (no new
  reading or count on a rerun); cell 8's counting sentence adds "an INCOMPLETE after P3".
- D7: cell 1 names the registration commit (P1, P2) and the second commit (the eight tries); cell 3's listing sentence
  and A5 and B6 say the listing and pins are filled by the second commit; cell 8's Order names the commit of each step,
  P2's "listing lines, which the second commit fills", and VOID's "one commit (the row's second commit)"; cell 11
  records both commits and the eight try numbers.
- D1: EITHER SIZE QUALIFYING judges clauses (a) to (g); the no-qualify sentence now reads "does not qualify against
  the named form and its controls at either size" (it was "does not beat the named form", false when (g) fails it);
  cell 9's placebo item names (g); A7 and A9 above name (g).
- D4 and D12: A11 above carries the once-only VOID clause and the new refusal.
- D11: section C's first point names the departure.

## E. The keyword filter: the rule agent's addendum to A13 (2026-11-23, after the database reviews)

A13, addendum (rule agent kwdecide; no class and no word added): the following are matched as spellings, notations or
abbreviations of the SCREEN paragraph's own words, and the filter must drop them: Tg (transgene); esiRNA and the sh/si
control names shctrl, sictrl, shscr (sh/si/shRNA; B8's refusal of a broad case-insensitive sh/si prefix still stands);
fl/- and f/f (flox); -/+ (mutant, the reversed +/- notation); ex3, delta-ex3/Δex3/dex3, del(ex3) (exon-3-deleted,
class ii; a bare ex3 drops as a doubt); mutagenized (mutant; a doubt drops); Cas9 and dCas9 (CRISPR, since Cas stands
for CRISPR-associated); point-mutation names of the residue-number-residue form (mutant/mutation, a notation like
+/-). Not matched: sgRNA alone (single guide RNA has no form of CRISPR; per A13 a named agent without the challenger's
word is kept). Reverse case: a metadata field whose entire value states no treatment ('treatment: none',
'treatment=untreated', 'no treatment') states Tier U's 'untreated' and names no class (iii) treatment, so the row may
be kept; only that whole value is exempt, every other occurrence of 'treat' still drops, and 'vehicle', 'DMSO',
'control', 'n/a' or an empty value stay dropped as a doubt. Named drugs, ligands and growth factors without a class
word remain kept (already settled). Under the committed filter (branch campaign/t3-chipq-2), all the matched tokens
above are kept, and all the reverse-case rows are dropped. Any placebo_kept.json written before the filter changes is
refused (fix item 4).

After two reviews of the patterns three exact fixes are being applied (Cre written straight after a capital-letter driver
name; point-mutation names with two or more position digits only, so that look-alike antigen, strain and agent names
stay kept; point-mutation names joined to a gene symbol). No word and no class is added by them.
