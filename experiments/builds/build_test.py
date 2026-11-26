#!/usr/bin/env python
"""Test-board candidate files, built with the harness's own builders in the test geometry (P3_RUNBOOK section 5).

Each entry is exactly the candidate the harness judged, run on the test's stages, then the test-only steps the
pseudo-tasks cannot run (the pre-registered sizes). Nothing is uploaded: the file is written with its build stamp and
checked locally (`vec check`); a server format check is a separate `vec submit ... --from-dgx` call from the local
machine, and an official upload needs the captain's go-ahead.

  t1      M1 (docs/PLAN.md 4.1; the captain's rulings of 2026-11-01): E10.5 cells re-weighted to E10.5's own chamber
          mix (or --cm-share, the cardiomyocyte share of heart cells), Feng's E9.5 -> E14.5 direction on the cells of
          Feng's groups at 2/5 of the bracket (E10.5 -> E12.5); the rest keep their E10.5 expression (R1 of
          experiments/harness/TRIES.md, the captain's ruling of 2026-11-05; --no-hold-unbracketed gives them B0's
          damped E9.5 -> E10.5 momentum, M1 as built before R1). R2 of TRIES.md also holds R2's gene modules flat in
          Feng's direction (the union N12 wrote, pinned in configs/t1_hygiene.yaml): on by default with R1 since the
          build commit after N12 (docs/PLAN.md decision log 2026-11-06 (k)); --no-feng-hygiene rebuilds the s1 of
          88b057f (M1 with R1, before R2). T12a and T12b (TRIES.md; the captain's rulings of 2026-11-10 (a) and (d)):
          --strong-null adds the strong-null down-tail prior after M1's steps (vec.t1_scope.apply_down_prior, genes by
          the reference's level alone), off by default and refused until configs/t1_strong_null.yaml says approved
          (a later build commit of the main session, after the stop rule and the captain's ruling). --g-mom (default
          off; only by the captain's ruling on try T11 of TRIES.md; refused together with --strong-null, whose order
          against it no row fixes) adds g_mom's within-type factors: fitted on E9.5 -> E10.5 and applied last, after
          Feng's direction, on the cells g_mom edits (metacell factors through apply_factors with keep_depth, R1's held
          labels left alone), so no share, cell count or scope changes. Row R3b (the captain's ruling of 2026-11-10
          (p)): the factors edit cardiomyocytes and endocardium only; epicardium and vascular endothelium keep M1's
          counts. --feng-groups none (default authors, every file as built; only the captain's ruling puts it in an
          upload) builds one file: the first upload's command (t1 --g-mom --base lib3) with Feng's direction taken
          out and nothing put in its place. Every cell keeps its E10.5 expression before g_mom's edit (R1 on every
          label), and no Feng file, R2 list or other outside data is read. g_mom's mask is the readout config's own
          (experiments/t1_direction/config.yaml): the factors edit the cm and endocardium cells whose label its
          regexes give a group, which on the labels configs/t1_chamber_classes.yaml places in those categories is
          every such cell but Early VIC, the same cells as in the entry. vec.t1_state_transport.R1_FENG_GROUPS, a
          list fixed in N12 by counting Feng cells, is not read under none
          (tested), and the hook's covered check is a self-check. Refused, before anything is read, without
          --g-mom and with --base pooled or mid, --cm-share, --blood-share, --lib3-pair, --herd,
          --no-hold-unbracketed, --feng-hygiene or --strong-null. The build record then carries a feng_groups field
          saying the direction step was removed, and its fraction is null. --released-bet share_trend | advance with
          --bet-strength (try T61 of experiments/harness/TRIES.md; off by default, --feng-groups none only): that file
          with its cells drawn at the chamber mix moved along the E9.5 -> E10.5-portions trend, or drawn inside each
          category towards the released momentum axis (vec.t1_bets). Both change composition, a lever the
          post-exposure block closes: refused unless --captain-ruling names the captain's ruling row for the arm.
          After a permitted-only relabelling whose record routes to recomputed, none is also refused unless
          --feng-none-ruling names a later decision-log entry of the captain (naming N54 and containing "sets the
          N54 route aside", not the run's own ruling, made after the record's commit) that sets the route aside;
          before any relabelling, none builds as before.
          --feng-groups recomputed (off unless the captain rules): the same one command, t1 --g-mom --base lib3, with
          Feng's direction formed on the groups experiments/external/feng_relabel.py recomputes from the permitted
          cells' raw counts (register rows X002 and X010), R1 and R2 as in the first upload. Refused, before
          anything is read, with any other base or option of that list (and with --no-feng-hygiene), and unless the
          labeller's record is committed at HEAD with the route recomputed; then refused unless the recomputed file
          passes feng_dirs.check_recomputed. R1's hold and g_mom's mask follow the groups with a recomputed direction
          after the fallback: a group without one moves no cell and g_mom leaves it alone. The record's feng_groups
          field names the file, its sha256, the yaml's sha256, the route, the fallbacks and the gate.
          The candidates fixed in docs/PLAN.md decision log 2026-11-20 (i) (off by default; no ruling; only the
          captain's ruling puts one in an upload), each a single change to that one file and refused on any other
          command or seed and on a dry run: --g-mom-damping 1.0 (candidate A: g_mom's factors finished at damping 1.0
          in place of B0_DAMPING 0.5, as practice arm gmom_s1 of experiments/t1_m1/bets_read.py finishes them) and
          --feng-fraction-mult 1.5 (candidate B: the direction step at 1.5 times the calendar fraction 2/5, as
          experiments/t1_m1/double_step_read.py multiplies it); both together are candidate C. Their defaults (0.5,
          1.0) build today's bytes; given at all, the record's t1_candidate field names both values.
          --feng-fraction-mult 2.0 --captain-ruling R20 (row R20 of experiments/harness/TRIES.md, the captain's
          order of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (x)): the same file with the direction step at 2.0
          times the calendar fraction, nothing else changed; refused without that flag, with any other row, unless
          row R20 is the captain's order naming the option, with --g-mom-damping other than 0.5, with --released-bet,
          on any other command or seed and on a dry run. Its record's t1_candidate field names row R20.
          --feng-fraction-mult 2.5 or 3.0 --captain-ruling R22 (row R22 of experiments/harness/TRIES.md, the captain's
          order of 2026-11-22, docs/PLAN.md decision log 2026-11-22 (j): the first search line): the same file with the
          direction step at 2.5 or 3.0 times the calendar fraction, nothing else changed; refused as the 2.0 value is,
          with R22 in place of R20 (and 2.0 refused with R22). Its record's t1_candidate field names row R22.
          --g-mom-damping 1.5 or 2.0 --captain-ruling R22 (row R22, the second T1 search line of docs/PLAN.md decision
          log 2026-11-22): the same file with g_mom's within-type damping at 1.5 or 2.0, nothing else changed;
          refused without that flag, with any other row, unless row R22 is the captain's order naming the option,
          with --feng-fraction-mult other than 1.0, with --released-bet, on any other command or seed and on a dry
          run. Its record's t1_candidate field names row R22; --g-mom-damping 1.0 (candidate A) needs no ruling, as
          before.
          --g-mom-damping 1.0 --feng-fraction-mult 2.0 or 2.5 --captain-ruling R22 (row R22, the third T1 search line
          of docs/PLAN.md decision log 2026-11-22): file C's command with the direction step at 2.0 or 2.5 times the
          calendar fraction, g_mom's damping held at 1.0, nothing else changed; refused without that flag, with any
          other row, unless row R22's first cell is the captain's order naming the combined option
          (`--g-mom-damping 1.0 --feng-fraction-mult 2.0` or `... 2.5`; the bare multiplier's option does not count),
          with --released-bet, on any other command or seed and on a dry run. The pair is tested before the
          multiplier's own routing; every other pair of the two values is taken or refused as before. Its record's
          t1_candidate field names row R22.
          --g-mom-damping 1.5 or 2.0 --feng-fraction-mult 2.0 --captain-ruling R22 (row R22, the fifth T1 search line,
          files T1_s2rcB2d15 and T1_s2rcB2d20): the 2.0 file's command with g_mom's within-type damping at 1.5 or 2.0,
          nothing else changed; refused without that flag, with any other row, unless row R22's first cell is the
          captain's order naming the combined option (`--g-mom-damping 1.5 --feng-fraction-mult 2.0` or `... 2.0
          --feng-fraction-mult 2.0`; the bare options do not count), with --released-bet, on any other command or seed
          and on a dry run. The pair is tested beside the third line's, before the multiplier's routing and the
          damping's; every other pair of the two values is taken or refused as before. Its record's t1_candidate field
          and g_mom's rule text name the fifth line.
          --g-mom-damping 1.5 or 2.0 --feng-fraction-mult 2.5 --captain-ruling R22 (row R22, the sixth T1 search line:
          T1 line 6 of docs/PLAN.md decision log 2026-11-24 (bb), the captain's order of ruling row R40; files
          T1_s2rcB25d15 and T1_s2rcB25d20): the 2.5 file's command with g_mom's within-type damping at 1.5 or 2.0,
          nothing else changed; refused without that flag, with any other row, unless row R22's first cell names the
          combined option (`--g-mom-damping 1.5 --feng-fraction-mult 2.5` or `... 2.0 --feng-fraction-mult 2.5`; the
          bare options do not count), with --released-bet, on any other command or seed and on a dry run. The pair is
          tested beside the fifth line's, before the multiplier's routing, which keeps refusing 2.5 with any damping
          but 0.5 (and 1.0, the third line's pair); every other pair is taken or refused as before. Its record's
          t1_candidate field and g_mom's rule text name line 6.
          --g-mom-damping 1.0 --feng-fraction-mult 3.0 --captain-ruling R22 (row R22, the seventh T1 search line: T1
          line 7 of docs/PLAN.md decision log 2026-11-24 (bb), the third line's one further half step, the captain's
          order of ruling row R39; file T1_s2rcCm30): file C's command with the direction step at 3.0 times the
          calendar fraction, g_mom's damping held at 1.0, nothing else changed; refused without that flag, with any
          other row, unless row R22's first cell names `--g-mom-damping 1.0 --feng-fraction-mult 3.0`, with
          --released-bet, on any other command or seed and on a dry run. The pair is tested before the multiplier's
          routing, which keeps refusing 3.0 with any damping but 0.5. Its record's t1_candidate field names line 7 and
          the captain's order; g_mom's rule text stays candidate A's.
          --direction-momentum all --captain-ruling R23 (row R23 of experiments/harness/TRIES.md, the captain's order
          of 2026-11-23, docs/PLAN.md decision log 2026-11-23 (c): candidate MS): R13's file with one change, the
          drawn cells whose label's group has a recomputed direction also get B0's per-category damped mean momentum
          (harness.B0_DAMPING, pseudocount 0.5, keep_depth, stream [seed, 1], the build's held genes; fitted from E9.5
          to E10.5's portion libraries 1 and 2, library 3 never in the fit), applied before the atlas step ([seed, 2],
          the calendar fraction) and g_mom ([seed, 5]), which stay as they are; cells without a group keep R1's hold.
          As the practice arm m1rc_gmom_mom of experiments/t1_m1/ms_read.py (row T65) builds it, the fit on the test's
          stages. Refused without --captain-ruling R23, with any other row, unless row R23's first cell is the
          captain's order naming the option, with --g-mom-damping, --feng-fraction-mult, --released-bet,
          --bet-strength or --no-hold-unbracketed, on any other command or seed and on a dry run. Its record's
          t1_candidate field names row R23, the momentum's constants and fit stages, and the command
          --feng-fraction-mult 1.5 --direction-momentum all --captain-ruling R26 (row R26 of
          experiments/harness/TRIES.md, the captain's order of 2026-11-23, "combine MS"; docs/PLAN.md decision log
          2026-11-23 (j), the T1 combined MS list's one file T1_s2rcMSB): candidate MS's step exactly as row R23 fixes
          it, on candidate B's command exactly as B has it (the atlas step at 1.5 times the calendar fraction, g_mom's
          damping 0.5), nothing else changed. Refused unless row R26's first cell is the captain's order naming
          `--feng-fraction-mult 1.5 --direction-momentum all` and `--captain-ruling R26`, with --g-mom-damping (any
          value), --released-bet, --bet-strength or --no-hold-unbracketed, on any other command or seed and on a dry
          run; the pair without R26 is refused as row R23's path refuses it, --direction-momentum all with any other
          multiplier stays refused, and R26 with any other option set is refused. Its record's t1_candidate field names
          row R26, the list's file, B's two values, the momentum's record and the command.
          --base pooled --feng-groups recomputed --captain-ruling R37 with --g-mom (row R37 of
          experiments/harness/TRIES.md, the captain's named exception of 2026-11-24, "Yes, T1 pooled base, my
          exception"; docs/PLAN.md decision log 2026-11-24 (z) and (bb); file T1_s2rc_pooled): R13's command on the
          pooled base, every E10.5 chamber cell at the stage's own mix in place of library 3's cells at library 3's
          mix, nothing else changed. --feng-groups recomputed keeps refusing --base pooled without --captain-ruling
          R37; with it, the base is refused unless row R37's first cell is the captain's named exception ("ruling, no
          run") naming `--base pooled` with `--feng-groups recomputed` and `--captain-ruling R37`, with
          --g-mom-damping, --feng-fraction-mult, --direction-momentum, --released-bet or --bet-strength, with another
          seed and on a dry run; R13's command's own guards hold every other option as before. Its record's
          base_exception field names row R37, the order text and the command. The pool's depth record
          (docs/DATA_POLICY.md 2026-11-04 (f)(3)) is experiments/builds/t1_pool_depth_record.py's, from released E10.5
          alone.
          --feng-fraction-mult 2.5 --held-momentum 0.5 | 1.0 | 2.0 --captain-ruling <ROW> (the T1 design round of
          2026-11-24, files F1, F7 and F8 as its cold challenge fixes them; files T1_s2rcB25H05, T1_s2rcB25H10 and
          T1_s2rcB25H20): the banked file T1_s2rcB25's command (R13's with the atlas step at 2.5 x the calendar
          fraction, g_mom's damping 0.5) with one change. The drawn cells labelled Blood, which R1 holds at their drawn
          counts, get Blood's own released change: B0's per-category damped mean momentum (vec.t1_scope.momentum_shift
          with _t1_ms_step's constants: eps 0.5, min_cells 20, keep_depth, the build's held genes, R2's union with the
          pre-registered list, at f = 1) fitted on Blood-labelled cells at both ends (within={"blood": ["Blood"]}) from
          E9.5 to E10.5's portion libraries 1 and 2 (library 3 never in the fit), at damping 0.5, 1.0 or 2.0, on its
          own stream [seed, 6], after the atlas step and before g_mom, which leaves Blood alone. Macrophage, Early VIC
          and every other cell are the banked file's, byte for byte. It reverses R1 for Blood in that one file, so it
          is taken only with --captain-ruling <ROW>, ROW a row ID given on the command line (R, digits, an optional
          letter), while experiments/harness/TRIES.md holds exactly one row of that ID whose first cell is the
          captain's ruling ("<ROW> (the captain's ...", "ruling, no run") naming the file's option string
          (`--held-momentum 0.5`, `--held-momentum 1.0` or `--held-momentum 2.0`) and `--captain-ruling <ROW>`, each
          in backticks. Refused before anything else is read with any other value (the literal strings 0.5, 1.0 and
          2.0 only), with --direction-momentum (Blood's step with MS, file T1_s2rcB25MSH05, is named but not coded),
          with --g-mom-damping, --released-bet or --bet-strength, without --feng-fraction-mult 2.5 or with another
          multiplier, with --no-hold-unbracketed, on any other command or seed and on a dry run; and refused, still
          before the draw, unless E10.5's libraries other than 3 are 1 and 2 and Blood's own cells reach min_cells at
          both ends of the fit (released stages only, so the step is never idle). Its record's t1_candidate field
          names the value, the row, the file, the order text, the command and the step's record (its fit: stages,
          libraries, label, status; no cell number).
          --feng-fraction-mult 2.5 --direction-momentum all --captain-ruling <ROW> (the T1 design round of 2026-11-24,
          file F2 as its cold challenge fixes it; file T1_s2rcB25MS): candidate MS's step exactly as row R23 fixes it
          (_t1_ms_step, unchanged) on the banked file's command (the atlas step at 2.5 x the calendar fraction, g_mom's
          damping 0.5), nothing else changed; a third MS file beside rows R23's and R26's. Taken only while TRIES.md
          holds exactly one row <ROW> whose first cell is the captain's ruling naming `--feng-fraction-mult 2.5
          --direction-momentum all` and `--captain-ruling <ROW>`; refused with --g-mom-damping, --held-momentum,
          --released-bet, --bet-strength or --no-hold-unbracketed, on any other command or seed and on a dry run.
          Without such a row the pair is refused as row R23's path refuses it, as before; R23's and R26's paths are
          unchanged. Its record's t1_candidate field names the row, the file, the banked file's two values, the
          momentum's record and the command.
          Row R51 of experiments/harness/TRIES.md (the captain's order of 2026-11-24, decision log 2026-11-24 (iii);
          T1's line 10): seven files, each one change to T1_s2rcMSB's command (row R26's), each taken only with
          --captain-ruling R51 while TRIES.md holds exactly one row R51 whose first cell is the captain's ruling
          ("ruling, no run") naming `--captain-ruling R51` and the file's exact option string, each as a whole
          backticked string (no string passes on another that contains it): --feng-fraction-mult 2.0 or 3.0
          --direction-momentum all (T1_s2rcMS20, T1_s2rcMS30: MS's step at that multiplier); --feng-fraction-mult 1.5
          --direction-momentum all --g-mom-damping 1.0 or 1.5 (T1_s2rcMSBd10, T1_s2rcMSBd15: g_mom's within-type
          damping); --feng-fraction-mult 1.5 --direction-momentum all --held-momentum 0.5 or 1.0 (T1_s2rcMSBH05,
          T1_s2rcMSBH10: row R44's step, Blood-labelled cells only, fitted on Blood-labelled cells at both ends from
          libraries 1 and 2, library 3 never in the fit, its own stream [seed, 6], after MS's step and before g_mom);
          --feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.0 (T1_s2rcMSBm10: MS's own
          damping passed to _t1_ms_step in place of harness.B0_DAMPING; the option takes the literal 1.0 only and is
          refused with every other command), all on t1 --g-mom --base lib3 --feng-groups recomputed at the default seed
          on the test board; nothing else changed. Tested before every other route: a command that is not exactly one of
          the seven under such a row (the struck pair --feng-fraction-mult 2.0 --g-mom-damping 1.0 with MS included) is
          taken or refused exactly as before. One of the seven is refused with --released-bet, --bet-strength or
          --no-hold-unbracketed, on any other command or seed and on a dry run. Its record's t1_candidate field names
          the row, the file, its option string, the base file, the multiplier, g_mom's damping, Blood's damping (or
          null), MS's damping, the order text, the command, MS's step's record and, in the H files, Blood's step's
          record (held_momentum_step); g_mom's rule text names row R51 at damping 1.0 and 1.5.
  heart B0 (copy of E10.5 plus damped E9.5 -> E10.5 momentum in count space), sized by the pre-registered heart
          procedure: --entry s1 = scope hypothesis H1; s2_trim = the s2 of docs/S2_RULES.md (H2 as a trim to half
          E10.5's scope); s2 = H2 as first built (size only, superseded) (approved rule 3, `size_or_hold`). --trend
          (both heart files since the captain's ruling of 2026-11-13, docs/PLAN.md decision log (w); row H11 of
          experiments/harness/TRIES.md): b0_trend, each lineage's momentum the least-squares trend over E8.25_late to
          E10.5 in place of the E9.5 -> E10.5 step alone. --expression copy (off by default; --entry s1 only, refused
          with --trend, --depth, --herd and --cfm): the copy baseline's cells and expression (the draw B0 starts
          from, `vec.harness._copy`), no momentum. --z-mode raw (off by default): the size step's isotropic mode,
          every axis scaled by one factor to the hypothesis's target RMS, in place of the mode the size JSON
          recorded. The captain's direction for the heart's second upload (2026-11-18, docs/PLAN.md decision log
          2026-11-18 (k)): heart --entry s1 --expression copy --z-mode raw. The scope bet (no ruling; off unless all
          three flags are given; scope is a lever of the post-exposure block, so only the captain's ruling puts it in
          an upload): heart --entry s2_trim --expression copy --z-mode raw, the same two fixes on the half trim (H2,
          docs/S2_RULES.md), sized to the raw mode's own logged H2 figure. The count is s1's: vec.t2_build.n_cells
          does not read the pool's size (it only annotates the record), so the trimmed pool gets the count the server
          refused on 2026-11-18, and the file is halved with halve.py as the ruled file is. The full trim (no ruling;
          off unless all three flags are given; refused with --trend, --depth, --herd and --cfm): heart --entry
          core_trim --expression copy --z-mode raw, the same draw from E10.5 under radial_trim_mask at fraction 0,
          sized to the size JSON's raw-mode core forecast (z_modes.raw.pred_rms_core, the cloud at the core's own
          scope); its record and the half trim's carry scope_option (scope rule, kept-share field, size field,
          halving). The organ cut, file D (row R21 of experiments/harness/TRIES.md, the captain's order of
          2026-11-22 to code and read it; off unless all three flags are given; refused with --trend, --depth, --herd
          and --cfm): heart --entry organ_scope --expression copy --z-mode raw, the same draw from E10.5 under
          vec.t2_build.organ_scope_mask (the chamber categories of configs/t2_chamber_classes.yaml joined to the core
          at L, the 99th percentile of the core's nearest-core distances; an unmapped label refuses), sized by
          vec.t2_build.size_organ_scope to z_modes.raw.pred_rms_core times rms(organ pool) / rms(core) on E10.5; then
          halve.py at seed 1000001 as for the trims. Its record names the scope, the label config's sha256, L, the
          pool's E10.5 rows and the size field. The heart search line (row R28 of experiments/harness/TRIES.md, the
          captain's order of 2026-11-23, "start the heart search"; docs/PLAN.md decision log 2026-11-23 (t); off by
          default): heart --entry core_trim --trend-strength 0.5 | 1.0 --z-mode raw --captain-ruling R28, files TS05
          and TS10, the full trim C's command with row H28's trend step at that strength (the b0_trend factors and
          momentum_shift's fallback damping both times the strength, a test budget of 0.5 x the strength day) before
          the size step, on C's draw, pool, count rule, size field and seed; then halve.py at seed 1000001 as for C.
          At strength 1.0 the expression is vec.harness._b0_t2(ctx, trend=True) on the core-trim context bit for bit.
          Refused before anything is read with any other value, entry or z mode, with --expression copy, --trend,
          --depth, --herd and --cfm, on any board but the test board (no dry run), without --captain-ruling R28 or
          with another row, and unless TRIES.md holds exactly one row R28 whose first cell is the captain's order
          naming `--trend-strength 0.5`, `--trend-strength 1.0` and `--captain-ruling R28`; --captain-ruling on the
          heart without --trend-strength is refused too. Its record's trend_strength field names the option, the
          strength, the file, the ruling row, the order text, the command and the budgets. File D's line (row R29 of
          experiments/harness/TRIES.md, the captain's order "i want it too"; docs/PLAN.md decision log 2026-11-23
          (w); off by default): heart --entry organ_scope --trend-strength 0.5 | 1.0 --z-mode raw --captain-ruling
          R29, files OTS05 and OTS10, D's command with --expression copy replaced by the same step at that strength
          before the size step, on D's draw, organ-cut pool, count rule, size procedure and seed; then halve.py at
          seed 1000001 as for D. The rows, coordinates, cell number, size field and half's rows are D's. Every refusal
          of row R28's line holds here too, and core_trim with R29 or organ_scope with R28 is refused; it is refused
          unless TRIES.md holds exactly one row R29 whose first cell is the captain's order ("ruling, no run") naming
          `--entry organ_scope`, `--trend-strength 0.5`, `--trend-strength 1.0` and `--captain-ruling R29`. Its
          record's trend_strength field names the option, the strength, the file, R29, the order text and the
          command; its scope_option field stays D's own. File D's two further strengths (the design round of
          2026-11-24, its FAST proposal for the heart with 0.25 dropped by the round's challenge; off by default):
          heart --entry organ_scope --trend-strength 0.75 | 1.25 --z-mode raw --captain-ruling <ROW>, files OTS075
          and OTS125, the same step on D at those strengths, with OTS05's draw, organ-cut pool, count rule, size
          procedure, seed and half (halve.py at seed 1000001); the rows and coordinates are OTS05's. <ROW> is given
          on the command line, a ruling row's id (R and its number); refused before anything is read as row R29's
          line refuses, and unless TRIES.md holds exactly one row <ROW> whose first cell is the captain's order
          ("ruling, no run") naming `--entry organ_scope`, `--trend-strength 0.75`, `--trend-strength 1.25` and
          `--captain-ruling <ROW>`. Rows R28's and R29's lines, and every other value's refusal (0.25's included),
          are as before. Its record's trend_strength field names the option, the strength, the file, <ROW>, the
          order text, the row's first cell as TRIES.md holds it and the command (the strength written as given,
          never with one decimal). Since row R52's change <ROW> is pinned: the two are row R49's alone, as 0.5 and
          1.0 are rows R28's and R29's (a row TRIES.md does not hold as their order is refused word for word as
          before, and a row that does but is not R49 is refused after that check). File D's strength 1.5 (row R52
          of experiments/harness/TRIES.md, the captain's order of 2026-11-24, ruling, no run: his named exception on
          the heart's scope and size for this file, H34's stop and the line's cap lifted for this file only; off by
          default): heart --entry organ_scope --trend-strength 1.5 --z-mode raw --captain-ruling R52, file OTS15,
          the same step on D at 1.5, with OTS05's draw, organ-cut pool, count rule, size procedure, seed and half
          (halve.py at seed 1000001); the rows and coordinates are OTS05's. Refused before anything is read with any
          other entry or z mode, with --expression copy, --trend, --depth, --herd and --cfm, on any board but the
          test board (no dry run), without --captain-ruling R52 or with another row, and unless TRIES.md holds
          exactly one row R52 whose first cell is the captain's order ("ruling, no run") naming
          `--entry organ_scope`, `--trend-strength 1.5` and `--captain-ruling R52`. The value is the literal 1.5:
          1.50, 1.5e0 and every other value stay refused as before, and rows R28's, R29's and R49's lines are
          unchanged. Its record's trend_strength field names the option, the strength, the file, R52, the order
          text, the row's first cell as TRIES.md holds it and the command. The trend-form line (row R55 of
          experiments/harness/TRIES.md, the captain's order of 2026-11-24, ruling, no run: his named exception on the
          heart's scope and size for OTK125, OTX125, OTW125 and OTP125, bound to their screen, row H35; off by
          default): heart --entry organ_scope --trend-strength 1.25 --trend-form category | contrast | whole |
          portion --z-mode raw --captain-ruling R55, files OTK125, OTX125, OTW125 and OTP125, OTS125's command with
          one change, the direction of the trend step (vec.trend_forms): the Endothelial lineage split by chamber
          category; each lineage's factor less the whole stage's shared change, its norm held; E10.5's point on the
          whole stage; the slope over the stages before E10.5, its norm held. OTS125's draw, rows, coordinates,
          count rule, size procedure, seed and half (halve.py at seed 1000001). Refused before anything is read
          on every other command, board or seed, without --captain-ruling R55 or with another row, and unless
          TRIES.md holds exactly one row R55 whose first cell is the captain's order naming `--entry organ_scope`,
          `--trend-strength 1.25`, the form's `--trend-form <F>` and `--captain-ruling R55` as whole backticked
          strings. Its record's trend_strength field names the option, the form, the file, R55, the order text, the
          row's first cell, the command, b0_trend's record and each type's status, and no cell number
  embryo  --method pseudobulk_shift | fgw | fgw_shape0 | fgw_shape1: E7.5 cells moved by the E7.25 -> E7.5 change, or
          the FGW interpolation between E7.5 and E8.0 at t = 1/2 (fgw_shape0, the s2 of docs/S2_RULES.md, keeps the
          mixture in E7.5's geometry); sized by the pre-registered embryo procedure (size_embryo_E7.75.json).
          --positions floor_exact (off by default; with --method fgw_shape0 only; tries E47a and E47b of
          experiments/harness/TRIES.md): every cell on a real E7.5 position, the harness's fgw_shape0x, or with
          --expression series its fgw_shape0x_series
          --row-order copy_parity (off by default; no try, no ruling; branch campaign/embryo-occupancy): the same
          cells, expression and positions, the rows re-ordered (vec.t2_build.copy_parity_order) so that the file's
          pre-fix SVD parity on the rows the scorer's _match_n takes is the floor copy's (vec.evaluate.subsample of
          E7.5 as released, 0.1, seed 0, the bake protocol); refused unless the copy and the file get the same truth
          rows (a gate on the test board's truth_cells, read by code; nothing about it is printed). On a test board it
          is refused unless --row-order-test-ok (the captain's yes) is also given: a refusal there, or a success,
          shows one bit of the test count. It assumes the server scores the prediction's rows as the file holds them
          (file order, full count, no subsample or sort); score_h5ad's `_prediction` does this in the pinned copy,
          but the live server's handling is not shown. The as-built uploaded files' pre-fix parities are never
          computed, by any script, and the returned occupancy of a copy_parity upload is never compared with an
          earlier upload's to infer a mode (limit 2)
          --shape-interp stretch_ratio | stretch_registered | overlay (off by default: none; no try read, no ruling;
          row E49 as drafted, branch campaign/embryo-shape-interp): with --method fgw_shape0 --positions floor_exact
          only, one form of vec.shape_interp applied to the exact-position cloud between E7.5 and E8.0 at t = 1/2,
          before the size rule (so the pre-registered size and the count rule are untouched); the record gains a
          shape_interp block. On a test board it is refused unless --shape-interp-test-ok (the captain's yes) is also
          given. With --row-order copy_parity the order is fixed: form, size, jitter, cloud guard, then the row order
          on the final coordinates. Row E55 (as specified, branch campaign/emb-osr) adds --shape-interp overlay_stretch
          (vec.shape_interp.COMPOSED_FORMS: the overlay's union stretched along its own principal axes to the brackets'
          interpolated profile), under the same refusals; its record block gains f1, f2, f3, g, the target SDs and the
          thin SD. Its placebos (overlay_stretch_anti, overlay_stretch_axes) are harness-only and refused by name
          --copy-rows (off by default; the investigator's Option A, the captain's order of 2026-11-21, docs/PLAN.md
          decision log 2026-11-21 (h); row E50 as specified): with --method fgw_shape0 --positions floor_exact and
          --expression mixture or series only, the file takes the floor copy's own rows, order, count and float32
          coordinates (floor_copy_coords) and the expression of vec.t2_build.copy_rows_cloud; no size rule, centring,
          jitter or row-order step. It departs from the count rule (the record says "floor_copy" for n_cells on a test
          board, and its count_rule says departed) and may hold fewer cells than the truth (allow_below_truth). The
          limits are checked here, and every refusal on this path shows no value; the coordinates are checked bit for
          bit before writing and on the file read back. On a test board it is refused unless --copy-rows-test-ok (the
          captain's yes) is also given, before anything loads
          --copy-rows-size (off by default; the copy-rows file at the team's pre-registered size, the captain's order
          of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (v); row E52 as specified): with --copy-rows and
          --expression mixture only. After the floor copy's bit check, the coordinates are scaled in float64 about the
          copy's own centroid to the embryo size JSON's preregistered.pred_rms_full (raw mode only, read as for the
          exact-position file; vec.t2_build.copy_rows_sized) and stored as float32; rows, order, count and expression
          stay the copy-rows file's, and nothing follows but writing (no jitter, centring, cloud guard or row-order
          step). Self-checks, each refusing with no value: the RMS within 1e-6 relative of the target, the prefix
          parity of the matched rows at metric seed 0 equal to the floor copy's for every count from 4 to the copy's,
          and the file read back equal to the float32 coordinates bit for bit (else the file and its stamp are
          removed). On a test board it needs --copy-rows-size-test-ok as well as --copy-rows-test-ok, before anything
          loads
          --expression series_back --captain-ruling R48 (off by default; row R48 of experiments/harness/TRIES.md, the
          captain's order of 2026-11-24, "Build T2embryo_s2x_ovsb.h5ad as named": the embryo's depth line reopened
          for this one file, an extra file beside O and never in its place): O's command (embryo --method fgw_shape0
          --positions floor_exact --shape-interp overlay --shape-interp-test-ok) with the E8.0 rows moved into E7.5's
          run profile and depth (vec.t2_build.series_back_expression, the mirror of the series rendering: SERIES_RULE's
          factor computed from E8.0 to E7.5 on the same matched extraembryonic labels, weights E8.0's label shares,
          applied in count space, then RENDER_RULE's lo thinning to E7.5's median depth) and the E7.5 rows as
          released; O's rows, order, positions, count, scope and size rule unchanged. Refused before anything is read
          on any other command (the test board, the default seed, nothing else changed), without --captain-ruling R48
          or with another row, and unless TRIES.md holds exactly one row R48 whose first cell is the captain's order
          ("ruling, no run") naming `--expression series_back` and `--captain-ruling R48`; embryo --captain-ruling
          without it is refused. Its record gains a series_back block (the option, the file, R48, the order text, the
          command, the rule) and the flat series_back_* fields, no cell number. Its step 0 is row E56
          (experiments/coords/ovsb_step0.py); its test builds, experiments/builds/embryo_ovsb.sh
          --mix-lean later --captain-ruling R59 (off by default; row R59 of experiments/harness/TRIES.md, the captain's
          order of 2026-11-24, "Yes to T2embryo_s2x_ovt23.h5ad (Ot23): my named exception on the embryo's overlay form
          for this file, bound to its gate G1."): O's command with the later bracket's share of the drawn cells at
          t_mix = 2/3 exactly, passed to vec.t2_build.interpolate_clouds (t = 1/2 still sets the moves): a third of the
          rows real E7.5 cells, two thirds real E8.0 cells, so the overlay's union and pool weights lean to E8.0's
          specimen (it moves the overlay's form). Count, size rule, scope, seed and every stream are O's. Refused
          before anything is read on any other command (the test board, the default seed, nothing else changed),
          without --captain-ruling R59 or with another row, and unless TRIES.md holds exactly one row R59 whose first
          cell is the captain's order ("ruling, no run") naming `--mix-lean later` and `--captain-ruling R59`; embryo
          --captain-ruling R59 without it is refused. Its record gains a mix_lean block (the option, the file, R59, the
          order text, the command, t_mix), no cell number. Its gate is G1, row E58
          (experiments/coords/e58_mix_lean_read.py); its test builds, experiments/builds/embryo_ot23.sh
          --neighbourhood run_sorted --captain-ruling <ROW> (off by default: as_built; row E57 of
          experiments/harness/TRIES.md, "RS: one run per side", the file T2embryo_s2x_ovrs.h5ad): O's command with one
          last step on the final coordinates, after the cloud guard (the place --row-order takes), that re-assigns
          expression rows only (vec.run_sort.sort_rows: each bracket's rows moved to its own side of the earlier
          bracket's midline plane, within fixed lineage pools, by linear_sum_assignment on mirror-homologous
          positions); O's coordinates (bytes and order), count, names and the expression rows as a multiset are
          unchanged. Refused before anything is read on any other command (the test board, the default seed, nothing
          else changed), without --captain-ruling, and unless TRIES.md holds exactly one row of the given ID whose
          first cell is the captain's ruling ("ruling, no run") naming `--neighbourhood run_sorted` and
          `--captain-ruling <ROW>` (_captain_ruling_row). Its record gains a neighbourhood block (the option, the
          file, the row, the command, the rule, the plane's normal and contrast), no cell number. Its practice read is
          row E57 (experiments/coords/e57_rs_read.py); its test builds, experiments/builds/embryo_rs.sh

The board is the test board from boards/ (after the intake, `vec boards update --from-file`); before that, --board
may name a validation board of the same task for a dry run, and the script refuses a dry run whose output file name
does not contain "dry" (such a file must never be uploaded). The build seed is never 0 (runbook section 5): default
vec.harness.BUILD_SEED_OFFSET + 1. T2 counts go through approved rule 7's committed decision
(`t2_build.n_cells_from_decision`, runbook section 5), checked against the count the harness's builder draws. On a
test board no cell count is printed or stored: the record says "max_cells" when the file holds the board's max_cells
(the count rule) and refuses otherwise (red-team review 2026-10-31; s1 build review 2026-11-01).

Usage (DGX, repo root, BLAS threads capped):
  uv run python experiments/builds/build_test.py t1 --out runs/builds/T1_s1.h5ad
  uv run python experiments/builds/build_test.py heart --entry s1 --trend --out runs/builds/T2heart_s1.h5ad
  uv run python experiments/builds/build_test.py embryo --method fgw --out runs/builds/T2embryo_s1.h5ad
  uv run python experiments/builds/build_test.py t1 --cm-share 0.60 --out runs/builds/T1_s2.h5ad
The commands of the entries as they stand are in experiments/builds/build_all.sh and docs/S2_RULES.md (since
2026-11-13 T1's two files are `t1 --g-mom --base lib3` and `t1 --g-mom --base mid --cm-share 0.60`, and both heart
files carry --trend: the captain's ruling of 2026-11-13, docs/PLAN.md decision log (w); row H11).
  uv run python experiments/builds/build_test.py heart --entry s2_trim --trend --out runs/builds/T2heart_s2.h5ad
  uv run python experiments/builds/build_test.py embryo --method fgw_shape0 --out runs/builds/T2embryo_s2.h5ad
Without --trend the heart builder gives the heart file as it was before that ruling (B0's E9.5 -> E10.5 step alone),
kept for reproducing the earlier reference set (builds_38616fe).

Rows R54 and R58 of experiments/harness/TRIES.md (the step-change round of 2026-11-24/25, T1 lines LM16 and m20):
`t1 --g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --late-mix e15_17
--captain-ruling R54` (T1_s2rcMSB_LM16) is T1_s2rcMSB's command with library 3's own chamber mix moved to the shares
row N64's step commits (experiments/t1_scope/late_mix.py; cardiomyocytes of heart cells, then blood of chamber cells,
by vec.t1_scope.with_share), taken only while row R54's first cell is the captain's ruling naming `--late-mix e15_17`
and `--captain-ruling R54` and the shares file committed at HEAD has the sha256 row N64 names; `t1 --g-mom --base lib3
--feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 2.0
--captain-ruling R58` (T1_s2rcMSBm20) is T1_s2rcMSB's command with MS's own damping 2.0, taken only while row R58's
first cell is the captain's order naming that option string. Both routes are tested before row R51's; without
--late-mix or --captain-ruling R58 every command is taken or refused exactly as before.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

from vec import harness, provenance, run_sort, shape_interp, t1_bets, t1_scope, t2_build
from vec.boards import REPO_ROOT, load_boards
from vec.evaluate import SUBSAMPLE_FRAC, Stage, subsample
from vec.io import COORD_KEY, dense, read_stage
from vec.prep import load_stage

sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_m1"))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_hygiene"))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_strong_null"))
import feng_dirs  # noqa: E402
import hygiene  # noqa: E402
import strong_null  # noqa: E402

TEST_BOARDS = {"t1": "T1:test", "heart": "T2:heart:test_extrap", "embryo": "T2:embryo:test_interp"}
DRY_BOARDS = {"t1": {"T1:val"}, "heart": {"T2:heart:val_extrap", "T2:heart:val_interp"},
              "embryo": {"T2:embryo:val_interp"}}
COUNT_DECISION = REPO_ROOT / "experiments/density/results/release/count_decision.json"
COUNT_FIELDS = ("n", "n_before", "n_after", "n_jittered")  # the interpolation record's counts
SIZES = {"heart": REPO_ROOT / "experiments/coords/results/p3/size_heart_E12.5.json",
         "embryo": REPO_ROOT / "experiments/coords/results/p3/size_embryo_E7.75.json"}
DATA = harness.DATA_DIR
# fgw and its s2 shape variants (docs/S2_RULES.md): the position fraction passed to interpolate_clouds
EMBRYO_T_POS = {"fgw": None, "fgw_shape0": 0.0, "fgw_shape1": 1.0}
# R3b of experiments/harness/TRIES.md (the captain's ruling of 2026-11-10 (p)): the only categories g_mom's factors edit
# at test. These are the two any backtest read (T25, T26); epicardium and vascular endothelium have factors but no read.
G_MOM_EDIT_CATEGORIES = ("cm", "endocardium")
# --feng-groups (off by default): the cell groups Feng's direction is formed on. "authors" is Feng's own labels, every
# file as built. "none" is no Feng direction on any group: the entry with that one step taken out and nothing put in
# its place: no Feng file is read, and neither is any list that was fixed by counting Feng's cells (docs/PLAN.md
# decision log 2026-11-19 (a) to (e): the authors' labels were made with every stage pooled, the in-window ones
# included, which an organiser's condition of 2026-11-16 excludes). Only the captain's ruling puts it in an upload.
# "recomputed" is the first upload's design on the groups experiments/external/feng_relabel.py recomputes from the
# permitted cells' raw counts alone (register rows X002 and X010; the captain's "utilize the atlases", decision log
# 2026-11-19 (k)): it builds one file, the first upload's command, and only when the labeller's committed record routes
# to it. Only the captain's ruling puts it in an upload.
FENG_GROUPS = feng_dirs.FENG_GROUPS
FENG_GROUPS_DEFAULT = FENG_GROUPS[0]
# R6 of experiments/harness/TRIES.md (2026-11-13: the AI assistant's decision under the captain's delegation,
# docs/PLAN.md decision log 2026-11-13 (d)): the pool T1's cells are drawn from. "pooled": every E10.5 cell (the entry
# until then, and still T1's second upload). "lib3": E10.5's dissected-heart library alone, whose cell names end in
# this suffix, at that library's own chamber mix (tries T27 to T31). Feng's direction, R1, R2 and the fit of g_mom
# are untouched.
# R7 (the captain's ruling of 2026-11-13, docs/PLAN.md decision log 2026-11-13 (k); try T33): "mid", the base half-way
# between the two, for T1's second upload: every library-3 cell and, per chamber category, a seeded random fraction
# s / (1 + s) of the other libraries' cells (s: library 3's share of that category's cells in the file), so that
# library 3 supplies (s + 1) / 2 of each category, at that pool's own chamber mix.
T1_BASES = ("pooled", "lib3", "mid")
T1_DISSECTED_LIBRARY = "3"
T1_MID_STREAM = 21  # the mid base's draw of other-library cells: default_rng([seed, 21]), apart from every build stream
T1_MID_CM_SHARE = 0.60  # the one cardiomyocyte share the mid base is ruled for (docs/S2_RULES.md: T1's s2)
# b0_trend's earlier input stages for the heart (--trend; row H11, the captain's ruling of 2026-11-13, docs/PLAN.md
# decision log (w)): role -> (file, day as vec.harness.STAGE_DAYS gives it); with E9.5 ("start") and E10.5 ("ref")
HEART_TREND_STAGES = {role: (f, harness.STAGE_DAYS[f]) for role, f in
                      (("e8_75", "T2/E8.75.h5ad"), ("e8_5", "T2/E8.5.h5ad"), ("e8_25_late", "T2/E8.25_late.h5ad"))}
# The heart's scope bet (campaign/heart-iso-trim; no ruling; only the captain's ruling puts it in an upload): the ruled
# second file's two fixes (the copy's expression, the size step's raw mode) on the pre-registered trimmed scope, H2 as
# the half trim of docs/S2_RULES.md, sized by the size file's own logged raw-mode figure for H2. Scope is a lever of the
# post-exposure block (docs/DATA_POLICY.md 2026-11-04): the option is off unless all three flags are given.
HEART_ISO_TRIM = "heart --entry s2_trim --expression copy --z-mode raw"
# The full trim (campaign/heart-scope-options; no ruling; off unless all three flags are given, only the captain's
# ruling puts it in an upload): the same two fixes on E10.5 trimmed by vec.t2_build.radial_trim_mask at fraction 0
# (surrounding cells removed farthest in plane first until the in-plane scope is at most zero: the core and any
# surrounding cells inside its in-plane spread), sized to the raw mode's own core forecast
# (z_modes.raw.pred_rms_core of the size JSON: the cloud at the core's own scope, full-over-core ratio one).
HEART_CORE_TRIM = "heart --entry core_trim --expression copy --z-mode raw"
# Per scope option: its trim fraction, the size step's hypothesis, the size JSON field its size comes from
HEART_SCOPE_OPTIONS = {
    "s2_trim": {"frac": 0.5, "hypothesis": "H2", "size_source": "z_modes.raw.scope.H2.pred_rms_full",
                "rule": "surrounding cells removed farthest in plane from the core first, to half the in-plane scope"},
    "core_trim": {"frac": 0.0, "hypothesis": t2_build.SIZE_CORE_SCOPE, "size_source": "z_modes.raw.pred_rms_core",
                  "rule": "surrounding cells removed farthest in plane from the core first, to an in-plane scope of "
                          "at most zero (the core's own)"},
}
# The organ cut, file D of the heart's order rule (row R21 of experiments/harness/TRIES.md, the captain's order of
# 2026-11-22 to code and read it; docs/S2_RULES.md, amendment of 2026-11-22; off unless all three flags are given, and
# only the captain's yes by name, after its gate and its read, puts it in an upload): the copy draw from E10.5 under
# vec.t2_build.organ_scope_mask (the chamber categories of configs/t2_chamber_classes.yaml joined to the core), sized
# to the raw-mode core forecast times rms(organ pool) / rms(core) on E10.5 (vec.t2_build.size_organ_scope); the count
# rule, the seed and the halving are the full trim's.
HEART_ORGAN_SCOPE = "heart --entry organ_scope --expression copy --z-mode raw"
HEART_ORGAN_RULING = ("row R21 of experiments/harness/TRIES.md: the captain's order to code and read file D; an "
                      "upload needs its gate, its read, two identical DGX builds and his yes by name")
# The heart search line (row R28 of experiments/harness/TRIES.md, the captain's order of 2026-11-23, "start the heart
# search"; docs/PLAN.md decision log 2026-11-23 (t)): its two files TS05 and TS10 are the full trim C's command with the
# trend step at strength 0.5 and 1.0 and nothing else changed. The strength is row H28's: the b0_trend factors
# (vec.harness.b0_trend_factors) and momentum_shift's fallback damping both multiplied by it, a test budget of 0.5 x
# the strength day, applied before the size step on C's own draw, pool, count rule, size field and seed. Off by
# default; taken only with --entry core_trim --z-mode raw --captain-ruling R28, and only while row R28 is the captain's
# order naming both options and the ruling.
HEART_TS_ROW = "R28"
HEART_TS_OPTIONS = {0.5: "`--trend-strength 0.5`", 1.0: "`--trend-strength 1.0`"}
HEART_TS_CHOICES = ("0.5", "1.0")  # the option's literal values, as row R28 names them
HEART_TS_FILES = {0.5: "TS05", 1.0: "TS10"}
HEART_TS_COMMANDS = {s: f"heart --entry core_trim --trend-strength {s:.1f} --z-mode raw --captain-ruling {HEART_TS_ROW}"
                     for s in HEART_TS_OPTIONS}
HEART_TS_ORDER = ("experiments/harness/TRIES.md row R28, the captain's order of 2026-11-23 "
                  "(\"start the heart search\"; docs/PLAN.md decision log 2026-11-23 (t)): the heart search line's "
                  "two files, TS05 and TS10, the full trim C's command with the trend step at strength 0.5 and 1.0 "
                  "(row H28's strength), nothing else changed")
HEART_TS_STEP = ("row H28's strength-scaled b0_trend step: the vec.harness.b0_trend_factors mapping and "
                 "momentum_shift's fallback damping (vec.harness.B0_DAMPING) both multiplied by the strength, before "
                 "the size step, on C's draw, pool, count rule, size field and seed")
# The heart search line on file D (row R29 of experiments/harness/TRIES.md, the captain's order "i want it too": the
# trend tried on D as well; recorded in docs/PLAN.md decision log 2026-11-23 (w)): its two files OTS05 and OTS10 are the
# organ cut D's command with --expression copy replaced by the same strength-scaled step (_trend_strength_step), before
# the size step on D's own draw, organ-cut pool, count rule, size procedure and seed; the rows, coordinates, cell
# number, size field and half's rows stay D's. Off by default; taken only with --entry organ_scope --z-mode raw
# --captain-ruling R29, and only while TRIES.md holds exactly one row R29 whose first cell is the captain's order
# ("ruling, no run") naming `--entry organ_scope`, both options and `--captain-ruling R29`. Row R28's line is unchanged:
# core_trim with R29 and organ_scope with R28 are both refused.
HEART_OTS_ROW = "R29"
HEART_OTS_ENTRY = "`--entry organ_scope`"
HEART_OTS_FILES = {0.5: "OTS05", 1.0: "OTS10"}
HEART_OTS_COMMANDS = {s: f"heart --entry organ_scope --trend-strength {s:.1f} --z-mode raw --captain-ruling "
                         f"{HEART_OTS_ROW}" for s in HEART_TS_OPTIONS}
HEART_OTS_ORDER = ("experiments/harness/TRIES.md row R29, the captain's order (\"i want it too\"; docs/PLAN.md "
                   "decision log 2026-11-23 (w)): the trend step tried on file D too, two files, OTS05 and OTS10, the "
                   "organ cut D's command with --expression copy replaced by the trend step at strength 0.5 and 1.0 "
                   "(row H28's strength), nothing else changed")
HEART_OTS_STEP = ("row H28's strength-scaled b0_trend step: the vec.harness.b0_trend_factors mapping and "
                  "momentum_shift's fallback damping (vec.harness.B0_DAMPING) both multiplied by the strength, before "
                  "the size step, on D's draw, organ-cut pool, count rule, size procedure and seed")
# File D's two further strengths (the design round of 2026-11-24, its FAST proposal for the heart with 0.25 dropped by
# the round's challenge; its ruling row, "R-F" there, is the row given on the command line): files OTS075 and OTS125,
# the organ cut D's command with --expression copy replaced by the same strength-scaled step (_trend_strength_step) at
# 0.75 and 1.25, before the size step on OTS05's draw, organ-cut pool, count rule, size procedure and seed (D's own);
# then halve.py at seed 1000001 as for D. Off by default; taken only with --entry organ_scope --z-mode raw
# --captain-ruling <ROW>, <ROW> a ruling row's id (R and its number), and only while TRIES.md holds exactly one row
# <ROW> whose first cell is the captain's order ("ruling, no run") naming `--entry organ_scope`, the two options and
# `--captain-ruling <ROW>` (_captain_order_names). Rows R28's and R29's lines are unchanged: 0.5 and 1.0 stay theirs,
# and every other value, 0.25 included, stays refused as before. The commands, records and file names carry the
# option's literal values (the R28 and R29 tables' one-decimal format would write 0.75 and 1.25 as 0.8 and 1.2).
HEART_FURTHER_CHOICES = ("0.75", "1.25")  # the option's literal values, as the ruling row names them
HEART_FURTHER_OPTIONS = {float(v): f"`--trend-strength {v}`" for v in HEART_FURTHER_CHOICES}
HEART_FURTHER_FILES = {0.75: "OTS075", 1.25: "OTS125"}
HEART_FURTHER_ROW_FORM = re.compile(r"R[0-9]+[a-z]?")  # a ruling row's id, e.g. R44
HEART_FURTHER_COMMANDS = {float(v): f"heart --entry organ_scope --trend-strength {v} --z-mode raw --captain-ruling "
                                    "{row}" for v in HEART_FURTHER_CHOICES}  # .format(row=<ROW>)
HEART_FURTHER_ORDER = ("experiments/harness/TRIES.md row {row}, the captain's order (\"ruling, no run\") naming "
                       "`--entry organ_scope`, `--trend-strength 0.75`, `--trend-strength 1.25` and "
                       "`--captain-ruling {row}` (the design round of 2026-11-24, its FAST proposal for the heart "
                       "with 0.25 dropped by the round's challenge): file D's two further strengths, files OTS075 and "
                       "OTS125, the organ cut D's command with --expression copy replaced by the trend step at "
                       "strength 0.75 and 1.25 (row H28's strength), nothing else changed")
HEART_FURTHER_STEP = ("row H28's strength-scaled b0_trend step: the vec.harness.b0_trend_factors mapping and "
                      "momentum_shift's fallback damping (vec.harness.B0_DAMPING) both multiplied by the strength, "
                      "before the size step, on D's draw, organ-cut pool, count rule, size procedure and seed "
                      "(OTS05's and OTS10's)")
# The two further strengths' ruling row, pinned since row R52's change (campaign/heart-ots15) as rows R28's and R29's
# are: row R49 of experiments/harness/TRIES.md, the captain's order of 2026-11-24 naming them. The row given on the
# command line is still checked as before (HEART_FURTHER_ROW_FORM, then _heart_further_ordered), so every refusal the
# line gave stands word for word; a row that passes that check but is not R49 is then refused too.
HEART_FURTHER_ROW = "R49"
# File D's strength 1.5, file OTS15 (row R52 of experiments/harness/TRIES.md, the captain's order of 2026-11-24, ruling,
# no run: "Yes to OTS15: my named exception on the heart's scope and size for this file, and H34's stop and the line's
# cap lifted for this file only."; decision log 2026-11-24 (kkk)): the organ cut D's command with --expression copy
# replaced by the same strength-scaled step (_trend_strength_step) at 1.5, before the size step on OTS05's draw,
# organ-cut pool, count rule, size procedure and seed (D's own); then halve.py at seed 1000001 as for D. Off by default;
# taken only with --entry organ_scope --z-mode raw --captain-ruling R52, and only while TRIES.md holds exactly one row
# R52 whose first cell is the captain's order ("ruling, no run") naming `--entry organ_scope`, `--trend-strength 1.5`
# and `--captain-ruling R52` (_heart_ots15_ordered). The value is the option's literal 1.5 (1.50, 1.5e0 and the like
# stay refused as before), and rows R28's, R29's and R49's lines are unchanged.
HEART_OTS15_ROW = "R52"
HEART_OTS15_CHOICE = "1.5"  # the option's literal value, as row R52 names it
HEART_OTS15_STRENGTH = float(HEART_OTS15_CHOICE)
HEART_OTS15_OPTION = f"`--trend-strength {HEART_OTS15_CHOICE}`"
HEART_OTS15_FILE = "OTS15"
HEART_OTS15_COMMAND = (f"heart --entry organ_scope --trend-strength {HEART_OTS15_CHOICE} --z-mode raw "
                       f"--captain-ruling {HEART_OTS15_ROW}")
HEART_OTS15_ORDER = ("experiments/harness/TRIES.md row R52, the captain's order of 2026-11-24 (ruling, no run: "
                     "\"Yes to OTS15: my named exception on the heart's scope and size for this file, and H34's stop "
                     "and the line's cap lifted for this file only.\"; decision log 2026-11-24 (kkk)): file D's "
                     "strength 1.5, file OTS15, the organ cut D's command with --expression copy replaced by the "
                     "trend step at strength 1.5 (row H28's strength), nothing else changed")
HEART_OTS15_STEP = ("row H28's strength-scaled b0_trend step: the vec.harness.b0_trend_factors mapping and "
                    "momentum_shift's fallback damping (vec.harness.B0_DAMPING) both multiplied by the strength, "
                    "before the size step, on D's draw, organ-cut pool, count rule, size procedure and seed "
                    "(OTS05's, OTS10's, OTS075's and OTS125's)")
# --shape-interp (row E49 as drafted; no ruling): off by default; the fixed order of the steps when it is given
SHAPE_INTERP_OFF = "none"
OSR_TEST_BLOCKS = "1|2,3"  # row E55: overlay_stretch on a test board only with these bracket blocks (seed-free)
SHAPE_INTERP_ORDER = ("floor_exact positions, the shape form, the size rule's RMS rescale, the duplicate jitter, the "
                      "cloud guard, then --row-order copy_parity if given (on the final coordinates)")
# --copy-rows (off by default; the investigator's Option A, the captain's order of 2026-11-21, docs/PLAN.md decision log
# 2026-11-21 (h); row E50 as specified): the floor copy's own rows, order, count and float32 coordinates with the
# exact-position file's expression rule (vec.t2_build.copy_rows_cloud)
COPY_ROWS_EXPRESSIONS = ("mixture", "series")
COPY_ROWS_N = "_copy_rows_n"  # the floor copy's count, handed to _hide_counts and never written
COPY_ROWS_COUNT_RULE = "departed: the floor copy's count (the captain's order of 2026-11-21)"
# --copy-rows-size (off by default; row E52 as specified, the captain's order of 2026-11-21, docs/PLAN.md decision log
# 2026-11-21 (v)): the copy-rows file (mixture only) rescaled about the copy's centroid to the pre-registered size
COPY_ROWS_SIZE_EXPRESSION = "mixture"
COPY_ROWS_SIZE_RULE = "applied: rescale_to_rms about the copy centroid to preregistered pred_rms_full (raw mode)"
COPY_ROWS_SIZE_RMS_TOL = 1e-6  # the RMS self-check, relative
COPY_ROWS_PARITY_FROM = 4  # the prefix parity self-check's smallest count (three rows have rank two: arbitrary)
COPY_ROWS_PARITY_SEED = 0  # its metric seed
# --expression series_back (row R48 of experiments/harness/TRIES.md, the captain's order of 2026-11-24, "Build
# T2embryo_s2x_ovsb.h5ad as named": the embryo's depth line reopened for this one file, an extra file beside O and never
# in its place): O's command with the E8.0 rows moved into E7.5's run profile and depth (vec.t2_build.SERIES_BACK_RULE,
# the mirror of the series rendering); O's rows, order, positions, count, scope and size rule are unchanged. Off by
# default; taken only on O's command with --captain-ruling R48 (the test board, the default seed, nothing else changed),
# and only while TRIES.md holds exactly one row R48 whose first cell is the captain's order ("ruling, no run") naming
# `--expression series_back` and `--captain-ruling R48`. Its step 0 is row E56 (experiments/coords/ovsb_step0.py).
EMBRYO_OVSB_ROW = "R48"
EMBRYO_OVSB_OPTION = "`--expression series_back`"
EMBRYO_OVSB_FILE = "T2embryo_s2x_ovsb.h5ad"
EMBRYO_O_COMMAND = "embryo --method fgw_shape0 --positions floor_exact --shape-interp overlay --shape-interp-test-ok"
EMBRYO_OVSB_COMMAND = ("embryo --method fgw_shape0 --positions floor_exact --expression series_back --shape-interp "
                       f"overlay --shape-interp-test-ok --captain-ruling {EMBRYO_OVSB_ROW}")
EMBRYO_OVSB_ORDER = ("experiments/harness/TRIES.md row R48, the captain's order of 2026-11-24 (\"Build "
                     "T2embryo_s2x_ovsb.h5ad as named\"; decision log 2026-11-24 (tt)): the embryo's depth line "
                     "reopened for this one file, an extra file beside O and never in its place; O's command with "
                     "--expression series_back, taken only with --captain-ruling R48")
EMBRYO_OVSB_STEP0 = "row E56 of experiments/harness/TRIES.md (experiments/coords/ovsb_step0.py)"
# --mix-lean later (row R59 of experiments/harness/TRIES.md, the captain's order of 2026-11-24, "Yes to
# T2embryo_s2x_ovt23.h5ad (Ot23): my named exception on the embryo's overlay form for this file, bound to its gate
# G1."): O's command with the later bracket's share of the drawn cells at t_mix = 2/3 exactly, passed to
# t2_build.interpolate_clouds (t = 1/2 still sets the moves). It moves the overlay's form: the union and its pool
# weights lean to E8.0's specimen. Off by default; taken only on O's command with --captain-ruling R59 (the test board,
# the default seed, nothing else changed), and only while TRIES.md holds exactly one row R59 whose first cell is the
# captain's order ("ruling, no run") naming `--mix-lean later` and `--captain-ruling R59`. Its gate is G1, row E58
# (experiments/coords/e58_mix_lean_read.py: the harness arm takes t_mix = min(1, t + 1/6), which is 2/3 at t = 1/2).
EMBRYO_OT23_ROW = "R59"
EMBRYO_OT23_OPTION = "`--mix-lean later`"
EMBRYO_OT23_FILE = "T2embryo_s2x_ovt23.h5ad"
MIX_LEANS = ("later",)  # the earlier lean (Ot13) was struck by the cold challenge: no value for it
EMBRYO_OT23_T_MIX = 2.0 / 3.0  # the nearest double to 2/3, the same float as min(1, 1/2 + 1/6)
EMBRYO_OT23_COMMAND = f"{EMBRYO_O_COMMAND} --mix-lean later --captain-ruling {EMBRYO_OT23_ROW}"
EMBRYO_OT23_ORDER = ("experiments/harness/TRIES.md row R59, the captain's order of 2026-11-24 (\"Yes to "
                     "T2embryo_s2x_ovt23.h5ad (Ot23): my named exception on the embryo's overlay form for this file, "
                     "bound to its gate G1.\"; decision log 2026-11-24 (www)): O's command with --mix-lean later, "
                     "taken only with --captain-ruling R59")
EMBRYO_OT23_GATE = "G1, row E58 of experiments/harness/TRIES.md (experiments/coords/e58_mix_lean_read.py)"
EMBRYO_OT23_RULE = ("t_mix = 2/3 in vec.t2_build.interpolate_clouds at t = 1/2: round(n / 3) rows real E7.5 cells, the "
                    "rest real E8.0 cells; the overlay's union and pool weights follow t_mix; count, size rule, scope, "
                    "seed and every stream as O's")
# Row R75 (the captain's order of 2026-11-26, ruling, no run, decision log 2026-11-26 (a)): the embryo's refinement
# round on O, five files sent at once with no gate, each O's command (EMBRYO_O_COMMAND) with exactly one option string
# of the table below, each taken only with --captain-ruling R75 (the test board, the default seed, nothing else
# changed), and only while TRIES.md holds exactly one row R75 whose first cell is the captain's order ("ruling, no
# run") naming the string and `--captain-ruling R75`, each as a whole backticked string, the string followed there by
# " (<file>)" as the row pairs each string with its file (_r75_row_names, row R73's check). Row R59's route: a lean
# passes its t_mix to t2_build.interpolate_clouds (t = 1/2 still sets the moves), which takes it as float(t_mix), so
# each exact fraction goes in as its nearest double (2/3 is EMBRYO_OT23_T_MIX's float); --expression series and depth
# are the builder's existing expressions on O's rows and positions. _r75_refusals runs first in build_embryo, before
# the older rows' refusals read the ruling; every other command is what it was. This table, in the row's order, is the
# one the row's first cell is checked against and the one experiments/builds/embryo_r75_line.sh's literal case table
# copies (tests/test_embryo_r75_launcher.py holds the two equal).
EMBRYO_R75_ROW = "R75"
EMBRYO_R75_TABLE = (
    ("T2embryo_s2x_ovt712.h5ad", "--mix-lean later_half"),
    ("T2embryo_s2x_ovt56.h5ad", "--mix-lean later2"),
    ("T2embryo_s2x_ovt512.h5ad", "--mix-lean earlier_half"),
    ("T2embryo_s2xs_ovt23.h5ad", "--expression series --mix-lean later"),
    ("T2embryo_s2x_ovd.h5ad", "--expression depth"),
)
# each lean's t_mix, exact: the later bracket's share of the drawn cells (row R75: 7/12, 5/6, 5/12 and 2/3)
EMBRYO_R75_LEANS = {"later_half": Fraction(7, 12), "later2": Fraction(5, 6), "earlier_half": Fraction(5, 12),
                    "later": Fraction(2, 3)}
EMBRYO_R75_NEW_LEANS = tuple(lean for lean in EMBRYO_R75_LEANS if lean not in MIX_LEANS)  # row R75's own values
EMBRYO_R75_EXPRESSIONS = ("series", "depth")  # the expressions a string of the table may name; O's is mixture
EMBRYO_R75_DEPTH = "depth"  # with --shape-interp overlay, row R75's file alone (the embryo's depth line, reopened)


def _r75_values(name: str, option: str) -> tuple:
    """(expression, lean or None, t_mix as an exact Fraction or None) of one entry of row R75's table, read off its
    string. Raises ValueError unless the string is one of the row's three forms, single-spaced: --mix-lean LEAN,
    --expression series --mix-lean LEAN (LEAN a lean of EMBRYO_R75_LEANS) or --expression depth (no lean), and the name
    is O's (T2embryo_s2x_ov.h5ad) with the expression's mark and the lean's fraction: s after s2x for series, d after ov
    for depth, t and the fraction's numerator and denominator after ov for a lean (T2embryo_s2x_ovt712.h5ad for 7/12,
    T2embryo_s2xs_ovt23.h5ad for series at 2/3)."""
    words = option.split()
    expression, lean = "mixture", None
    if words[:1] == ["--expression"]:
        if len(words) < 2 or words[1] not in EMBRYO_R75_EXPRESSIONS:
            raise ValueError(f"not a row {EMBRYO_R75_ROW} option string: {option}")
        expression, words = words[1], words[2:]
    if words[:1] == ["--mix-lean"]:
        if len(words) != 2 or words[1] not in EMBRYO_R75_LEANS:
            raise ValueError(f"not a row {EMBRYO_R75_ROW} option string: {option}")
        lean, words = words[1], []
    # a lean on every form but depth's, which takes none (no depth with a lean, no series or mixture without one)
    if words or (lean is None) != (expression == EMBRYO_R75_DEPTH) or option != " ".join(option.split()):
        raise ValueError(f"not a row {EMBRYO_R75_ROW} option string: {option}")
    t_mix = None if lean is None else EMBRYO_R75_LEANS[lean]
    want = ("T2embryo_s2x" + ("s" if expression == "series" else "") + "_ov"
            + ("" if t_mix is None else f"t{t_mix.numerator}{t_mix.denominator}")
            + ("d" if expression == EMBRYO_R75_DEPTH else "") + ".h5ad")
    if name != want:
        raise ValueError(f"row {EMBRYO_R75_ROW}'s file {name} does not match its option string {option}")
    return expression, lean, t_mix


def _r75_rule(expression: str, t_mix: Fraction | None) -> str:
    """What row R75's file of this expression and t_mix is, in its record's words."""
    if t_mix is None:
        return ("O's rows, order and positions (interpolate_clouds called as O's, t_mix unset), each bracket's drawn "
                "cells thinned to the predicted depth by vec.t2_build.depth_matched_expression at t = 1/2, stream "
                "[seed, 7]: " + t2_build.DEPTH_RULE)
    rows = (f"t_mix = {t_mix.numerator}/{t_mix.denominator} in vec.t2_build.interpolate_clouds at t = 1/2 (its "
            "nearest double): int(round((1 - t_mix) n)) rows real E7.5 cells, the rest real E8.0 cells; the overlay's "
            "union and pool weights follow t_mix; ")
    expr = ("every row's expression as released (O's mixture)" if expression == "mixture" else
            "the expression SO's: vec.t2_build.series_rendered_expression at t = 1/2, stream [seed, 7], the E7.5 rows "
            "rendered into E8.0's run profile and depth (SERIES_RULE)")
    return rows + expr + "; count, size rule, scope, seed and every other stream as O's"


EMBRYO_R75_FILES = {name: (option, *_r75_values(name, option)) for name, option in EMBRYO_R75_TABLE}
EMBRYO_R75_COMMAND = EMBRYO_O_COMMAND + " {option} --captain-ruling " + EMBRYO_R75_ROW
EMBRYO_R75_ORDER = ("experiments/harness/TRIES.md row R75, the captain's order of 2026-11-26 (\"Yes to "
                    "T2embryo_s2x_ovt712.h5ad, T2embryo_s2x_ovt56.h5ad, T2embryo_s2x_ovt512.h5ad, "
                    "T2embryo_s2xs_ovt23.h5ad and T2embryo_s2x_ovd.h5ad under R22, sent at once in one ruling row, "
                    "with no gate\"; \"ruling, no run\"; decision log 2026-11-26 (a), the embryo's refinement round on "
                    "O) naming `{option}` ({file}) and `--captain-ruling R75`: file {file}, O's command with that one "
                    "option string, taken only with --captain-ruling R75; {what}")
EMBRYO_R75_OVERLAY = "the captain's named exception on the embryo's overlay form and occupancy term for this file"
EMBRYO_R75_WHAT = {
    "T2embryo_s2x_ovt712.h5ad": EMBRYO_R75_OVERLAY,
    "T2embryo_s2x_ovt56.h5ad": EMBRYO_R75_OVERLAY,
    "T2embryo_s2x_ovt512.h5ad": (EMBRYO_R75_OVERLAY + ", and the lean-early mixing line declined on 2026-11-06 (e) "
                                 "reopened for this file alone"),
    "T2embryo_s2xs_ovt23.h5ad": EMBRYO_R75_OVERLAY,
    "T2embryo_s2x_ovd.h5ad": "the embryo's depth line reopened for this file alone, beside O and never in its place",
}
EMBRYO_R75_RULES = {name: _r75_rule(expression, t_mix) for name, (_, expression, _, t_mix) in EMBRYO_R75_FILES.items()}
# Row R78 (the captain's order of 2026-11-26, "ruling, no run"; decision log 2026-11-26 (h)): the embryo's fourth
# refinement round on O (EMBRYO_REF4_design.md with EMBRYO_REF4_challenge.md's fixes), both files sent at once with no
# gate. Row R75's route with this one literal table, in the row's order: O's command with exactly one option string
# more, each only with --captain-ruling R78. later_11_12 leans the draw to t_mix 11/12 (both brackets' tissue kept:
# int(round(n / 12)) rows real E7.5 cells); earlier, the lean-early line declined on 2026-11-06 (e) and reopened for
# its file alone, to 1/3. Depth moves only through the shares; the E8.0 twin (t_mix 1) stays declined and has no lean;
# the design's second file (series at 5/6) is struck as a re-draw in substance and is not built.
EMBRYO_R78_ROW = "R78"
EMBRYO_R78_TABLE = (
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean later_11_12"),
    ("T2embryo_s2x_ovt13.h5ad", "--mix-lean earlier"),
)
EMBRYO_R78_LEANS = {"later_11_12": Fraction(11, 12), "earlier": Fraction(1, 3)}  # t_mix as exact fractions
# both leans are new, neither row R59's nor row R75's: the parser takes them, and row R78's check refuses them under
# every ruling but R78 before anything is read
EMBRYO_R78_NEW_LEANS = tuple(lean for lean in EMBRYO_R78_LEANS
                             if lean not in MIX_LEANS and lean not in EMBRYO_R75_LEANS)


def _r78_values(name: str, option: str) -> tuple:
    """Row R78's pair (file name, option string) read by the row's form: the option string exactly `--mix-lean <lean>`
    with one of EMBRYO_R78_LEANS (single spaces, nothing before, after or between), the name O's with that lean's exact
    t_mix as its mark, T2embryo_s2x_ovt<numerator><denominator>.h5ad (row R75's naming rule). Returns the expression
    (O's mixture), the lean and the exact t_mix; raises ValueError for anything else."""
    words = option.split()
    if len(words) != 2 or words[0] != "--mix-lean" or words[1] not in EMBRYO_R78_LEANS or option != " ".join(words):
        raise ValueError(f"not a row {EMBRYO_R78_ROW} option string: {option}")
    t_mix = EMBRYO_R78_LEANS[words[1]]
    if name != f"T2embryo_s2x_ovt{t_mix.numerator}{t_mix.denominator}.h5ad":
        raise ValueError(f"row {EMBRYO_R78_ROW}'s file {name} does not match its option string {option}")
    return "mixture", words[1], t_mix


EMBRYO_R78_FILES = {name: (option, *_r78_values(name, option)) for name, option in EMBRYO_R78_TABLE}
EMBRYO_R78_COMMAND = EMBRYO_O_COMMAND + " {option} --captain-ruling " + EMBRYO_R78_ROW
EMBRYO_R78_WORDS = {  # the captain's words for each file, as row R78's first cell records them
    "T2embryo_s2x_ovt1112.h5ad": ("Yes to T2embryo_s2x_ovt1112.h5ad under R22, sent at once in one ruling row, with "
                                  "no gate"),
    "T2embryo_s2x_ovt13.h5ad": "Yes to T2embryo_s2x_ovt13.h5ad under R22 in the same row, sent at once, with no gate",
}
EMBRYO_R78_WHAT = {
    "T2embryo_s2x_ovt1112.h5ad": EMBRYO_R75_OVERLAY,
    "T2embryo_s2x_ovt13.h5ad": (EMBRYO_R75_OVERLAY + ", and the lean-early mixing line declined on 2026-11-06 (e) "
                                "reopened for this file alone"),
}
EMBRYO_R78_ORDER = ("experiments/harness/TRIES.md row R78, the captain's order of 2026-11-26 (\"{words}\"; \"ruling, "
                    "no run\"; decision log 2026-11-26 (h), the embryo's fourth refinement round on O) naming "
                    "`{option}` ({file}) and `--captain-ruling R78`: file {file}, O's command with that one option "
                    "string, taken only with --captain-ruling R78; {what}; both brackets' tissue kept, depth moved "
                    "only through the shares, the E8.0 twin (t_mix 1) still declined")
EMBRYO_R78_RULES = {name: _r75_rule(expression, t_mix) for name, (_, expression, _, t_mix) in EMBRYO_R78_FILES.items()}
# --neighbourhood run_sorted (row E57 of experiments/harness/TRIES.md, "RS: one run per side": the score-blind ERROR
# designer's proposal P1 of the embryo design round of 2026-11-24/25 as its cold challenger fixed it; off by default):
# O's command with one last step on the final coordinates, after the cloud guard (the place --row-order takes), that
# re-assigns expression rows only (vec.run_sort.sort_rows); O's coordinates, their order, the count and the expression
# rows as a multiset are unchanged. Taken only on O's command (the test board, the default seed, nothing else changed)
# with --captain-ruling <ROW>, ROW a ruling row's ID, while TRIES.md holds exactly one row of that ID whose first cell
# is the captain's ruling ("<ROW> (the captain's ...", "ruling, no run") naming `--neighbourhood run_sorted` and
# `--captain-ruling <ROW>` (_captain_ruling_row, the builder's existing check). The placebo run_sorted_same is
# harness-only and never built.
NEIGHBOURHOODS = run_sort.BUILDABLE  # as_built (the default: every file as built), run_sorted
EMBRYO_RS_OPTION = "`--neighbourhood run_sorted`"
EMBRYO_RS_FILE = "T2embryo_s2x_ovrs.h5ad"
EMBRYO_RS_COMMAND = f"{EMBRYO_O_COMMAND} --neighbourhood run_sorted --captain-ruling {{row}}"
EMBRYO_RS_READ = ("row E57 of experiments/harness/TRIES.md (experiments/coords/e57_rs_read.py, launcher "
                  "experiments/coords/run_e57.sh)")
HEART_HALVING = ("experiments/builds/halve.py <file> <half>: a uniform random half of the rows (halve.DEFAULT_SEED, "
                 "order kept), as for the ruled heart2 file; the halving is a separate step, not this build")


# Row R71's files (the captain's order of 2026-11-25, decision log 2026-11-25 (qqq)): the embryo's consensus-form files
# BF on O and on SO only (vec.shape_interp.BAYES_FORMS), each taken only on its base's command with --shape-interp
# overlay_bayes --captain-ruling R71, and built (experiments/builds/embryo_bf.sh) only after its gate row reads PASS.
BUILDABLE_SHAPES = (*shape_interp.FORMS, *shape_interp.COMPOSED_FORMS, *shape_interp.BAYES_FORMS)


class _OnceAction(argparse.Action):
    """The embryo command's --captain-ruling, taken once: argparse would keep the last of two, so a second refuses."""

    def __call__(self, parser, namespace, values, option_string=None):
        if getattr(namespace, self.dest, None) is not None:
            parser.error(f"{option_string} is given once")
        setattr(namespace, self.dest, values)


EMBRYO_BF_ROW = "R71"
EMBRYO_BF_GATE_ROW = "E59"
EMBRYO_BF_WORDS = "BF on O and on SO only"  # words of row R71's first cell that the builder's check requires
EMBRYO_BF_COMMANDS = {  # file -> its command (O's or SO's command with the form and the ruling)
    "T2embryo_s2x_ovbf.h5ad": ("embryo --method fgw_shape0 --positions floor_exact --shape-interp overlay_bayes "
                               "--shape-interp-test-ok --captain-ruling R71"),
    "T2embryo_s2xs_ovbf.h5ad": ("embryo --method fgw_shape0 --positions floor_exact --expression series --shape-interp "
                                "overlay_bayes --shape-interp-test-ok --captain-ruling R71"),
}
EMBRYO_BF_BASES = {"T2embryo_s2x_ovbf.h5ad": "T2embryo_s2x_ov.h5ad (O)",
                   "T2embryo_s2xs_ovbf.h5ad": "T2embryo_s2xs_ov.h5ad (SO)"}
EMBRYO_BF_ORDER = ("experiments/harness/TRIES.md row R71, the captain's order of 2026-11-25 (\"My named exception on the "
                   "embryo's overlay form and the z-extent frozen rule for each BF file, as R59, with the challenger's "
                   "fixes.\"; decision log 2026-11-25 (qqq)): the consensus-form files BF on O and on SO only, an extra "
                   "file beside each base, taken only with --captain-ruling R71")
EMBRYO_BF_GATE = (f"row {EMBRYO_BF_GATE_ROW} of experiments/harness/TRIES.md (experiments/coords/bf_read.py, launcher "
                  "experiments/coords/run_bf.sh); a file is built only after the row's outcome cell reads PASS for it "
                  "(experiments/builds/embryo_bf.sh)")
EMBRYO_BF_RULE = ("vec.shape_interp.overlay_bayes on the base's overlay union at t = 1/2: q is the equal-weight mean of "
                  "the canonical densities of the released specimens " + ", ".join(shape_interp.BAYES_PRIOR_STAGES)
                  + " (duplicates dropped, eight flip copies, 48-cube over [-3, 3], Gaussian SD 0.1875); 4,000 rows of "
                  "the union and 4,000 draws from q at every count (streams [seed, 62] and [seed, 61]), the exact "
                  "assignment, each row the displacement of its nearest drawn row, the Gaussian 30-nearest-cell "
                  "average; the pool moved alike; rows, order, count, expression, size rule, scope, seed and every "
                  "other stream as the base's")


def _view(key: str, task: str, kind: str, setting: str | None, test_board: str, roles: dict[str, tuple[str, float]],
          reference_role: str, target_day: float) -> harness.TaskView:
    """A TaskView for the test geometry (the harness's view of a pseudo-task, with the test's stages)."""
    return harness.TaskView(key=key, task=task, kind=kind, setting=setting, test_board=test_board,
                            roles=tuple(roles), role_days={r: d for r, (_, d) in roles.items()},
                            role_stems={r: harness.stem(f) for r, (f, _) in roles.items()},
                            reference_role=reference_role, target_day=target_day)


def _board(which: str, override: str | None):
    boards = load_boards()
    key = override or TEST_BOARDS[which]
    if key not in boards:
        sys.exit(f"{key} is not in boards/: run the board intake first, or pass --board <validation board> for a dry "
                 "run")
    dry = key != TEST_BOARDS[which]
    return boards[key], dry


def _t2_count(board, pool_n: int) -> dict:
    """Approved rule 7 through its committed decision file (runbook section 5), checked against the count the harness's
    builder draws (harness.count_for on the candidate board); returns the decision's provenance, no counts."""
    n_rule = harness.count_for(harness._candidate_board(board), pool_n)
    n_dec, rec = t2_build.n_cells_from_decision(board, pool_n, COUNT_DECISION)
    if n_dec != n_rule:
        sys.exit("the committed count decision gives another count than the builder's count rule (values not shown); "
                 "build through n_cells_from_decision before any upload")
    d = rec["count_decision"]
    return {k: d.get(k) for k in ("decision_file", "decision_sha256", "kind", "validation_target", "override",
                                  "decision")}


def _hide_counts(rec: dict, board) -> dict:
    """On a test board, the record without cell counts: n_cells becomes "max_cells" (the count rule), and anything
    else refuses. T1's realised composition is left out (below)."""
    copy_n = rec.get(COPY_ROWS_N)  # --copy-rows only: the floor copy's own count, checked and then dropped
    if copy_n is not None:
        rec = {k: v for k, v in rec.items() if k != COPY_ROWS_N}
    if not board.counts_hidden:
        return rec
    if copy_n is not None:  # off by default: the copy-rows file holds the floor copy's count, not max_cells
        if rec["n_cells"] != copy_n:
            sys.exit("the built file does not hold the floor copy's count (values not shown)")
        out = dict(rec, n_cells="floor_copy")
        if isinstance(out.get("record"), dict):
            out["record"] = {k: v for k, v in out["record"].items() if k not in COUNT_FIELDS}
        return out
    if rec["n_cells"] != board.max_cells:
        sys.exit("the built file does not hold the board's max_cells, the approved count (values not shown)")
    out = dict(rec, n_cells="max_cells")
    # Stripped, not rounded (docs/PLAN.md decision log 2026-11-06 (g)): T1's realised composition is quota / n_cells per
    # category (t1_scope.composition of the drawn cells), so at full precision its denominators give the closed test
    # count. Rounding would not make it safe: the quotas are vec.core.quotas of the target shares, which the record
    # keeps at full precision ("shares"), so a rounded realised share still moves with n_cells wherever a target sits
    # within 1/n_cells of a rounding boundary. And it adds nothing: the quotas meet the targets to within 1/n_cells by
    # construction. The h5ad is untouched: this is the record written beside it.
    out.pop("realised", None)
    if isinstance(out.get("record"), dict):
        out["record"] = {k: v for k, v in out["record"].items() if k not in COUNT_FIELDS}
    return out


def _check_seed(seed: int) -> int:
    if seed < harness.BUILD_SEED_OFFSET:
        sys.exit(f"build seed {seed} is below {harness.BUILD_SEED_OFFSET} (runbook section 5: never seed 0)")
    return seed


def _t1_momentum(hold: bool, load=load_stage) -> dict:
    """build_submission's momentum arguments for T1. R1 (experiments/harness/TRIES.md; the captain's ruling of
    2026-11-05): with hold, cells without a Feng group keep their E10.5 expression, so no momentum is built and E9.5
    is not loaded; without it they get B0's damped E9.5 -> E10.5 momentum (M1 as built before R1)."""
    if hold:
        return {}
    return {"start": load(DATA / "T1/E9.5_RNA.h5ad", lineage=False), "damping": harness.B0_DAMPING,
            "start_source": "E9.5_RNA"}


def _feng_hygiene_on(flag: bool | None, hold_unbracketed: bool) -> bool:
    """R2's build flag (experiments/harness/TRIES.md R2). Since the build commit after N12 (docs/PLAN.md decision log
    2026-11-06 (k)), --feng-hygiene not given (None) follows --hold-unbracketed: on with R1, the default, and off with
    --no-hold-unbracketed, which keeps rebuilding the archived s1 of ae6d6b3 byte-identically. --no-feng-hygiene
    rebuilds the s1 of 88b057f (M1 with R1, before R2). An explicit --feng-hygiene with --no-hold-unbracketed is
    refused (R2 departs from M1 with R1)."""
    if flag and not hold_unbracketed:
        sys.exit("--feng-hygiene is refused with --no-hold-unbracketed: R2 departs from M1 with R1 "
                 "(experiments/harness/TRIES.md R2)")
    return bool(hold_unbracketed) if flag is None else bool(flag)


def _no_blood(direction_groups) -> None:
    """R2's guard, on group names only: --feng-hygiene is refused while Blood is among Feng's direction groups."""
    if "Blood" in set(direction_groups):
        sys.exit("--feng-hygiene is refused: Blood is among Feng's direction groups (experiments/harness/TRIES.md R2)")


def _feng_directions(genes, feng_groups: str = FENG_GROUPS_DEFAULT):
    """feng_dirs.feng_directions on the groups --feng-groups names. The default makes exactly the call the entries were
    built with (no keyword), so the default build is unchanged; "none" asks for no direction and reads no Feng file;
    "recomputed" asks for the permitted-only groups of experiments/external/feng_relabel.py (refused there unless that
    file is written through the labeller's gate, hashed in its manifest, usable, made from the yaml as it stands and
    its manifest committed as the labeller's RECORD), and is refused here unless the yaml's terminal rule routes that
    relabelling to recomputed."""
    if feng_groups == FENG_GROUPS_DEFAULT:
        return feng_dirs.feng_directions(genes)
    if feng_groups != "recomputed":
        return feng_dirs.feng_directions(genes, groups=feng_groups)
    import feng_relabel

    try:
        directions, ext = feng_dirs.feng_directions(genes, groups=feng_groups)
    except (ValueError, FileNotFoundError, StopIteration, feng_relabel.Refused) as e:
        sys.exit(f"--feng-groups recomputed is refused: {e}")
    if ext.get("route") != "recomputed":
        sys.exit(f"--feng-groups recomputed is refused: feng_relabel.yaml's terminal rule routes this relabelling to "
                 f"'{ext.get('route')}' (fallbacks {ext.get('fallbacks')}); the rule's file is --feng-groups "
                 f"{ext.get('route')}")
    return directions, ext


def _g_mom_covered(groups: dict, directions, feng_groups: str = FENG_GROUPS_DEFAULT) -> set[str]:
    """The labels g_mom's hook checks its edits against. Default: the labels the build moved along Feng's direction
    (as before). Under --feng-groups none: the labels the readout config's regexes give a group (`_labels_without_group`
    names the others), the same rule as the hook's own mask under none, so the check is a self-check there.
    vec.t1_state_transport.R1_FENG_GROUPS, a list fixed in N12 by counting Feng cells, is not read. Under recomputed:
    the default rule, the labels whose group has a recomputed direction after the fallback (`directions` holds only
    those groups), which is also the hook's own mask there (`_labels_without_direction`)."""
    if feng_groups == "none":
        return {lab for lab, group in groups.items() if not group.startswith(NO_GROUP)}
    return {lab for lab, group in groups.items() if group in directions}


def _labels_without_direction(labels, with_direction) -> set[str]:
    """Under --feng-groups recomputed, the labels g_mom holds back inside the categories it edits: those whose group
    (the readout config's regexes, `feng_dirs.label_groups`) has no recomputed direction after the fallback, labels
    with no group included. It takes the place of R1's mask (vec.t1_state_transport.r1_unbracketed_labels), whose list
    R1_FENG_GROUPS was fixed in N12 by counting cells under the authors' labels, so R1's mask in the build (cells of a
    group without a direction keep their E10.5 counts) and g_mom's mask follow the same groups."""
    groups = feng_dirs.label_groups(sorted({str(x) for x in labels}))
    return {lab for lab, group in groups.items() if group not in set(with_direction)}


NO_GROUP = "other:"  # what the readout config's harmonisation gives a label that none of its regexes matches


def _labels_without_group(labels) -> set[str]:
    """The labels the readout config's regexes map to no harmonised group (`feng_dirs.label_groups`: "other:..."). Under
    --feng-groups none this is the set g_mom holds back inside the categories it edits: the readout config's own
    rule, in place of R1's mask (vec.t1_state_transport.r1_unbracketed_labels), which reads R1_FENG_GROUPS, a list
    fixed in N12 by counting Feng cells. On every label configs/t1_chamber_classes.yaml can place in cm or
    endocardium the two masks agree: both hold back Early VIC and nothing else (tests/test_t1_feng_groups_none.py
    pins this on the committed configs)."""
    groups = feng_dirs.label_groups(sorted({str(x) for x in labels}))
    return {lab for lab, group in groups.items() if group.startswith(NO_GROUP)}


FENG_NONE_COMMAND = "t1 --g-mom --base lib3"  # the first upload's command, the one design --feng-groups none builds


def _feng_none_refusals(a) -> None:
    """--feng-groups none takes Feng's direction out and puts nothing in its place, and builds one file: the first
    upload's command (FENG_NONE_COMMAND). Refused before anything is read: without --g-mom, and with --base pooled or
    mid, --cm-share, --blood-share, --lib3-pair or --herd, which no row prices with it; then the options that act on or
    beside that direction: --no-hold-unbracketed (B0's momentum on every cell, a file no row prices), an explicit
    --feng-hygiene (R2 holds genes flat in a direction that is not there) and --strong-null (its prior acts on the
    cells Feng's direction moves)."""
    if (not a.g_mom or a.base != "lib3" or a.cm_share is not None or a.blood_share is not None or a.lib3_pair
            or a.herd):
        sys.exit(f"--feng-groups none builds one file: the first upload's command ({FENG_NONE_COMMAND}) with Feng's "
                 "direction taken out; it is refused with --base pooled or mid, --cm-share, --blood-share, "
                 "--lib3-pair, --herd or without --g-mom, which no row prices with it")
    if not a.hold_unbracketed:
        sys.exit("--feng-groups none is refused with --no-hold-unbracketed: with no Feng group every cell would get "
                 "B0's momentum, a file no row of experiments/harness/TRIES.md prices")
    if a.feng_hygiene:
        sys.exit("--feng-groups none is refused with --feng-hygiene: R2 holds genes flat in Feng's direction, and "
                 "this file has none")
    if a.strong_null:
        sys.exit("--feng-groups none is refused with --strong-null: its prior acts on the cells Feng's direction moves")
    _none_after_relabelling(getattr(a, "feng_none_ruling", None))


def _relabel_record() -> dict | None:
    """The permitted-only relabelling's record (experiments/external/feng_relabel.py RECORD), as in the tree or else
    as committed at HEAD; None while no relabelling exists. Its route is the yaml's terminal rule's outcome."""
    import feng_relabel

    raw = feng_relabel.RECORD.read_bytes() if feng_relabel.RECORD.exists() else feng_relabel.committed_at_head(
        feng_relabel.RECORD)
    return None if raw is None else json.loads(raw)


def _none_after_relabelling(ruling: str | None) -> None:
    """The final review's lesser point 1, for --feng-groups none. Before any relabelling exists, none builds as on main
    (and --feng-none-ruling is refused: it has nothing to override). After a relabelling whose recorded route is
    recomputed, none is refused unless --feng-none-ruling names a later docs/PLAN.md decision-log entry of the captain,
    naming N54 and containing the fixed phrase "sets the N54 route aside" (feng_relabel.NONE_RULING_PHRASE): not the
    ruling that authorised the run (the record's gate), dated on or after the record's commit, and not yet in
    docs/PLAN.md at that commit, which must be readable there (feng_relabel.later_ruling; the review of the port, lesser
    point 2; the final check of the fix round, items E and F). A switch to none made after N54's listing (which compares
    with the authors' groups) would otherwise be a choice by similarity to the authors' labels."""
    import feng_relabel

    rec = _relabel_record()
    if rec is None or rec.get("route") != "recomputed":
        if ruling is not None:
            sys.exit("--feng-none-ruling is refused: no relabelling routes to recomputed, so --feng-groups none needs "
                     "no ruling")
        return
    why = ("no --feng-none-ruling" if ruling is None else
           feng_relabel.later_ruling(ruling, feng_relabel.RECORD, feng_relabel.ROW,
                                     excluded=[(rec.get("gate") or {}).get("ruling")],
                                     phrase=feng_relabel.NONE_RULING_PHRASE))
    if why is not None:
        sys.exit("--feng-groups none is refused: the permitted-only relabelling's record routes T1's file to "
                 "recomputed (feng_relabel.yaml's terminal rule); only a later ruling of the captain, a docs/PLAN.md "
                 "decision-log entry naming N54 and containing \"sets the N54 route aside\", given as "
                 "--feng-none-ruling, other than the run's own ruling, dated on or after the record's commit and not "
                 f"yet in docs/PLAN.md then, sets that aside ({why})")


def _feng_recomputed_refusals(a) -> None:
    """--feng-groups recomputed builds one file: the first upload's command (FENG_NONE_COMMAND, t1 --g-mom --base lib3)
    on the recomputed groups, with R1 and R2 as in that upload. Refused before anything is read: without --g-mom, with
    --base pooled or mid, --cm-share, --blood-share, --lib3-pair, --herd, --no-hold-unbracketed, --no-feng-hygiene,
    --strong-null or --feng-none-ruling; then unless the labeller's record is committed at HEAD with the route
    recomputed (the yaml's terminal rule). Row R37's ruled exception: --base pooled is taken when the command is the
    one row R37 orders (_pooled_ordered: --captain-ruling R37, no candidate's, MS's or bet's option, and row R37's
    first cell the captain's named exception naming the pair), every other check as it is; without it --base pooled
    is refused here as before."""
    import feng_relabel

    pooled = a.base == "pooled" and _pooled_ordered(a)  # row R37's one file; nothing else in R13's command changes
    if (not a.g_mom or (a.base != "lib3" and not pooled) or a.cm_share is not None or a.blood_share is not None
            or a.lib3_pair or a.herd or not a.hold_unbracketed or a.feng_hygiene is False or a.strong_null
            or getattr(a, "feng_none_ruling", None) is not None):
        if pooled:
            sys.exit(f"--base pooled with --feng-groups recomputed builds row {T1_R37_ROW}'s one file, R13's command "
                     f"on the pooled base ({T1_POOLED_COMMAND}), nothing else changed; it is refused with --cm-share, "
                     "--blood-share, --lib3-pair, --herd, --no-hold-unbracketed, --no-feng-hygiene, --strong-null, "
                     "--feng-none-ruling or without --g-mom")
        sys.exit(f"--feng-groups recomputed builds one file: the first upload's command ({FENG_NONE_COMMAND}) on the "
                 "recomputed groups; it is refused with --base pooled or mid, --cm-share, --blood-share, --lib3-pair, "
                 "--herd, --no-hold-unbracketed, --no-feng-hygiene, --strong-null, --feng-none-ruling or without "
                 "--g-mom")
    head = feng_relabel.committed_at_head(feng_relabel.RECORD)
    if head is None:
        sys.exit("--feng-groups recomputed is refused: no permitted-only relabelling is committed "
                 f"({feng_relabel.RECORD.relative_to(REPO_ROOT)} at HEAD)")
    route = json.loads(head).get("route")
    if route != "recomputed":
        sys.exit(f"--feng-groups recomputed is refused: the committed relabelling's route is '{route}' "
                 "(feng_relabel.yaml's terminal rule); the rule's file is --feng-groups none")


# --released-bet (try T61 of experiments/harness/TRIES.md; off by default): a variant of T1's last file (the command
# FENG_NONE_COMMAND with --feng-groups none, row R8) that changes which cells are drawn. Both bets touch the
# post-exposure block (docs/DATA_POLICY.md 2026-11-04 (d): T1's composition, except by its legal route): each is
# refused unless --captain-ruling names a row of experiments/harness/TRIES.md that is the captain's ruling and names
# the arm. Read like for like: the trend and the axes run from E9.5 (embryo portions) to E10.5's portion libraries
# (every library but the dissected heart, T1_DISSECTED_LIBRARY), never to library 3, so no step is fitted on library
# 3's cells alone (row R6).
#   share_trend  library 3's own chamber mix x (E10.5 portions' share / E9.5's share) ** strength, renormalised
#                (vec.t1_bets.trend_shares; strength in t1_bets.LAMBDAS)
#   advance      library 3's cells drawn inside each category by quotas over five z-quantile bins, in proportion to
#                each bin's sum of exp(strength z), z the cell's projection on its category's E9.5 -> E10.5-portions
#                momentum axis (vec.t1_bets.advance_reweight; quotas, never importance weights: docs/DATA_POLICY.md
#                2026-11-04 (f)(2); strength in t1_bets.TAUS)
RELEASED_BETS = {"share_trend": ("l", t1_bets.LAMBDAS), "advance": ("t", t1_bets.TAUS)}
TRIES_FILE = Path(__file__).resolve().parents[2] / "experiments" / "harness" / "TRIES.md"


def _bet_arm(bet: str, strength: float) -> str:
    """The arm's name in row T61 (e.g. share_trend_l0.5, advance_t1)."""
    return f"{bet}_{RELEASED_BETS[bet][0]}{strength:g}"


def _captain_ruling_names(row: str, arm: str, path: Path | None = None) -> bool:
    """Whether experiments/harness/TRIES.md (TRIES_FILE) holds a ruling row `row` (its ID starts with R) whose first
    cell holds "(the captain's ruling" and whose text puts the arm in with the dedicated phrase "rules `<arm>` in". A
    row that merely mentions the arm, refuses it or awaits a ruling does not pass."""
    if not row.startswith("R"):
        return False
    for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines():
        if not ln.startswith(f"| {row} "):
            continue
        first = ln[2:].split("|", 1)[0]
        if "(the captain's ruling" in first and f"rules `{arm}` in" in ln:
            return True
    return False


def _bet_refusals(a) -> None:
    """--released-bet's refusals, before anything is read: only with --feng-groups none (whose own refusals then hold
    the rest of the command to the last file's); the strength must be one row T61 reads; and --captain-ruling must name
    the captain's ruling row for that arm. --captain-ruling or --bet-strength without --released-bet is refused too."""
    if _l12_ordered(a):  # rows R60 and R46 (_t1_l12 refused a released bet's option before this)
        return
    bet = getattr(a, "released_bet", None)
    if bet is None:
        ruling = getattr(a, "captain_ruling", None)
        if ruling is not None and (_ruled_fraction_ordered(a) or _ruled_damping_ordered(a)  # rows R20, R22, R23, R26:
                                   or _ruled_line3_ordered(a) or _ruled_line5_ordered(a) or _ms_ordered(a)
                                   or _msb_ordered(a) or _ruled_line6_ordered(a) or _ruled_line7_ordered(a)
                                   or _pooled_ordered(a)  # and R37; then the rows named on the command line for
                                   or _hm_ordered(a) or _ms25_ordered(a)  # Blood's momentum and MS at 2.5; and
                                   or _r51_ordered(a)  # row R51's seven files; and rows R54's and R58's files
                                   or _r54_ordered(a) or _r58_ordered(a)
                                   or _r67_ordered(a)):  # and row R67's three files
            return  # _t1_candidate, _t1_pooled, _t1_hm, _t1_ms and _t1_r51 refused every other use
        if ruling is not None or getattr(a, "bet_strength", None) is not None:
            sys.exit("--captain-ruling and --bet-strength go with --released-bet only")
        return
    if getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "none":
        sys.exit("--released-bet is a variant of T1's last file: it needs --feng-groups none "
                 f"({FENG_NONE_COMMAND} --feng-groups none)")
    allowed = RELEASED_BETS[bet][1]
    if a.bet_strength is None or a.bet_strength not in allowed:
        sys.exit(f"--released-bet {bet} needs --bet-strength, one of {list(allowed)} (row T61's arms)")
    arm = _bet_arm(bet, a.bet_strength)
    if not a.captain_ruling or not _captain_ruling_names(a.captain_ruling, arm):
        sys.exit(f"--released-bet {bet} is refused: it changes T1's composition, a lever closed by the post-exposure "
                 "block, and only the captain's own ruling opens it; --captain-ruling must name a row R<n> of "
                 "experiments/harness/TRIES.md whose first cell holds \"(the captain's ruling\" and whose text says "
                 f"\"rules `{arm}` in\"")


# The T1 candidates of docs/PLAN.md decision log 2026-11-20 (i), each a single change to the recomputed last file
# (FENG_NONE_COMMAND with --feng-groups recomputed). Off by default (None: the option not given, today's path and record
# exactly). The first value of each tuple is the default and today's behaviour; the second is the candidate's.
G_MOM_DAMPINGS = (harness.B0_DAMPING, 1.0)  # A: g_mom's factors finished at 1.0 (practice arm gmom_s1, rows T49, T61)
FENG_FRACTION_MULTS = (1.0, 1.5)  # B: the direction at 1.5 x the calendar fraction (rows T50, T58, T63)
T1_CANDIDATE_COMMAND = f"{FENG_NONE_COMMAND} --feng-groups recomputed"
# Row R20 of experiments/harness/TRIES.md (the captain's order of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (x),
# the amendment of docs/S2_RULES.md that names R20): one more file, T1_CANDIDATE_COMMAND with the direction at 2.0 x the
# calendar fraction and nothing else changed. Kept apart from FENG_FRACTION_MULTS, whose two values need no ruling: this
# one is taken only with --captain-ruling R20, and only while that row is the captain's order naming the option.
FENG_FRACTION_MULT_R20 = 2.0
T1_R20_ROW = "R20"
T1_R20_OPTION = "`--feng-fraction-mult 2.0`"
# Row R22 of experiments/harness/TRIES.md (the captain's order of 2026-11-22, docs/PLAN.md decision log 2026-11-22 (j),
# the amendment of docs/S2_RULES.md "a score-guided search of round steps may run"): the first search line, the next
# half steps of the same scale, 2.5 and 3.0 x the calendar fraction on T1_CANDIDATE_COMMAND with nothing else changed.
# Each is taken only with --captain-ruling R22, and only while that row is the captain's order naming the option.
FENG_FRACTION_MULTS_R22 = (2.5, 3.0)
T1_R22_ROW = "R22"
T1_R22_OPTIONS = {2.5: "`--feng-fraction-mult 2.5`", 3.0: "`--feng-fraction-mult 3.0`"}
# Every value taken only under a ruling row: value -> (the row, the option as the row's first cell names it, the
# record's order text). R20's entry is the 2.0 file exactly as before.
FENG_FRACTION_MULTS_RULED = {
    FENG_FRACTION_MULT_R20: (T1_R20_ROW, T1_R20_OPTION,
                             f"experiments/harness/TRIES.md row {T1_R20_ROW}, the captain's order of 2026-11-21 "
                             "(docs/PLAN.md decision log 2026-11-21 (x)): the direction at 2.0 x the calendar "
                             "fraction, nothing else changed"),
    **{v: (T1_R22_ROW, T1_R22_OPTIONS[v],
           f"experiments/harness/TRIES.md row {T1_R22_ROW}, the captain's order of 2026-11-22 (docs/PLAN.md decision "
           f"log 2026-11-22 (j)), the first search line: the direction at {v:.1f} x the calendar fraction, nothing "
           "else changed") for v in FENG_FRACTION_MULTS_R22},
}
# Row R22's second T1 search line (the lines added in docs/PLAN.md decision log 2026-11-22): g_mom's within-type
# damping at 1.5 and at 2.0 on T1_CANDIDATE_COMMAND, nothing else changed (fraction multiplier 1.0). Kept apart from
# G_MOM_DAMPINGS, whose two values need no ruling (1.0 is candidate A, as before): each is taken only with
# --captain-ruling R22, and only while that row is the captain's order naming the option.
G_MOM_DAMPINGS_R22 = (1.5, 2.0)
T1_R22_DAMPING_OPTIONS = {1.5: "`--g-mom-damping 1.5`", 2.0: "`--g-mom-damping 2.0`"}
G_MOM_DAMPING_ORDER_R22 = ("experiments/harness/TRIES.md row {row}, the captain's order of 2026-11-22 (docs/PLAN.md "
                           "decision log 2026-11-22), the second T1 search line: g_mom's within-type damping at "
                           "{v:.1f}, nothing else changed")
# Row R22's third T1 search line (the line added in docs/PLAN.md decision log 2026-11-22): the direction at 2.0 and at
# 2.5 x the calendar fraction on file C's command, g_mom's within-type damping held at 1.0 (candidate A's value).
# Kept apart from FENG_FRACTION_MULTS_RULED, whose 2.0 (row R20's file) and 2.5 (the first line's) are taken at damping
# 0.5 only: each pair is taken only with --captain-ruling R22, and only while row R22's first cell names the combined
# option exactly; the bare multiplier's option, which the first line's 2.5 names, never orders the pair.
G_MOM_DAMPING_LINE3 = G_MOM_DAMPINGS[1]
FENG_FRACTION_MULTS_LINE3 = (2.0, 2.5)
T1_R22_LINE3_OPTIONS = {2.0: "`--g-mom-damping 1.0 --feng-fraction-mult 2.0`",
                        2.5: "`--g-mom-damping 1.0 --feng-fraction-mult 2.5`"}
T1_LINE3_ORDER_R22 = ("experiments/harness/TRIES.md row {row}, the captain's order of 2026-11-22 (docs/PLAN.md "
                      "decision log 2026-11-22), the third T1 search line: on file C's command, g_mom's within-type "
                      "damping held at 1.0 and the direction at {v:.1f} x the calendar fraction, nothing else changed")
# Row R22's fifth T1 search line (T1 line 5 of docs/PLAN.md's decision log): g_mom's within-type damping at 1.5 and at
# 2.0 on the 2.0 file's command (the direction at 2.0 x the calendar fraction), nothing else changed (files
# T1_s2rcB2d15 and T1_s2rcB2d20). Kept apart from FENG_FRACTION_MULTS_RULED, whose 2.0 (row R20's file) is taken at
# damping 0.5 only, and from G_MOM_DAMPINGS_R22, whose 1.5 and 2.0 (the second line's) are taken at multiplier 1.0
# only: each pair is taken only with --captain-ruling R22, and only while row R22's first cell names the combined
# option exactly; the bare options of the second line and of row R20 never order the pair.
FENG_FRACTION_MULT_LINE5 = FENG_FRACTION_MULT_R20
G_MOM_DAMPINGS_LINE5 = G_MOM_DAMPINGS_R22
T1_R22_LINE5_OPTIONS = {1.5: "`--g-mom-damping 1.5 --feng-fraction-mult 2.0`",
                        2.0: "`--g-mom-damping 2.0 --feng-fraction-mult 2.0`"}
T1_LINE5_ORDER_R22 = ("experiments/harness/TRIES.md row {row}, the captain's order (docs/PLAN.md decision log, T1 line "
                      "5), the fifth T1 search line: on the 2.0 file's command (the direction at 2.0 x the calendar "
                      "fraction), g_mom's within-type damping at {v:.1f}, nothing else changed")
# Row R22's sixth T1 search line (T1 line 6 of docs/PLAN.md decision log 2026-11-24 (bb); the captain's order of
# 2026-11-24, ruling row R40): g_mom's within-type damping at 1.5 and at 2.0 on the 2.5 file's command (the direction at
# 2.5 x the calendar fraction, the first line's T1_s2rcB25), nothing else changed (files T1_s2rcB25d15 and
# T1_s2rcB25d20). Kept apart from FENG_FRACTION_MULTS_RULED, whose 2.5 is taken at damping 0.5 only, from
# G_MOM_DAMPINGS_R22 (the second line's, at multiplier 1.0 only) and from the fifth line's pairs (multiplier 2.0): each
# pair is taken only with --captain-ruling R22, and only while row R22's first cell names the combined option exactly;
# the bare options never order the pair.
FENG_FRACTION_MULT_LINE6 = FENG_FRACTION_MULTS_R22[0]
G_MOM_DAMPINGS_LINE6 = G_MOM_DAMPINGS_R22
T1_R22_LINE6_OPTIONS = {1.5: "`--g-mom-damping 1.5 --feng-fraction-mult 2.5`",
                        2.0: "`--g-mom-damping 2.0 --feng-fraction-mult 2.5`"}
T1_LINE6_FILES = {1.5: "T1_s2rcB25d15", 2.0: "T1_s2rcB25d20"}
T1_LINE6_ORDER_R22 = ("experiments/harness/TRIES.md row {row}, the captain's order of 2026-11-24 (ruling row R40, "
                      "\"yes for two extra Task 1 steps\"; docs/PLAN.md decision log 2026-11-24 (bb), T1 line 6), the "
                      "sixth T1 search line: on the 2.5 file's command (the direction at 2.5 x the calendar fraction), "
                      "g_mom's within-type damping at {v:.1f}, nothing else changed")
# Row R22's seventh T1 search line (T1 line 7 of docs/PLAN.md decision log 2026-11-24 (bb), the third line's one further
# half step; the captain's order of 2026-11-24, ruling row R39): file C's command with the direction at 3.0 x the
# calendar fraction, g_mom's within-type damping held at 1.0, nothing else changed (file T1_s2rcCm30). Kept apart from
# FENG_FRACTION_MULTS_RULED, whose 3.0 (the first line's T1_s2rcB30) is taken at damping 0.5 only, and from the third
# line's pairs (2.0 and 2.5): the pair is taken only with --captain-ruling R22, and only while row R22's first cell
# names the combined option exactly; the bare `--feng-fraction-mult 3.0` never orders it.
G_MOM_DAMPING_LINE7 = G_MOM_DAMPING_LINE3
FENG_FRACTION_MULT_LINE7 = FENG_FRACTION_MULTS_R22[1]
T1_R22_LINE7_OPTION = "`--g-mom-damping 1.0 --feng-fraction-mult 3.0`"
T1_LINE7_FILE = "T1_s2rcCm30"
T1_LINE7_ORDER_R22 = ("experiments/harness/TRIES.md row {row}, the captain's order of 2026-11-24 (ruling row R39, "
                      "\"yes for two extra Task 1 steps\", lifting the third line's cap by one half step; docs/PLAN.md "
                      "decision log 2026-11-24 (bb), T1 line 7), the seventh T1 search line: on file C's command, "
                      "g_mom's within-type damping held at 1.0 and the direction at {v:.1f} x the calendar fraction, "
                      "nothing else changed")
# Row R37 of experiments/harness/TRIES.md (the captain's named exception of 2026-11-24, "Yes, T1 pooled base, my
# exception"; docs/PLAN.md decision log 2026-11-24 (z) and (bb)): T1's depth, closed after the exposure of 4 November,
# changes for one file, R13's command (T1_CANDIDATE_COMMAND) on the pooled base: every E10.5 chamber cell at the
# stage's own mix in place of library 3's cells at library 3's mix (file T1_s2rc_pooled). --feng-groups recomputed takes
# --base pooled only with --captain-ruling R37, and only while row R37's first cell is the captain's named exception
# naming the pair as T1_R37_OPTION writes it; nothing else in R13's command changes.
T1_R37_ROW = "R37"
T1_R37_KIND = "named exception"  # the row's first cell starts "R37 (the captain's named exception"
T1_R37_OPTION = "`--base pooled` with `--feng-groups recomputed`"
T1_POOLED_FILE = "T1_s2rc_pooled.h5ad"
T1_POOLED_COMMAND = f"t1 --g-mom --base pooled --feng-groups recomputed --captain-ruling {T1_R37_ROW}"
T1_POOLED_ORDER = ("experiments/harness/TRIES.md row R37, the captain's named exception of 2026-11-24 (\"Yes, T1 "
                   "pooled base, my exception\"; docs/PLAN.md decision log 2026-11-24 (z) and (bb)): T1's depth, "
                   "closed after the exposure of 4 November, changes for one file, R13's command on the pooled base "
                   "(every E10.5 chamber cell at the stage's own mix in place of library 3's cells at library 3's "
                   "mix); nothing else changed")
# Row R23 of experiments/harness/TRIES.md (the captain's order of 2026-11-23, docs/PLAN.md decision log 2026-11-23 (c)):
# candidate MS, T1_CANDIDATE_COMMAND with one change. The drawn cells whose label's group has a recomputed direction
# (build_submission's covered cells) also get B0's per-category damped mean momentum (vec.t1_scope.momentum_shift at
# harness.B0_DAMPING with its own pseudocount eps 0.5, min_cells 20 and keep_depth, on build_submission's momentum
# stream [seed, 1], held genes the build's held list: R2's union joined with the pre-registered list), fitted from E9.5
# to E10.5's portion libraries 1 and 2 (library 3, the dissected heart, never enters the fit), applied before the atlas
# step ([seed, 2], the calendar fraction) and g_mom ([seed, 5]), which stay as they are. Cells without a group keep R1's
# hold. As the practice arm m1rc_gmom_mom of experiments/t1_m1/ms_read.py (row T65) builds it, with the fit's stages
# the test's. Off by default; taken only with --captain-ruling R23 while that row is the captain's order naming the
# option, never with candidate A, B or C's options or a released bet.
DIRECTION_MOMENTUM_ARMS = ("all",)  # the one arm row R23 names
T1_R23_ROW = "R23"
T1_R23_OPTION = "`--direction-momentum all`"
T1_MS_COMMAND = f"{T1_CANDIDATE_COMMAND} --direction-momentum all --captain-ruling {T1_R23_ROW}"
T1_MS_STREAM = 1  # build_submission's momentum stream, [seed, 1]
T1_MS_DIRECTION_STREAM = 2  # build_submission's direction stream, [seed, 2]
T1_MS_FIT_LIBRARIES = ("1", "2")  # E10.5's portion libraries, the fit's later end; library 3 never enters it
T1_MS_ORDER = ("experiments/harness/TRIES.md row R23, the captain's order of 2026-11-23 (docs/PLAN.md decision log "
               "2026-11-23 (c)): candidate MS, R13's file with B0's released per-category damped mean momentum added "
               "on the cells whose label's group has a recomputed direction, before the atlas step and g_mom (arm "
               "'all'); nothing else changed")
# Row R26 of experiments/harness/TRIES.md (the captain's order of 2026-11-23, "combine MS"; docs/PLAN.md decision log
# 2026-11-23 (j), the T1 combined MS list): one file, T1_s2rcMSB, candidate MS's step exactly as row R23 fixes it
# (_t1_ms_step, unchanged) put on candidate B's command exactly as B has it (the atlas step at 1.5 x the calendar
# fraction, g_mom's within-type damping 0.5 passed as B's build passes it), nothing else changed. Kept apart from row
# R23's constants, which stay as they are: the pair is taken only with --captain-ruling R26, and only while row R26's
# first cell is the captain's order naming the pair exactly. No constant is new: the multiplier is candidate B's.
T1_R26_ROW = "R26"
T1_MSB_MULT = FENG_FRACTION_MULTS[1]  # candidate B's multiplier, 1.5
T1_R26_OPTION = "`--feng-fraction-mult 1.5 --direction-momentum all`"
T1_MSB_FILE = "T1_s2rcMSB.h5ad"
T1_MSB_COMMAND = (f"{T1_CANDIDATE_COMMAND} --feng-fraction-mult 1.5 --direction-momentum all "
                  f"--captain-ruling {T1_R26_ROW}")
T1_MSB_ORDER = ("experiments/harness/TRIES.md row R26, the captain's order of 2026-11-23 (\"combine MS\"; docs/PLAN.md "
                "decision log 2026-11-23 (j), the T1 combined MS list): its one file, T1_s2rcMSB, candidate MS's step "
                "exactly as row R23 fixes it put on candidate B's command exactly as B has it (the atlas step at 1.5 x "
                "the calendar fraction, g_mom's within-type damping 0.5); nothing else changed")
# The two steps of the T1 design round of 2026-11-24 put on the banked file, T1_s2rcB25 (T1_CANDIDATE_COMMAND with the
# atlas step at 2.5 x the calendar fraction, row R22's first line, g_mom's damping 0.5), each a single change, each
# taken only with --captain-ruling <ROW>: ROW is given on the command line (a ruling row's ID: R, digits, an optional
# letter), and experiments/harness/TRIES.md must hold exactly one row of that ID whose first cell is the captain's
# ruling ("<ROW> (the captain's ...", "ruling, no run"; any kind: order, yes, named exception, words) naming the file's
# option string and `--captain-ruling <ROW>`, each in backticks (_captain_ruling_row).
T1_RULING_ROW_ID = r"R[0-9]+[a-z]?"  # ASCII digits, as the launcher's own check
T1_BANKED_FILE = "T1_s2rcB25"
T1_BANKED_MULT = FENG_FRACTION_MULTS_R22[0]  # 2.5
# Files F1, F7 and F8 as the round's cold challenge fixes them (its FIX (a)): the drawn cells labelled Blood, which R1
# holds at their drawn counts, get Blood's own released change, B0's per-category damped mean momentum
# (vec.t1_scope.momentum_shift with _t1_ms_step's constants: eps 0.5, min_cells 20, keep_depth, the build's held genes)
# fitted on Blood-labelled cells at both ends (within={"blood": ["Blood"]}: E10.5's first Macrophages never enter
# Blood's factors) from E9.5 to E10.5's portion libraries 1 and 2 (T1_MS_FIT_LIBRARIES; library 3 never in the fit), at
# damping d, on its own stream [seed, 6] (apart from the draw's, MS's [seed, 1], the atlas step's [seed, 2] and g_mom's
# [seed, 5]), after the atlas step and before g_mom, which leaves Blood alone. Macrophage, Early VIC and every other
# cell keep the banked file's path, so their rows are the banked file's byte for byte. It reverses R1 for Blood in one
# file.
T1_HM_CHOICES = ("0.5", "1.0", "2.0")  # the option's literal values, as a ruling row names them
T1_HM_OPTIONS = {0.5: "`--held-momentum 0.5`", 1.0: "`--held-momentum 1.0`", 2.0: "`--held-momentum 2.0`"}
T1_HM_FILES = {0.5: "T1_s2rcB25H05", 1.0: "T1_s2rcB25H10", 2.0: "T1_s2rcB25H20"}
T1_HM_LABEL = "Blood"  # the one label the step moves
T1_HM_CATEGORY = "blood"  # its chamber category (configs/t1_chamber_classes.yaml)
T1_HM_STREAM = 6
T1_HM_EPS = 0.5  # vec.t1_scope.momentum_shift's defaults, the constants _t1_ms_step's calls take
T1_HM_MIN_CELLS = 20
T1_HM_COMMAND = T1_CANDIDATE_COMMAND + " --feng-fraction-mult 2.5 --held-momentum {value} --captain-ruling {row}"
T1_HM_ORDER = ("experiments/harness/TRIES.md row {row}, the captain's ruling naming {option} (the T1 design round of "
               "2026-11-24, files F1, F7 and F8 as its cold challenge fixes them): the banked file T1_s2rcB25's "
               "command (the atlas step at 2.5 x the calendar fraction, g_mom's within-type damping 0.5) with Blood's "
               "own released E9.5 -> E10.5 change at damping {value} on the drawn Blood cells, after the atlas step "
               "and before g_mom; R1 reversed for Blood in this one file, Macrophage and Early VIC held; nothing else "
               "changed")
# File F2 as the round's cold challenge fixes it: candidate MS's step exactly as row R23 fixes it (_t1_ms_step,
# unchanged) on the banked file's command, nothing else changed; a third MS file beside rows R23's and R26's.
T1_MS25_OPTION = "`--feng-fraction-mult 2.5 --direction-momentum all`"
T1_MS25_FILE = "T1_s2rcB25MS"
T1_MS25_COMMAND = T1_CANDIDATE_COMMAND + " --feng-fraction-mult 2.5 --direction-momentum all --captain-ruling {row}"
T1_MS25_ORDER = ("experiments/harness/TRIES.md row {row}, the captain's ruling naming " + T1_MS25_OPTION + " (the T1 "
                 "design round of 2026-11-24, file F2 as its cold challenge fixes it): its file, T1_s2rcB25MS, "
                 "candidate MS's step exactly as row R23 fixes it put on the banked file T1_s2rcB25's command (the "
                 "atlas step at 2.5 x the calendar fraction, g_mom's within-type damping 0.5); nothing else changed")
# The options a held-momentum file never takes beside its own and the banked file's: flag -> attribute
HM_STACKED = (("--g-mom-damping", "g_mom_damping"), ("--released-bet", "released_bet"),
              ("--bet-strength", "bet_strength"))
# Row R51 of experiments/harness/TRIES.md (the captain's order of 2026-11-24, "ruling, no run"; decision log 2026-11-24
# (iii); T1's line 10, the second design pass of the round of 2026-11-24/25 as its cold challenge fixes it, P1 struck):
# seven files, each one change to T1_s2rcMSB's command (row R26's: candidate MS's step exactly as row R23 fixes it on
# candidate B's command, the atlas step at 1.5 x the calendar fraction, g_mom's within-type damping 0.5, MS's damping
# harness.B0_DAMPING), each taken only with --captain-ruling R51 while TRIES.md holds exactly one row R51 whose first
# cell is the captain's ruling ("R51 (the captain's ...", "ruling, no run") naming the file's option string and
# `--captain-ruling R51`, each as a whole backticked string (_r51_row_names: a backticked string that holds the option
# inside a longer one, or one the option holds, does not name it). The option strings, after T1_CANDIDATE_COMMAND:
#   T1_s2rcMS20, T1_s2rcMS30      MS's step at the atlas multiplier 2.0 or 3.0 in place of 1.5
#   T1_s2rcMSBd10, T1_s2rcMSBd15  g_mom's within-type damping 1.0 or 1.5 in place of 0.5 (the same fit finished at it)
#   T1_s2rcMSBH05, T1_s2rcMSBH10  Blood's own released change at damping 0.5 or 1.0: row R44's step (_t1_hm_step: the
#                                 drawn Blood cells only, fitted on Blood-labelled cells at both ends from E9.5 to
#                                 E10.5's libraries 1 and 2, library 3 never in the fit, its own stream [seed, 6]),
#                                 after MS's step and before g_mom; R1 reversed for Blood in these two files
#   T1_s2rcMSBm10                 MS's own damping 1.0 (--direction-momentum-damping, a new option) passed to
#                                 _t1_ms_step in place of harness.B0_DAMPING
# Nothing else changes in any of them. No other combination is taken: every other command is routed as before.
T1_R51_ROW = "R51"
T1_R51_BASE_FILE = "T1_s2rcMSB"
# file -> (its option string as row R51 names it, then the values its command gives for --feng-fraction-mult,
# --g-mom-damping, --held-momentum and --direction-momentum-damping, None where the option is not given)
T1_R51_FILES = {
    "T1_s2rcMS20": ("--feng-fraction-mult 2.0 --direction-momentum all", 2.0, None, None, None),
    "T1_s2rcMS30": ("--feng-fraction-mult 3.0 --direction-momentum all", 3.0, None, None, None),
    "T1_s2rcMSBd10": ("--feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.0", 1.5, 1.0, None, None),
    "T1_s2rcMSBd15": ("--feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.5", 1.5, 1.5, None, None),
    "T1_s2rcMSBH05": ("--feng-fraction-mult 1.5 --direction-momentum all --held-momentum 0.5", 1.5, None, "0.5", None),
    "T1_s2rcMSBH10": ("--feng-fraction-mult 1.5 --direction-momentum all --held-momentum 1.0", 1.5, None, "1.0", None),
    "T1_s2rcMSBm10": ("--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.0", 1.5, None,
                      None, "1.0"),
}
T1_R51_M10_FILE = "T1_s2rcMSBm10"  # the one file that takes --direction-momentum-damping
T1_R51_DMD_CHOICES = ("1.0",)  # --direction-momentum-damping's one literal value, as row R51 names it
T1_R51_COMMAND = T1_CANDIDATE_COMMAND + " {option} --captain-ruling " + T1_R51_ROW
T1_R51_CHANGES = {  # file -> its one change, in the record's order text
    "T1_s2rcMS20": "the atlas step at 2.0 x the calendar fraction",
    "T1_s2rcMS30": "the atlas step at 3.0 x the calendar fraction",
    "T1_s2rcMSBd10": "g_mom's within-type damping 1.0",
    "T1_s2rcMSBd15": "g_mom's within-type damping 1.5",
    "T1_s2rcMSBH05": ("Blood's own released E9.5 -> E10.5 change at damping 0.5 on the drawn Blood cells, after MS's "
                      "step and before g_mom (R1 reversed for Blood in this one file, Macrophage and Early VIC held)"),
    "T1_s2rcMSBH10": ("Blood's own released E9.5 -> E10.5 change at damping 1.0 on the drawn Blood cells, after MS's "
                      "step and before g_mom (R1 reversed for Blood in this one file, Macrophage and Early VIC held)"),
    "T1_s2rcMSBm10": "MS's own damping 1.0 in place of harness.B0_DAMPING",
}
T1_R51_ORDER = ("experiments/harness/TRIES.md row R51, the captain's order of 2026-11-24 (\"ruling, no run\"; decision "
                "log 2026-11-24 (iii), T1's line 10) naming `{option}` and `--captain-ruling R51`: file {file}, "
                "T1_s2rcMSB's command (candidate MS's step exactly as row R23 fixes it on candidate B's command: the "
                "atlas step at 1.5 x the calendar fraction, g_mom's within-type damping 0.5, MS's damping "
                "harness.B0_DAMPING) with one change, {change}; nothing else changed")
# The options a row R51 file never takes beside its own option set: flag -> attribute
R51_STACKED = (("--released-bet", "released_bet"), ("--bet-strength", "bet_strength"))
# Row R58 of experiments/harness/TRIES.md (the captain's order of 2026-11-24, "ruling, no run": "Yes, T1_s2rcMSBm20
# under R22, bound to its gate G0."; T1's line 11, the step-change round's T1 PORTFOLIO design as its cold challenge
# fixes it): one file, T1_s2rcMSBm20, T1_s2rcMSB's command with MS's own damping 2.0 in place of harness.B0_DAMPING
# (--direction-momentum-damping 2.0: the released day's per-category change at full pace on both days), nothing else
# changed. Taken only with --captain-ruling R58 while TRIES.md holds exactly one row R58 whose first cell is the
# captain's order ("ruling, no run") naming T1_R58_OPTION and `--captain-ruling R58` (_captain_order_names), on
# T1_CANDIDATE_COMMAND at the default seed on the test board. Its route is tested before row R51's; row R51's value 1.0
# (T1_s2rcMSBm10) and its refusals stay as they are, so --direction-momentum-damping 2.0 under any other ruling is
# refused by row R51's guard as before. Its build is gated on row T73's read (experiments/builds/t1_lm_m20_line.sh).
T1_R58_ROW = "R58"
T1_R58_FILE = "T1_s2rcMSBm20"
T1_R58_DMD = "2.0"  # --direction-momentum-damping's literal value under row R58
T1_R58_OPTION_STRING = "--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 2.0"
T1_R58_OPTION = f"`{T1_R58_OPTION_STRING}`"
T1_R58_COMMAND = f"{T1_CANDIDATE_COMMAND} {T1_R58_OPTION_STRING} --captain-ruling {T1_R58_ROW}"
T1_R58_ORDER = ("experiments/harness/TRIES.md row R58, the captain's order of 2026-11-24 (\"Yes, T1_s2rcMSBm20 under "
                "R22, bound to its gate G0.\"; T1's line 11) naming " + T1_R58_OPTION + " and `--captain-ruling R58`: "
                "file T1_s2rcMSBm20, T1_s2rcMSB's command with one change, MS's own damping 2.0 in place of "
                "harness.B0_DAMPING; nothing else changed")
T1_R58_NOTE = ("--direction-momentum-damping 2.0 (row R58's file T1_s2rcMSBm20): MS's own damping in place of "
               "harness.B0_DAMPING (0.5); nothing else changed")
# The options row R58's file never takes beside its own option set: flag -> attribute
R58_STACKED = (("--g-mom-damping", "g_mom_damping"), ("--held-momentum", "held_momentum"),
               ("--released-bet", "released_bet"), ("--bet-strength", "bet_strength"), ("--late-mix", "late_mix"))
# Row R67 of experiments/harness/TRIES.md (the captain's order of 2026-11-25, "ruling, no run": "Yes to T1_s2rcMSBm125,
# T1_s2rcMSBm15 and T1_s2rcMSBm075 under R22, sent at once, in one ruling row: MSB's command with MS's own damping at
# 1.25, 1.5 and 0.75."; decision log 2026-11-25 (kkk)): three files, each T1_s2rcMSB's command with one change, MS's
# own damping in place of harness.B0_DAMPING, nothing else changed. Row R58's route over three literals: T1_R67_TABLE
# is the one table, the row's three (file, option string) pairs in the row's order; each file's value is the literal
# string its option string ends in (_r67_value, which also holds the name to the value: T1_s2rcMSBm + the value's
# digits without the point, the first digit the units), and its record's note (T1_R67_NOTES) carries that literal,
# never a float formatted again. A file is taken only with --captain-ruling R67 while TRIES.md holds exactly one row
# R67 whose first cell is the captain's order ("ruling, no run") naming the file's option string and
# `--captain-ruling R67`, each as a whole backticked string, the option string followed there by " (<file>)" as the
# row pairs each string with its file (_r67_row_names), on T1_CANDIDATE_COMMAND at the default seed on the test board.
# Its route is tested after row R58's and before row R51's; every other command is routed as before (the three values
# under any other ruling or none are refused by row R51's guard as before). Built by experiments/builds/t1_r67_line.sh.
T1_R67_ROW = "R67"
T1_R67_TABLE = (
    ("T1_s2rcMSBm125", "--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.25"),
    ("T1_s2rcMSBm15", "--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.5"),
    ("T1_s2rcMSBm075", "--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 0.75"),
)
# MSB's option set (row R26's) and the flag every row R67 option string ends in, before its value
T1_R67_PREFIX = "--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping "


def _r67_value(name: str, option: str) -> str:
    """A row R67 option string's literal value, as the builder's parser gives it to --direction-momentum-damping: the
    string must be T1_R67_PREFIX followed by one decimal (digits, a point, digits) that is not B0's damping, row R51's
    1.0 or row R58's 2.0, and the file's name must be T1_s2rcMSBm followed by that decimal's digits without the point;
    else ValueError."""
    value = option.removeprefix(T1_R67_PREFIX)
    if value == option or not re.fullmatch(r"\d+\.\d+", value):
        raise ValueError(f"not a row {T1_R67_ROW} option string: {option}")
    if float(value) in (harness.B0_DAMPING, 1.0, 2.0) or name != "T1_s2rcMSBm" + value.replace(".", ""):
        raise ValueError(f"row {T1_R67_ROW}'s file {name} does not match its option string {option}")
    return value


# file -> (its option string as row R67 names it, the literal value --direction-momentum-damping takes for it)
T1_R67_FILES = {name: (option, _r67_value(name, option)) for name, option in T1_R67_TABLE}
T1_R67_COMMAND = T1_CANDIDATE_COMMAND + " {option} --captain-ruling " + T1_R67_ROW
T1_R67_ORDER = ("experiments/harness/TRIES.md row R67, the captain's order of 2026-11-25 (\"Yes to T1_s2rcMSBm125, "
                "T1_s2rcMSBm15 and T1_s2rcMSBm075 under R22, sent at once, in one ruling row: MSB's command with MS's "
                "own damping at 1.25, 1.5 and 0.75.\"; decision log 2026-11-25 (kkk)) naming `{option}` and "
                "`--captain-ruling R67`: file {file}, T1_s2rcMSB's command with one change, MS's own damping {value} "
                "in place of harness.B0_DAMPING; nothing else changed")
# file -> MS's record note for it (the damping_option field), as T1_R58_NOTE is row R58's file's
T1_R67_NOTES = {name: (f"--direction-momentum-damping {value} (row {T1_R67_ROW}'s file {name}): MS's own damping in "
                       f"place of harness.B0_DAMPING ({float(harness.B0_DAMPING):g}); nothing else changed")
                for name, (_, value) in T1_R67_FILES.items()}
# The options a row R67 file never takes beside its own option set: row R58's, flag -> attribute
R67_STACKED = R58_STACKED
# Row R54 of experiments/harness/TRIES.md (the captain's named exception of 2026-11-24, written as his order, "ruling,
# no run": "Yes, T1_s2rcMSB_LM16, my named exception, bound to T72: ..."; T1's late-mix line, the step-change round's
# T1 CLOSED design as its cold challenge fixes it, the E14.5 file struck): one file, T1_s2rcMSB_LM16, T1_s2rcMSB's
# command with library 3's own chamber mix moved to the shares row N64's step commits
# (experiments/t1_scope/late_mix.py, LATE_MIX_SHARES): cardiomyocytes of heart cells, then blood of chamber cells, by
# vec.t1_scope.with_share in --cm-share's and --blood-share's order, the other categories in proportion; quotas at the
# count rule, never weights; every step after the draw is T1_s2rcMSB's. Taken only with --captain-ruling R54 while
# TRIES.md holds exactly one row R54 whose first cell is the captain's order ("ruling, no run") naming T1_R54_OPTION and
# `--captain-ruling R54` (_captain_order_names), on T1_s2rcMSB's command at the default seed on the test board, and
# only while the shares file committed at HEAD (the tree's copy the same bytes) has the sha256 row N64 names as
# `shares.json sha256 <hex>` (exactly one such string in that row) and records the pool's own mix it was computed
# from. Its route is tested before every other, row R51's included; for this one file it takes the place of row R6's
# refusal of moved shares on --base lib3 (--cm-share and --blood-share stay refused there). Its build is gated on row
# T72's read (experiments/builds/t1_lm_m20_line.sh).
T1_R54_ROW = "R54"
LATE_MIX_CHOICES = ("e15_17",)  # the anchor row R54 rules: Feng's permitted CD1 cells at E15.5 to E17.5 (row N57)
T1_R54_OPTION = "`--late-mix e15_17`"
T1_R54_FILE = "T1_s2rcMSB_LM16"
T1_R54_COMMAND = (f"{T1_CANDIDATE_COMMAND} --feng-fraction-mult 1.5 --direction-momentum all --late-mix e15_17 "
                  f"--captain-ruling {T1_R54_ROW}")
LATE_MIX_ROW = "N64"  # the step-0 row whose cells name the shares file's sha256
LATE_MIX_SHARES_PATH = "experiments/t1_scope/results/n64_late_mix/shares.json"  # row N64's output, as committed
LATE_MIX_SHARES = REPO_ROOT / LATE_MIX_SHARES_PATH
LATE_MIX_SHA_NAMED = re.compile(r"`shares\.json sha256 ([0-9a-f]{64})`")  # how row N64 names it
LATE_MIX_TOL = 1e-12  # the pool's own mix against the recorded one; the moved shares against the recorded ones
T1_R54_ORDER = ("experiments/harness/TRIES.md row R54, the captain's named exception of 2026-11-24 (\"ruling, no "
                "run\"; T1's late-mix line) naming `--late-mix e15_17` and `--captain-ruling R54`: file "
                "T1_s2rcMSB_LM16, T1_s2rcMSB's command with library 3's own chamber mix moved to row N64's shares "
                "(cardiomyocytes of heart cells, then blood of chamber cells, the other categories in proportion); "
                "nothing else changed")
# The options row R54's file never takes beside T1_s2rcMSB's own: flag -> attribute (and the three switches below)
R54_STACKED = (("--g-mom-damping", "g_mom_damping"), ("--held-momentum", "held_momentum"),
               ("--direction-momentum-damping", "direction_momentum_damping"), ("--released-bet", "released_bet"),
               ("--bet-strength", "bet_strength"), ("--cm-share", "cm_share"), ("--blood-share", "blood_share"))
R54_SWITCHES = (("--lib3-pair", "lib3_pair"), ("--herd", "herd"), ("--strong-null", "strong_null"))


def _captain_order_names(row: str, option: str, path: Path | None = None, kind: str = "order") -> bool:
    """Whether experiments/harness/TRIES.md (TRIES_FILE) holds exactly one row `row`, a ruling row: its first cell
    starts with "<row> (the captain's order" and says "ruling, no run", and that first cell names `option` (in
    backticks) and "`--captain-ruling <row>`". A row that only mentions the option elsewhere does not pass. `kind`
    (default "order", every call before row R37) is the ruling's kind as the first cell opens with it: row R37's is
    "named exception" (T1_R37_KIND), so its first cell must start with "R37 (the captain's named exception"."""
    lines = [ln for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {row} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    return (first.startswith(f"{row} (the captain's {kind}") and "ruling, no run" in first and option in first
            and f"`--captain-ruling {row}`" in first)


def _ruled_fraction_ordered(a) -> bool:
    """The commands the ruling rows order, as _t1_candidate holds them: --feng-fraction-mult 2.0 with --captain-ruling
    R20, or 2.5 or 3.0 with --captain-ruling R22, no --released-bet or --bet-strength, and the value's row the
    captain's order naming its option."""
    ruled = FENG_FRACTION_MULTS_RULED.get(getattr(a, "feng_fraction_mult", None))
    return (ruled is not None and getattr(a, "captain_ruling", None) == ruled[0]
            and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None
            and _captain_order_names(ruled[0], ruled[1]))


def _ruled_damping_ordered(a) -> bool:
    """The commands row R22's second line orders, as _t1_candidate holds them: --g-mom-damping 1.5 or 2.0 with
    --captain-ruling R22, --feng-fraction-mult absent or 1.0, no --released-bet or --bet-strength, and row R22 the
    captain's order naming the option."""
    damping = getattr(a, "g_mom_damping", None)
    return (damping in T1_R22_DAMPING_OPTIONS and getattr(a, "captain_ruling", None) == T1_R22_ROW
            and getattr(a, "feng_fraction_mult", None) in (None, FENG_FRACTION_MULTS[0])
            and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None
            and _captain_order_names(T1_R22_ROW, T1_R22_DAMPING_OPTIONS[damping]))


def _is_line3(a) -> bool:
    """Whether the command gives the third line's pair: --g-mom-damping 1.0 with --feng-fraction-mult 2.0 or 2.5."""
    return (getattr(a, "g_mom_damping", None) == G_MOM_DAMPING_LINE3
            and getattr(a, "feng_fraction_mult", None) in T1_R22_LINE3_OPTIONS)


def _ruled_line3_ordered(a) -> bool:
    """The commands row R22's third line orders, as _t1_candidate holds them: --g-mom-damping 1.0 with
    --feng-fraction-mult 2.0 or 2.5 and --captain-ruling R22, no --released-bet or --bet-strength, and row R22 the
    captain's order naming the combined option."""
    return (_is_line3(a) and getattr(a, "captain_ruling", None) == T1_R22_ROW
            and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None
            and _captain_order_names(T1_R22_ROW, T1_R22_LINE3_OPTIONS[a.feng_fraction_mult]))


def _t1_ruled_line3(a, mult: float) -> dict:
    """--g-mom-damping 1.0 with --feng-fraction-mult 2.0 or 2.5 (row R22's third T1 search line): refused before
    anything is read without --captain-ruling R22, with any other row, unless row R22 is the captain's order naming the
    combined option, with --released-bet or --bet-strength, and on any command but T1_CANDIDATE_COMMAND at the default
    seed on the test board. Returns the record's candidate values, naming the row."""
    row, option = T1_R22_ROW, T1_R22_LINE3_OPTIONS[mult]
    ruling = getattr(a, "captain_ruling", None)
    pair = f"--g-mom-damping {G_MOM_DAMPING_LINE3:g} with --feng-fraction-mult {mult:g}"
    if ruling is None:
        sys.exit(f"{pair} is the captain's order of row {row} alone and needs --captain-ruling {row}")
    if ruling != row:
        sys.exit(f"{pair} is refused with --captain-ruling {ruling}: only row {row} orders it")
    if not _captain_order_names(row, option):
        sys.exit(f"{pair} is refused: experiments/harness/TRIES.md has no single row {row} whose first cell is the "
                 f"captain's order (\"ruling, no run\") naming {option} and `--captain-ruling {row}`")
    if getattr(a, "released_bet", None) is not None or getattr(a, "bet_strength", None) is not None:
        sys.exit(f"{pair} is refused with --released-bet and --bet-strength: row {row}'s line is one file with one "
                 "change")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"{pair} changes T1's recomputed last file only: it needs the command {T1_CANDIDATE_COMMAND} at the "
                 f"default seed on the test board (row {row})")
    return {"g_mom_damping": float(G_MOM_DAMPING_LINE3), "feng_fraction_mult": float(mult), "captain_ruling": row,
            "order": T1_LINE3_ORDER_R22.format(row=row, v=float(mult))}


def _is_line5(a) -> bool:
    """Whether the command gives the fifth line's pair: --g-mom-damping 1.5 or 2.0 with --feng-fraction-mult 2.0."""
    return (getattr(a, "g_mom_damping", None) in T1_R22_LINE5_OPTIONS
            and getattr(a, "feng_fraction_mult", None) == FENG_FRACTION_MULT_LINE5)


def _ruled_line5_ordered(a) -> bool:
    """The commands row R22's fifth line orders, as _t1_candidate holds them: --g-mom-damping 1.5 or 2.0 with
    --feng-fraction-mult 2.0 and --captain-ruling R22, no --released-bet or --bet-strength, and row R22 the captain's
    order naming the combined option."""
    return (_is_line5(a) and getattr(a, "captain_ruling", None) == T1_R22_ROW
            and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None
            and _captain_order_names(T1_R22_ROW, T1_R22_LINE5_OPTIONS[a.g_mom_damping]))


def _t1_ruled_line5(a, damping: float) -> dict:
    """--g-mom-damping 1.5 or 2.0 with --feng-fraction-mult 2.0 (row R22's fifth T1 search line): refused before
    anything is read without --captain-ruling R22, with any other row, unless row R22 is the captain's order naming the
    combined option, with --released-bet or --bet-strength, and on any command but T1_CANDIDATE_COMMAND at the default
    seed on the test board. Returns the record's candidate values, naming the row."""
    row, option = T1_R22_ROW, T1_R22_LINE5_OPTIONS[damping]
    ruling = getattr(a, "captain_ruling", None)
    pair = f"--g-mom-damping {damping:g} with --feng-fraction-mult {FENG_FRACTION_MULT_LINE5:g}"
    if ruling is None:
        sys.exit(f"{pair} is the captain's order of row {row} alone and needs --captain-ruling {row}")
    if ruling != row:
        sys.exit(f"{pair} is refused with --captain-ruling {ruling}: only row {row} orders it")
    if not _captain_order_names(row, option):
        sys.exit(f"{pair} is refused: experiments/harness/TRIES.md has no single row {row} whose first cell is the "
                 f"captain's order (\"ruling, no run\") naming {option} and `--captain-ruling {row}`")
    if getattr(a, "released_bet", None) is not None or getattr(a, "bet_strength", None) is not None:
        sys.exit(f"{pair} is refused with --released-bet and --bet-strength: row {row}'s line is one file with one "
                 "change")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"{pair} changes T1's recomputed last file only: it needs the command {T1_CANDIDATE_COMMAND} at the "
                 f"default seed on the test board (row {row})")
    return {"g_mom_damping": float(damping), "feng_fraction_mult": float(FENG_FRACTION_MULT_LINE5),
            "captain_ruling": row, "order": T1_LINE5_ORDER_R22.format(row=row, v=float(damping))}


def _is_line6(a) -> bool:
    """Whether the command gives the sixth line's pair: --g-mom-damping 1.5 or 2.0 with --feng-fraction-mult 2.5."""
    return (getattr(a, "g_mom_damping", None) in T1_R22_LINE6_OPTIONS
            and getattr(a, "feng_fraction_mult", None) == FENG_FRACTION_MULT_LINE6)


def _ruled_line6_ordered(a) -> bool:
    """The commands row R22's sixth line orders, as _t1_candidate holds them: --g-mom-damping 1.5 or 2.0 with
    --feng-fraction-mult 2.5 and --captain-ruling R22, no --released-bet or --bet-strength, and row R22 the captain's
    order naming the combined option."""
    return (_is_line6(a) and getattr(a, "captain_ruling", None) == T1_R22_ROW
            and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None
            and _captain_order_names(T1_R22_ROW, T1_R22_LINE6_OPTIONS[a.g_mom_damping]))


def _is_line7(a) -> bool:
    """Whether the command gives the seventh line's pair: --g-mom-damping 1.0 with --feng-fraction-mult 3.0."""
    return (getattr(a, "g_mom_damping", None) == G_MOM_DAMPING_LINE7
            and getattr(a, "feng_fraction_mult", None) == FENG_FRACTION_MULT_LINE7)


def _ruled_line7_ordered(a) -> bool:
    """The command row R22's seventh line orders, as _t1_candidate holds it: --g-mom-damping 1.0 with
    --feng-fraction-mult 3.0 and --captain-ruling R22, no --released-bet or --bet-strength, and row R22 the captain's
    order naming the combined option."""
    return (_is_line7(a) and getattr(a, "captain_ruling", None) == T1_R22_ROW
            and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None
            and _captain_order_names(T1_R22_ROW, T1_R22_LINE7_OPTION))


def _r22_pair_refusals(a, pair: str, option: str) -> None:
    """The refusals of a combined option of row R22 (the sixth and seventh T1 search lines, as the third and fifth
    lines' own): before anything is read, without --captain-ruling R22, with any other row, unless row R22 is the
    captain's order naming `option`, with --released-bet or --bet-strength, and on any command but
    T1_CANDIDATE_COMMAND at the default seed on the test board."""
    row = T1_R22_ROW
    ruling = getattr(a, "captain_ruling", None)
    if ruling is None:
        sys.exit(f"{pair} is the captain's order of row {row} alone and needs --captain-ruling {row}")
    if ruling != row:
        sys.exit(f"{pair} is refused with --captain-ruling {ruling}: only row {row} orders it")
    if not _captain_order_names(row, option):
        sys.exit(f"{pair} is refused: experiments/harness/TRIES.md has no single row {row} whose first cell is the "
                 f"captain's order (\"ruling, no run\") naming {option} and `--captain-ruling {row}`")
    if getattr(a, "released_bet", None) is not None or getattr(a, "bet_strength", None) is not None:
        sys.exit(f"{pair} is refused with --released-bet and --bet-strength: row {row}'s line is one file with one "
                 "change")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"{pair} changes T1's recomputed last file only: it needs the command {T1_CANDIDATE_COMMAND} at the "
                 f"default seed on the test board (row {row})")


def _t1_ruled_line6(a, damping: float) -> dict:
    """--g-mom-damping 1.5 or 2.0 with --feng-fraction-mult 2.5 (row R22's sixth T1 search line, T1 line 6, files
    T1_s2rcB25d15 and T1_s2rcB25d20): refused before anything is read as _r22_pair_refusals says. Returns the record's
    candidate values, naming the row and line 6."""
    _r22_pair_refusals(a, f"--g-mom-damping {damping:g} with --feng-fraction-mult {FENG_FRACTION_MULT_LINE6:g}",
                       T1_R22_LINE6_OPTIONS[damping])
    return {"g_mom_damping": float(damping), "feng_fraction_mult": float(FENG_FRACTION_MULT_LINE6),
            "captain_ruling": T1_R22_ROW, "order": T1_LINE6_ORDER_R22.format(row=T1_R22_ROW, v=float(damping))}


def _t1_ruled_line7(a) -> dict:
    """--g-mom-damping 1.0 with --feng-fraction-mult 3.0 (row R22's seventh T1 search line, T1 line 7, file
    T1_s2rcCm30): refused before anything is read as _r22_pair_refusals says. Returns the record's candidate values,
    naming the row, line 7 and the captain's order."""
    _r22_pair_refusals(a, f"--g-mom-damping {G_MOM_DAMPING_LINE7:g} with --feng-fraction-mult "
                          f"{FENG_FRACTION_MULT_LINE7:g}", T1_R22_LINE7_OPTION)
    return {"g_mom_damping": float(G_MOM_DAMPING_LINE7), "feng_fraction_mult": float(FENG_FRACTION_MULT_LINE7),
            "captain_ruling": T1_R22_ROW,
            "order": T1_LINE7_ORDER_R22.format(row=T1_R22_ROW, v=float(FENG_FRACTION_MULT_LINE7))}


def _t1_ruled_damping(a, damping: float) -> dict:
    """--g-mom-damping 1.5 or 2.0 (row R22's second T1 search line): refused before anything is read without
    --captain-ruling R22, with any other row, unless row R22 is the captain's order naming the option, with
    --released-bet or --bet-strength, with --feng-fraction-mult other than 1.0, and on any command but
    T1_CANDIDATE_COMMAND at the default seed on the test board. Returns the record's candidate values, naming the
    row."""
    row, option = T1_R22_ROW, T1_R22_DAMPING_OPTIONS[damping]
    mult, ruling = getattr(a, "feng_fraction_mult", None), getattr(a, "captain_ruling", None)
    if ruling is None:
        sys.exit(f"--g-mom-damping takes {G_MOM_DAMPINGS[0]:g} (the default) or {G_MOM_DAMPINGS[1]:g} only; "
                 f"{damping:g} is the captain's order of row {row} alone and needs --captain-ruling {row}")
    if ruling != row:
        sys.exit(f"--g-mom-damping {damping:g} is refused with --captain-ruling {ruling}: only row {row} orders it")
    if not _captain_order_names(row, option):
        sys.exit(f"--g-mom-damping {damping:g} is refused: experiments/harness/TRIES.md has no single row {row} whose "
                 f"first cell is the captain's order (\"ruling, no run\") naming {option} and "
                 f"`--captain-ruling {row}`")
    if getattr(a, "released_bet", None) is not None or getattr(a, "bet_strength", None) is not None:
        sys.exit(f"--g-mom-damping {damping:g} is refused with --released-bet and --bet-strength: row {row}'s line "
                 "is one file with one change")
    if mult is not None and mult != FENG_FRACTION_MULTS[0]:
        sys.exit(f"--g-mom-damping {damping:g} is refused with --feng-fraction-mult other than "
                 f"{FENG_FRACTION_MULTS[0]:g}: row {row}'s line is one file with one change")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"--g-mom-damping {damping:g} changes T1's recomputed last file only: it needs the command "
                 f"{T1_CANDIDATE_COMMAND} at the default seed on the test board (row {row})")
    return {"g_mom_damping": float(damping), "feng_fraction_mult": FENG_FRACTION_MULTS[0], "captain_ruling": row,
            "order": G_MOM_DAMPING_ORDER_R22.format(row=row, v=float(damping))}


def _t1_ruled(a, mult: float) -> dict:
    """--feng-fraction-mult 2.0 (row R20's file) or 2.5 or 3.0 (row R22's first search line): refused before anything
    is read without --captain-ruling naming the value's row, with any other row, unless that row is the captain's order
    naming the option, with --released-bet or --bet-strength, with --g-mom-damping other than the default, and on any
    command but T1_CANDIDATE_COMMAND at the default seed on the test board. Returns the record's candidate values,
    naming the row. At 2.0 the messages and the record are row R20's as before."""
    row, option, order = FENG_FRACTION_MULTS_RULED[mult]
    damping, ruling = getattr(a, "g_mom_damping", None), getattr(a, "captain_ruling", None)
    if ruling is None:
        sys.exit(f"--feng-fraction-mult takes {FENG_FRACTION_MULTS[0]:g} (the default) or {FENG_FRACTION_MULTS[1]:g} "
                 f"only; {mult:g} is the captain's order of row {row} alone and needs --captain-ruling {row}")
    if ruling != row:
        sys.exit(f"--feng-fraction-mult {mult:g} is refused with --captain-ruling {ruling}: only row {row} orders it")
    if not _captain_order_names(row, option):
        sys.exit(f"--feng-fraction-mult {mult:g} is refused: experiments/harness/TRIES.md has no single row {row} "
                 f"whose first cell is the captain's order (\"ruling, no run\") naming {option} and "
                 f"`--captain-ruling {row}`")
    if getattr(a, "released_bet", None) is not None or getattr(a, "bet_strength", None) is not None:
        sys.exit(f"--feng-fraction-mult {mult:g} is refused with --released-bet and --bet-strength: row {row} is one "
                 "file with one change")
    if damping is not None and damping != G_MOM_DAMPINGS[0]:
        sys.exit(f"--feng-fraction-mult {mult:g} is refused with --g-mom-damping other than {G_MOM_DAMPINGS[0]:g}: "
                 f"row {row} is one file with one change")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"--feng-fraction-mult {mult:g} changes T1's recomputed last file only: it needs the command "
                 f"{T1_CANDIDATE_COMMAND} at the default seed on the test board (row {row})")
    return {"g_mom_damping": G_MOM_DAMPINGS[0], "feng_fraction_mult": float(mult), "captain_ruling": row,
            "order": order}


def _t1_candidate(a) -> dict | None:
    """--g-mom-damping and --feng-fraction-mult: refused before anything is read with a value outside their two, on any
    command but T1_CANDIDATE_COMMAND (--feng-groups recomputed, whose own refusals then hold every other option to the
    first upload's), with a seed other than the default and on a dry run. Returns None when neither is given (today's
    file and record), else the two values as the record names them. --feng-fraction-mult 2.0 is row R20's file alone,
    2.5 and 3.0 row R22's line alone (_t1_ruled); --g-mom-damping 1.5 and 2.0 row R22's second T1 line alone
    (_t1_ruled_damping); --g-mom-damping 1.0 with --feng-fraction-mult 2.0 or 2.5 row R22's third T1 line alone
    (_t1_ruled_line3, tested before the multiplier's own routing, which keeps every other pair as before);
    --g-mom-damping 1.5 or 2.0 with --feng-fraction-mult 2.0 row R22's fifth T1 line alone (_t1_ruled_line5, tested
    beside the third line's, before the multiplier's routing and the damping's); --g-mom-damping 1.5 or 2.0 with
    --feng-fraction-mult 2.5 row R22's sixth T1 line alone (_t1_ruled_line6, tested beside the fifth line's, before
    the multiplier's routing, which keeps refusing 2.5 with any damping but 0.5); --g-mom-damping 1.0 with
    --feng-fraction-mult 3.0 row R22's seventh T1 line alone (_t1_ruled_line7, tested before the multiplier's routing,
    which keeps refusing 3.0 with any damping but 0.5); with the flag --captain-ruling, any other value is refused
    here or by _bet_refusals. With --held-momentum (_t1_hm's refusals first) and with MS's pair at 2.5 under a ruling
    row naming it (_t1_ms25's), the banked file's two values, naming that row: both are tested before every other
    route, and no other command reaches them."""
    if getattr(a, "held_momentum", None) is not None:  # Blood's momentum on the banked file
        hm = _t1_hm(a)
        return {k: hm[k] for k in ("g_mom_damping", "feng_fraction_mult", "captain_ruling", "order")}
    damping, mult = getattr(a, "g_mom_damping", None), getattr(a, "feng_fraction_mult", None)
    if damping is None and mult is None:
        return None
    if _is_ms25(a):  # MS on the banked file under its ruling row
        ms = _t1_ms25(a)
        return {k: ms[k] for k in ("g_mom_damping", "feng_fraction_mult", "captain_ruling", "order")}
    if _is_line3(a):
        return _t1_ruled_line3(a, mult)
    if _is_line5(a):
        return _t1_ruled_line5(a, damping)
    if _is_line6(a):
        return _t1_ruled_line6(a, damping)
    if _is_line7(a):
        return _t1_ruled_line7(a)
    if mult in FENG_FRACTION_MULTS_RULED:
        return _t1_ruled(a, mult)
    if damping in T1_R22_DAMPING_OPTIONS:
        return _t1_ruled_damping(a, damping)
    if damping is not None and damping not in G_MOM_DAMPINGS:
        sys.exit(f"--g-mom-damping takes {G_MOM_DAMPINGS[0]:g} (the default) or {G_MOM_DAMPINGS[1]:g} only "
                 f"({G_MOM_DAMPINGS_R22[0]:g} and {G_MOM_DAMPINGS_R22[1]:g} only with --captain-ruling {T1_R22_ROW})")
    if mult is not None and mult not in FENG_FRACTION_MULTS:
        sys.exit(f"--feng-fraction-mult takes {FENG_FRACTION_MULTS[0]:g} (the default) or {FENG_FRACTION_MULTS[1]:g} "
                 f"only ({FENG_FRACTION_MULT_R20:g} only with --captain-ruling {T1_R20_ROW}; "
                 f"{FENG_FRACTION_MULTS_R22[0]:g} and {FENG_FRACTION_MULTS_R22[1]:g} only with --captain-ruling "
                 f"{T1_R22_ROW})")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit("--g-mom-damping and --feng-fraction-mult change T1's recomputed last file only: they need the "
                 f"command {T1_CANDIDATE_COMMAND} at the default seed on the test board (docs/PLAN.md decision log "
                 "2026-11-20 (i))")
    return {"g_mom_damping": G_MOM_DAMPINGS[0] if damping is None else float(damping),
            "feng_fraction_mult": FENG_FRACTION_MULTS[0] if mult is None else float(mult),
            "order": "docs/PLAN.md decision log 2026-11-20 (i): A = damping 1.0, B = 1.5 x the calendar fraction, C = "
                     "both; no ruling"}


MS_STACKED = (("--g-mom-damping", "g_mom_damping"), ("--feng-fraction-mult", "feng_fraction_mult"),
              ("--released-bet", "released_bet"), ("--bet-strength", "bet_strength"))


def _ms_ordered(a) -> bool:
    """The command row R23 orders, as _t1_ms holds it: --direction-momentum all with --captain-ruling R23, none of
    candidate A, B or C's options or a released bet's, and row R23 the captain's order naming the option."""
    return (getattr(a, "direction_momentum", None) in DIRECTION_MOMENTUM_ARMS
            and getattr(a, "captain_ruling", None) == T1_R23_ROW
            and all(getattr(a, attr, None) is None for _, attr in MS_STACKED)
            and _captain_order_names(T1_R23_ROW, T1_R23_OPTION))


def _is_msb(a) -> bool:
    """Whether the command gives row R26's pair with its ruling: --direction-momentum all, --feng-fraction-mult 1.5 and
    --captain-ruling R26. Any other use of the pair (no ruling, R23 or another row) keeps row R23's path and
    refusals."""
    return (getattr(a, "direction_momentum", None) in DIRECTION_MOMENTUM_ARMS
            and getattr(a, "feng_fraction_mult", None) == T1_MSB_MULT
            and getattr(a, "captain_ruling", None) == T1_R26_ROW)


def _msb_ordered(a) -> bool:
    """The command row R26 orders, as _t1_msb holds it: the pair with --captain-ruling R26, no --g-mom-damping, no
    released bet's option, and row R26 the captain's order naming the pair."""
    return (_is_msb(a) and getattr(a, "g_mom_damping", None) is None
            and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None
            and _captain_order_names(T1_R26_ROW, T1_R26_OPTION))


def _t1_msb(a) -> dict:
    """--feng-fraction-mult 1.5 --direction-momentum all --captain-ruling R26 (row R26, the T1 combined MS list's one
    file): refused before anything is read unless row R26 is the captain's order naming the pair, with --g-mom-damping
    (any value: B's command has none), --released-bet, --bet-strength or --no-hold-unbracketed, and on any command but
    T1_CANDIDATE_COMMAND at the default seed on the test board (never a dry run). Returns the record's candidate values,
    naming the row, the list's file and B's two values."""
    row, option = T1_R26_ROW, T1_R26_OPTION
    what = f"--feng-fraction-mult {T1_MSB_MULT:g} with --direction-momentum all"
    if not _captain_order_names(row, option):
        sys.exit(f"{what} is refused: experiments/harness/TRIES.md has no single row {row} whose first cell is the "
                 f"captain's order (\"ruling, no run\") naming {option} and `--captain-ruling {row}`")
    stacked = [flag for flag, attr in MS_STACKED if attr != "feng_fraction_mult" and getattr(a, attr, None) is not None]
    if stacked:
        sys.exit(f"{what} is refused with {', '.join(stacked)}: row {row}'s file is candidate MS's step on candidate "
                 "B's command, nothing else changed")
    if not getattr(a, "hold_unbracketed", True):
        sys.exit(f"{what} is refused with --no-hold-unbracketed: the cells without a group keep R1's hold (row {row})")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"{what} changes T1's recomputed last file only: it needs the command {T1_CANDIDATE_COMMAND} at the "
                 f"default seed on the test board (row {row})")
    return {"direction_momentum": a.direction_momentum, "captain_ruling": row, "g_mom_damping": G_MOM_DAMPINGS[0],
            "feng_fraction_mult": float(T1_MSB_MULT), "list_file": T1_MSB_FILE, "order": T1_MSB_ORDER}


def _t1_ms(a) -> dict | None:
    """--direction-momentum all (candidate MS, row R23): refused before anything is read without --captain-ruling R23,
    with any other row, unless row R23 is the captain's order naming the option (`_captain_order_names`, as the R20 and
    R22 options), with --g-mom-damping, --feng-fraction-mult, --released-bet or --bet-strength (MS is never stacked on
    A, B or C), with --no-hold-unbracketed, and on any command but T1_CANDIDATE_COMMAND at the default seed on the test
    board (never a dry run). Returns None when the option is not given (today's path and record exactly), else the
    record's candidate values, naming the row. With --feng-fraction-mult 1.5 and --captain-ruling R26 it is row R26's
    combined file instead (_t1_msb); with --feng-fraction-mult 2.5 and --captain-ruling <ROW> whose first cell is the
    captain's ruling naming that pair, the banked command's MS file (_t1_ms25); every other use, R23's included, is
    taken or refused as before."""
    arm = getattr(a, "direction_momentum", None)
    if arm is None:
        return None
    row, option = T1_R23_ROW, T1_R23_OPTION
    what = f"--direction-momentum {arm}"
    if arm not in DIRECTION_MOMENTUM_ARMS:
        sys.exit(f"{what} is refused: row {row} names the arm 'all' only")
    if _is_ms25(a):  # the pair at 2.5 under a ruling row naming it (file T1_s2rcB25MS); every other use as before
        return _t1_ms25(a)
    if _is_msb(a):  # row R26's pair with its ruling (the combined MS list); every other use goes on as before
        return _t1_msb(a)
    ruling = getattr(a, "captain_ruling", None)
    if ruling is None:
        sys.exit(f"{what} is the captain's order of row {row} alone and needs --captain-ruling {row}")
    if ruling != row:
        sys.exit(f"{what} is refused with --captain-ruling {ruling}: only row {row} orders it")
    if not _captain_order_names(row, option):
        sys.exit(f"{what} is refused: experiments/harness/TRIES.md has no single row {row} whose first cell is the "
                 f"captain's order (\"ruling, no run\") naming {option} and `--captain-ruling {row}`")
    stacked = [flag for flag, attr in MS_STACKED if getattr(a, attr, None) is not None]
    if stacked:
        sys.exit(f"{what} is refused with {', '.join(stacked)}: candidate MS is one change to R13's file and is never "
                 f"stacked on candidate A, B or C or a released bet (row {row})")
    if not getattr(a, "hold_unbracketed", True):
        sys.exit(f"{what} is refused with --no-hold-unbracketed: the cells without a group keep R1's hold (row {row})")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"{what} changes T1's recomputed last file only: it needs the command {T1_CANDIDATE_COMMAND} at the "
                 f"default seed on the test board (row {row})")
    return {"direction_momentum": arm, "captain_ruling": row, "order": T1_MS_ORDER}


# The options row R37's file never takes (R13's command on the pooled base, nothing else changed): flag -> attribute
POOLED_STACKED = (("--g-mom-damping", "g_mom_damping"), ("--feng-fraction-mult", "feng_fraction_mult"),
                  ("--direction-momentum", "direction_momentum"), ("--released-bet", "released_bet"),
                  ("--bet-strength", "bet_strength"))


def _is_pooled_rc(a) -> bool:
    """Whether the command gives row R37's pair with its ruling: --base pooled, --feng-groups recomputed and
    --captain-ruling R37. Without the ruling --feng-groups recomputed refuses --base pooled as before."""
    return (getattr(a, "base", None) == "pooled" and getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) == "recomputed"
            and getattr(a, "captain_ruling", None) == T1_R37_ROW)


def _pooled_ordered(a) -> bool:
    """The command row R37 orders, as _t1_pooled holds it: the pair with --captain-ruling R37, none of POOLED_STACKED's
    options, and row R37's first cell the captain's named exception naming the pair (`_captain_order_names` with kind
    T1_R37_KIND)."""
    return (_is_pooled_rc(a) and all(getattr(a, attr, None) is None for _, attr in POOLED_STACKED)
            and _captain_order_names(T1_R37_ROW, T1_R37_OPTION, kind=T1_R37_KIND))


def _t1_pooled(a) -> dict | None:
    """--base pooled with --feng-groups recomputed and --captain-ruling R37 (row R37, the captain's named exception on
    T1's depth): refused before anything is read unless row R37's first cell is the captain's named exception ("ruling,
    no run") naming `--base pooled` with `--feng-groups recomputed` and `--captain-ruling R37`, with --g-mom-damping,
    --feng-fraction-mult, --direction-momentum, --released-bet or --bet-strength (any value), with another seed and on
    a dry run; R13's command's own guards (_feng_recomputed_refusals, which then takes the pooled base) hold every
    other option. Returns None unless the pair is given with --captain-ruling R37 (every other command as before:
    without the ruling --feng-groups recomputed refuses --base pooled), else the record's values, naming the row."""
    if not _is_pooled_rc(a):
        return None
    row, option, what = T1_R37_ROW, T1_R37_OPTION, "--base pooled with --feng-groups recomputed"
    if not _captain_order_names(row, option, kind=T1_R37_KIND):
        sys.exit(f"{what} is refused: experiments/harness/TRIES.md has no single row {row} whose first cell is the "
                 f"captain's {T1_R37_KIND} (\"ruling, no run\") naming {option} and `--captain-ruling {row}`")
    stacked = [flag for flag, attr in POOLED_STACKED if getattr(a, attr, None) is not None]
    if stacked:
        sys.exit(f"{what} is refused with {', '.join(stacked)}: row {row}'s file is R13's command on the pooled base, "
                 "nothing else changed")
    if a.seed != harness.BUILD_SEED_OFFSET + 1 or a.board is not None:
        sys.exit(f"{what} is row {row}'s one file: it needs the command {T1_POOLED_COMMAND} at the default seed on the "
                 "test board")
    return {"base": "pooled", "captain_ruling": row, "order": T1_POOLED_ORDER, "command": T1_POOLED_COMMAND}


def _captain_ruling_row(row, option: str, path: Path | None = None) -> bool:
    """Whether `row`, a row ID given on the command line, is a ruling row's ID (T1_RULING_ROW_ID: R, digits, an optional
    letter) and experiments/harness/TRIES.md holds exactly one row of that ID whose first cell is a ruling of the
    captain naming `option` (in backticks) and `--captain-ruling <row>`: `_captain_order_names` with any kind of the
    captain's ruling (its first cell starts "<row> (the captain's " and says "ruling, no run"). A row that names the
    option outside its first cell, names another row's flag or is not the captain's ruling does not pass."""
    if not isinstance(row, str) or not re.fullmatch(T1_RULING_ROW_ID, row):
        return False
    return _captain_order_names(row, option, path, kind="")


def _hm_value(a) -> float | None:
    """--held-momentum's value as a number when it is one of the literal strings T1_HM_CHOICES, else None."""
    raw = getattr(a, "held_momentum", None)
    return float(raw) if raw in T1_HM_CHOICES else None


def _hm_ordered(a) -> bool:
    """The commands a ruling row orders for Blood's momentum, as _t1_hm holds them: --held-momentum 0.5, 1.0 or 2.0
    with --captain-ruling <ROW>, no released bet's option, and row ROW's first cell the captain's ruling naming the
    value's option string."""
    v = _hm_value(a)
    return (v is not None and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None
            and _captain_ruling_row(getattr(a, "captain_ruling", None), T1_HM_OPTIONS[v]))


def _t1_hm(a) -> dict | None:
    """--held-momentum 0.5, 1.0 or 2.0 (Blood's released momentum on the banked file; files T1_s2rcB25H05, H10 and
    H20): refused before anything but the ledger is read with any other value (the literal strings T1_HM_CHOICES
    only), with --direction-momentum (Blood's step with MS, file T1_s2rcB25MSH05, is named but not coded), without
    --captain-ruling, unless experiments/harness/TRIES.md holds exactly one row of the given ID whose first cell is the
    captain's ruling naming the value's option string and `--captain-ruling <ROW>` (_captain_ruling_row), with
    --g-mom-damping, --released-bet or --bet-strength, without --feng-fraction-mult 2.5 or with another multiplier,
    with --no-hold-unbracketed, and on any command but the banked file's (T1_CANDIDATE_COMMAND, whose own refusals
    then hold the rest, with --feng-fraction-mult 2.5) at the default seed on the test board (never a dry run).
    Returns None when the option is not given (every other command exactly as before), else the record's values,
    naming the value, the row, the file, the order text and the command."""
    raw = getattr(a, "held_momentum", None)
    if raw is None:
        return None
    what = f"--held-momentum {raw}"
    if raw not in T1_HM_CHOICES:
        sys.exit(f"{what} is refused: it takes {', '.join(T1_HM_CHOICES)} only, as a ruling row names them")
    if getattr(a, "direction_momentum", None) is not None:
        sys.exit(f"{what} is refused with --direction-momentum: Blood's step with MS (file T1_s2rcB25MSH05) is named "
                 "but not coded")
    v = float(raw)
    option, ruling = T1_HM_OPTIONS[v], getattr(a, "captain_ruling", None)
    if ruling is None:
        sys.exit(f"{what} reverses R1 for Blood in one file and is the captain's ruling alone: it needs "
                 "--captain-ruling <ROW>, the row of experiments/harness/TRIES.md whose first cell is his ruling "
                 f"naming {option} and `--captain-ruling <ROW>`")
    if not _captain_ruling_row(ruling, option):
        sys.exit(f"{what} is refused: experiments/harness/TRIES.md has no single row {ruling} whose first cell is the "
                 f"captain's ruling (\"ruling, no run\") naming {option} and `--captain-ruling {ruling}`")
    stacked = [flag for flag, attr in HM_STACKED if getattr(a, attr, None) is not None]
    if stacked:
        sys.exit(f"{what} is refused with {', '.join(stacked)}: its file is the banked file {T1_BANKED_FILE}'s command "
                 "with one change")
    if getattr(a, "feng_fraction_mult", None) != T1_BANKED_MULT:
        sys.exit(f"{what} is refused without --feng-fraction-mult 2.5 or with another multiplier: its file is the "
                 f"banked file {T1_BANKED_FILE}'s command with one change")
    if not getattr(a, "hold_unbracketed", True):
        sys.exit(f"{what} is refused with --no-hold-unbracketed: every cell but Blood keeps R1's hold")
    command = T1_HM_COMMAND.format(value=raw, row=ruling)
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"{what} changes the banked file {T1_BANKED_FILE} only: it needs the command {command} at the default "
                 "seed on the test board")
    return {"held_momentum": v, "captain_ruling": ruling, "g_mom_damping": G_MOM_DAMPINGS[0],
            "feng_fraction_mult": float(T1_BANKED_MULT), "file": T1_HM_FILES[v],
            "order": T1_HM_ORDER.format(row=ruling, option=option, value=raw), "command": command}


def _is_ms25(a) -> bool:
    """Whether the command gives MS's option at the banked file's multiplier under a ruling row naming the pair:
    --direction-momentum all, --feng-fraction-mult 2.5 and --captain-ruling <ROW> whose first cell is the captain's
    ruling naming T1_MS25_OPTION (_captain_ruling_row). Any other use of the pair keeps row R23's path and refusals."""
    return (getattr(a, "direction_momentum", None) in DIRECTION_MOMENTUM_ARMS
            and getattr(a, "feng_fraction_mult", None) == T1_BANKED_MULT
            and _captain_ruling_row(getattr(a, "captain_ruling", None), T1_MS25_OPTION))


def _ms25_ordered(a) -> bool:
    """The command a ruling row orders for MS on the banked command, as _t1_ms25 holds it: the pair under its ruling
    row, no --g-mom-damping, --held-momentum or released bet's option."""
    return (_is_ms25(a) and getattr(a, "g_mom_damping", None) is None and getattr(a, "held_momentum", None) is None
            and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None)


def _t1_ms25(a) -> dict:
    """--feng-fraction-mult 2.5 --direction-momentum all --captain-ruling <ROW> (file T1_s2rcB25MS): refused before
    anything is read with --held-momentum (Blood's step with MS is named but not coded), --g-mom-damping (any value),
    --released-bet, --bet-strength or --no-hold-unbracketed, and on any command but T1_CANDIDATE_COMMAND at the
    default seed on the test board (never a dry run); _is_ms25 has checked the row. Returns the record's values,
    naming the row, the file, the banked file's two values, the order text and the command."""
    row = a.captain_ruling
    what = f"--feng-fraction-mult {T1_BANKED_MULT:g} with --direction-momentum all"
    if getattr(a, "held_momentum", None) is not None:
        sys.exit(f"{what} is refused with --held-momentum: MS with Blood's step (file T1_s2rcB25MSH05) is named but "
                 "not coded")
    stacked = [flag for flag, attr in MS_STACKED if attr != "feng_fraction_mult" and getattr(a, attr, None) is not None]
    if stacked:
        sys.exit(f"{what} is refused with {', '.join(stacked)}: row {row}'s file is candidate MS's step on the banked "
                 f"file {T1_BANKED_FILE}'s command, nothing else changed")
    if not getattr(a, "hold_unbracketed", True):
        sys.exit(f"{what} is refused with --no-hold-unbracketed: the cells without a group keep R1's hold (row {row})")
    command = T1_MS25_COMMAND.format(row=row)
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"{what} changes the banked file {T1_BANKED_FILE} only: it needs the command {command} at the default "
                 "seed on the test board")
    return {"direction_momentum": a.direction_momentum, "captain_ruling": row, "g_mom_damping": G_MOM_DAMPINGS[0],
            "feng_fraction_mult": float(T1_BANKED_MULT), "file": T1_MS25_FILE,
            "order": T1_MS25_ORDER.format(row=row), "command": command}


# Row R60 (the captain's order of 2026-11-25, T1's line 12): seven files, each one change to T1_s2rcMS20's command (row
# R51's file), each taken only with --captain-ruling R60 while row R60's first cell is the captain's order naming the
# file's option string (_captain_order_names). Row R46 (the captain's order of 2026-11-24): the pair of T1's lines 8
# and 9, T1_s2rcB25MSH05, row R44's Blood step on T1_s2rcB25MS's command (row R45's), taken only with --captain-ruling
# R46 while row R46's first cell is the captain's and names the file. Both go through row R51's route (_t1_r51), whose
# record's shape they share.
T1_R60_ROW = "R60"
T1_R60_BASE_FILE = "T1_s2rcMS20"
T1_R60_BASE_OPTION = "--feng-fraction-mult 2.0 --direction-momentum all"
T1_R60_FILES = {  # file -> its option string, the multiplier, g_mom's damping, Blood's damping, MS's damping (as R51's)
    "T1_s2rcMS20d10": (f"{T1_R60_BASE_OPTION} --g-mom-damping 1.0", 2.0, 1.0, None, None),
    "T1_s2rcMS20H05": (f"{T1_R60_BASE_OPTION} --held-momentum 0.5", 2.0, None, "0.5", None),
    "T1_s2rcMS20d20": (f"{T1_R60_BASE_OPTION} --g-mom-damping 2.0", 2.0, 2.0, None, None),
    "T1_s2rcMS20m10": (f"{T1_R60_BASE_OPTION} --direction-momentum-damping 1.0", 2.0, None, None, "1.0"),
    "T1_s2rcMS20d15": (f"{T1_R60_BASE_OPTION} --g-mom-damping 1.5", 2.0, 1.5, None, None),
    "T1_s2rcMS20H10": (f"{T1_R60_BASE_OPTION} --held-momentum 1.0", 2.0, None, "1.0", None),
    "T1_s2rcMS20m20": (f"{T1_R60_BASE_OPTION} --direction-momentum-damping 2.0", 2.0, None, None, "2.0"),
}
T1_R60_COMMAND = T1_CANDIDATE_COMMAND + " {option} --captain-ruling " + T1_R60_ROW
T1_L12_BLOOD = ("Blood's own released E9.5 -> E10.5 change at damping {:.1f} on the drawn Blood cells, after MS's step "
                "and before g_mom (R1 reversed for Blood in this one file, Macrophage and Early VIC held)")
T1_R60_CHANGES = {
    "T1_s2rcMS20d10": "g_mom's within-type damping 1.0", "T1_s2rcMS20d15": "g_mom's within-type damping 1.5",
    "T1_s2rcMS20d20": "g_mom's within-type damping 2.0",
    "T1_s2rcMS20H05": T1_L12_BLOOD.format(0.5), "T1_s2rcMS20H10": T1_L12_BLOOD.format(1.0),
    "T1_s2rcMS20m10": "MS's own damping 1.0 in place of harness.B0_DAMPING",
    "T1_s2rcMS20m20": "MS's own damping 2.0 in place of harness.B0_DAMPING",
}
T1_R60_ORDER = ("experiments/harness/TRIES.md row R60, the captain's order of 2026-11-25 (\"ruling, no run\"; decision "
                "log 2026-11-25 (i), T1's line 12) naming `{option}` and `--captain-ruling R60`: file {file}, "
                "T1_s2rcMS20's command (row R51's file: candidate MS's step exactly as row R23 fixes it, the "
                "atlas step at 2.0 x the calendar fraction, g_mom's within-type damping 0.5, MS's damping "
                "harness.B0_DAMPING) with one change, {change}; nothing else changed")
T1_R46_ROW = "R46"
T1_R46_FILE = "T1_s2rcB25MSH05"
T1_R46_BASE_FILE = T1_MS25_FILE
T1_R46_OPTION = "--feng-fraction-mult 2.5 --direction-momentum all --held-momentum 0.5"
T1_R46_VALUES = (2.5, None, "0.5", None)  # the multiplier, g_mom's damping, Blood's damping, MS's damping, as given
T1_R46_COMMAND = f"{T1_CANDIDATE_COMMAND} {T1_R46_OPTION} --captain-ruling {T1_R46_ROW}"
T1_R46_ORDER = ("experiments/harness/TRIES.md row R46, the captain's order of 2026-11-24 (\"ruling, no run\"; decision "
                "log 2026-11-24 (ss)) naming T1_s2rcB25MSH05, the pair of T1's lines 8 and 9: T1_s2rcB25MS's command "
                "(row R45's: candidate MS's step exactly as row R23 fixes it, the atlas step at 2.5 x the calendar "
                "fraction, g_mom's within-type damping 0.5) with row R44's change, " + T1_L12_BLOOD.format(0.5)
                + "; nothing else changed")
# Row R66 (the captain's order of 2026-11-25, ruling, no run; decision log 2026-11-25 (hh)): T1's third-round step-on
# slots on T1_s2rcMS20's command (row R51's file, the base of row R60's), thirty files, each taken only with
# --captain-ruling R66 while row R66's first cell is the captain's order naming the file's option string
# (_captain_order_names), as row R60's files are taken: for K1, K2 and K3, two or three of line 12's steps together
# (g_mom's damping d, Blood's step H, MS's own damping m); for C1, g_mom's damping 2.5; for C3, Blood's step at 2.0.
# Blood's step goes as in row R60's H files (R1 reversed for Blood and T64's close lifted in each file with
# --held-momentum). The one table: (file, option string) in the row's own order, each name T1_s2rcMS20 plus the
# suffixes in d, H, m order; every value the builder routes on is read off the option string (_r66_values), never
# written a second time, so the row's strings, the names and the routes cannot drift apart.
T1_R66_ROW = "R66"
T1_R66_BASE_FILE = T1_R60_BASE_FILE
T1_R66_TABLE = (
    ("T1_s2rcMS20d10H05", f"{T1_R60_BASE_OPTION} --g-mom-damping 1.0 --held-momentum 0.5"),
    ("T1_s2rcMS20d10H10", f"{T1_R60_BASE_OPTION} --g-mom-damping 1.0 --held-momentum 1.0"),
    ("T1_s2rcMS20d15H05", f"{T1_R60_BASE_OPTION} --g-mom-damping 1.5 --held-momentum 0.5"),
    ("T1_s2rcMS20d15H10", f"{T1_R60_BASE_OPTION} --g-mom-damping 1.5 --held-momentum 1.0"),
    ("T1_s2rcMS20d20H05", f"{T1_R60_BASE_OPTION} --g-mom-damping 2.0 --held-momentum 0.5"),
    ("T1_s2rcMS20d20H10", f"{T1_R60_BASE_OPTION} --g-mom-damping 2.0 --held-momentum 1.0"),
    ("T1_s2rcMS20d10m10", f"{T1_R60_BASE_OPTION} --g-mom-damping 1.0 --direction-momentum-damping 1.0"),
    ("T1_s2rcMS20d10m20", f"{T1_R60_BASE_OPTION} --g-mom-damping 1.0 --direction-momentum-damping 2.0"),
    ("T1_s2rcMS20d15m10", f"{T1_R60_BASE_OPTION} --g-mom-damping 1.5 --direction-momentum-damping 1.0"),
    ("T1_s2rcMS20d15m20", f"{T1_R60_BASE_OPTION} --g-mom-damping 1.5 --direction-momentum-damping 2.0"),
    ("T1_s2rcMS20d20m10", f"{T1_R60_BASE_OPTION} --g-mom-damping 2.0 --direction-momentum-damping 1.0"),
    ("T1_s2rcMS20d20m20", f"{T1_R60_BASE_OPTION} --g-mom-damping 2.0 --direction-momentum-damping 2.0"),
    ("T1_s2rcMS20H05m10", f"{T1_R60_BASE_OPTION} --held-momentum 0.5 --direction-momentum-damping 1.0"),
    ("T1_s2rcMS20H05m20", f"{T1_R60_BASE_OPTION} --held-momentum 0.5 --direction-momentum-damping 2.0"),
    ("T1_s2rcMS20H10m10", f"{T1_R60_BASE_OPTION} --held-momentum 1.0 --direction-momentum-damping 1.0"),
    ("T1_s2rcMS20H10m20", f"{T1_R60_BASE_OPTION} --held-momentum 1.0 --direction-momentum-damping 2.0"),
    ("T1_s2rcMS20d10H05m10",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 1.0 --held-momentum 0.5 --direction-momentum-damping 1.0"),
    ("T1_s2rcMS20d10H05m20",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 1.0 --held-momentum 0.5 --direction-momentum-damping 2.0"),
    ("T1_s2rcMS20d10H10m10",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 1.0 --held-momentum 1.0 --direction-momentum-damping 1.0"),
    ("T1_s2rcMS20d10H10m20",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 1.0 --held-momentum 1.0 --direction-momentum-damping 2.0"),
    ("T1_s2rcMS20d15H05m10",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 1.5 --held-momentum 0.5 --direction-momentum-damping 1.0"),
    ("T1_s2rcMS20d15H05m20",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 1.5 --held-momentum 0.5 --direction-momentum-damping 2.0"),
    ("T1_s2rcMS20d15H10m10",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 1.5 --held-momentum 1.0 --direction-momentum-damping 1.0"),
    ("T1_s2rcMS20d15H10m20",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 1.5 --held-momentum 1.0 --direction-momentum-damping 2.0"),
    ("T1_s2rcMS20d20H05m10",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 2.0 --held-momentum 0.5 --direction-momentum-damping 1.0"),
    ("T1_s2rcMS20d20H05m20",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 2.0 --held-momentum 0.5 --direction-momentum-damping 2.0"),
    ("T1_s2rcMS20d20H10m10",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 2.0 --held-momentum 1.0 --direction-momentum-damping 1.0"),
    ("T1_s2rcMS20d20H10m20",
     f"{T1_R60_BASE_OPTION} --g-mom-damping 2.0 --held-momentum 1.0 --direction-momentum-damping 2.0"),
    ("T1_s2rcMS20d25", f"{T1_R60_BASE_OPTION} --g-mom-damping 2.5"),
    ("T1_s2rcMS20H20", f"{T1_R60_BASE_OPTION} --held-momentum 2.0"),
)
# the steps a row R66 file can change, in d, H, m order: (the option's flag, the name's letter)
T1_R66_STEPS = (("--g-mom-damping", "d"), ("--held-momentum", "H"), ("--direction-momentum-damping", "m"))


def _r66_values(option: str, base_option: str = T1_R60_BASE_OPTION, row: str = T1_R66_ROW) -> tuple:
    """A row R66 option string's values in T1_R60_FILES's order and types, as the builder's parser gives them: the
    multiplier, g_mom's damping (floats), Blood's damping and MS's damping (the literal strings), None where the string
    does not give the option. Raises ValueError unless the string is T1_R60_BASE_OPTION followed by one to three of
    T1_R66_STEPS' flags, each once with one value, in d, H, m order."""
    # base_option, row: the strings' head and the row a refusal names (row R73's four strings start from MSB's
    # pair, row R26's; the defaults are row R66's, whose table reads as before)
    words, base = option.split(), base_option.split()
    flags = [flag for flag, _ in T1_R66_STEPS]
    tail = words[len(base):]
    given = tail[::2]
    if (words[:len(base)] != base or not tail or len(tail) % 2 or any(f not in flags for f in given)
            or [flags.index(f) for f in given] != sorted({flags.index(f) for f in given})):
        raise ValueError(f"not a row {row} option string: {option}")
    values = dict(zip(given, tail[1::2], strict=True))
    damping = values.get("--g-mom-damping")
    return (float(base[1]), None if damping is None else float(damping), values.get("--held-momentum"),
            values.get("--direction-momentum-damping"))


# file -> its option string, the multiplier, g_mom's damping, Blood's damping, MS's damping (T1_R60_FILES's shape)
T1_R66_FILES = {name: (option, *_r66_values(option)) for name, option in T1_R66_TABLE}
T1_R66_COMMAND = T1_CANDIDATE_COMMAND + " {option} --captain-ruling " + T1_R66_ROW
# file -> its changes to T1_s2rcMS20's command in d, H, m order, each in row R60's words for that step
T1_R66_CHANGES = {
    name: tuple(text for value, text in (
        (damping, None if damping is None else f"g_mom's within-type damping {damping:.1f}"),
        (held, None if held is None else T1_L12_BLOOD.format(float(held))),
        (ms_damping, f"MS's own damping {ms_damping} in place of harness.B0_DAMPING")) if value is not None)
    for name, (_, _, damping, held, ms_damping) in T1_R66_FILES.items()}
# file -> its slot in the round's design as row R66 groups them: two or three changes are K1, K2 and K3's, the one
# change of g_mom's damping C1's and the one change of Blood's step C3's
T1_R66_SLOTS = {name: ("slots K1, K2 and K3" if len(changes) > 1 else
                       "slot C1" if T1_R66_FILES[name][2] is not None else "slot C3")
                for name, changes in T1_R66_CHANGES.items()}
T1_R66_COUNT = {1: "one change", 2: "two changes", 3: "three changes"}
T1_R66_ORDER = ("experiments/harness/TRIES.md row R66, the captain's order of 2026-11-25 (\"ruling, no run\"; decision "
                "log 2026-11-25 (hh), T1's third-round step-on {slot}) naming `{option}` and `--captain-ruling R66`: "
                "file {file}, T1_s2rcMS20's command (row R51's file: candidate MS's step exactly as row R23 fixes it, "
                "the atlas step at 2.0 x the calendar fraction, g_mom's within-type damping 0.5, MS's damping "
                "harness.B0_DAMPING) with {changes}; nothing else changed{blood}")
T1_R66_BLOOD = ("; R1 reversed for Blood and T64's close lifted in this file, as row R66 orders for each file with a "
                "Blood step and as row R60 did in its H files")


def _r66_order(name: str) -> str:
    """T1_R66_ORDER for a row R66 file: its slot, option string, name, its changes counted and listed, and Blood's
    ruling where the file has Blood's step."""
    changes = T1_R66_CHANGES[name]
    listed = changes[0] if len(changes) == 1 else ", ".join(changes[:-1]) + " and " + changes[-1]
    return T1_R66_ORDER.format(slot=T1_R66_SLOTS[name], option=T1_R66_FILES[name][0], file=name,
                               changes=f"{T1_R66_COUNT[len(changes)]}, {listed}",
                               blood=T1_R66_BLOOD if T1_R66_FILES[name][3] is not None else "")


# Row R73 (the captain's order of 2026-11-25, "All four", ruling, no run): T1's second refinement round, four files,
# each T1_s2rcMSBm10's command (row R51's file: T1_s2rcMSB's, row R26's, with MS's own damping 1.0) with one further
# change, g_mom's within-type damping (d) or Blood's own step (H), each taken only with --captain-ruling R73. Row R66's
# route (_t1_l12) on MSB's pair: _r66_values reads each string with MSB's pair as its head, and _r73_values checks each
# name against its string (T1_s2rcMSBm10, then the step's tag and its value's digits, the point dropped). This table, in
# the row's order, is the one the row's first cell is checked against and the one experiments/builds/t1_r73_line.sh's
# literal case table copies (tests/test_t1_r73_launcher.py holds the two equal).
T1_R73_ROW = "R73"
T1_R73_BASE_FILE = T1_R51_M10_FILE  # T1_s2rcMSBm10, row R51's file
T1_R73_BASE_OPTION = "--feng-fraction-mult 1.5 --direction-momentum all"  # MSB's pair (row R26's), each string's head
T1_R73_TABLE = (
    ("T1_s2rcMSBm10d10", f"{T1_R73_BASE_OPTION} --g-mom-damping 1.0 --direction-momentum-damping 1.0"),
    ("T1_s2rcMSBm10H05", f"{T1_R73_BASE_OPTION} --held-momentum 0.5 --direction-momentum-damping 1.0"),
    ("T1_s2rcMSBm10d15", f"{T1_R73_BASE_OPTION} --g-mom-damping 1.5 --direction-momentum-damping 1.0"),
    ("T1_s2rcMSBm10H10", f"{T1_R73_BASE_OPTION} --held-momentum 1.0 --direction-momentum-damping 1.0"),
)
T1_R73_STEPS = {"--g-mom-damping": "d", "--held-momentum": "H"}  # the one further step's flag -> its tag in the name


def _r73_values(name: str, option: str) -> tuple:
    """(multiplier, g_mom's damping or None, Blood's damping or None, MS's damping) of one entry of row R73's table, as
    _r66_values reads them off the string with MSB's pair as its head. Raises ValueError unless the string is MSB's
    pair, exactly one further step (--g-mom-damping or --held-momentum, its value one digit, a point and one digit) and
    --direction-momentum-damping 1.0, in that order, and the name is T1_s2rcMSBm10 plus the step's tag and its value's
    digits with the point dropped (T1_s2rcMSBm10d10 for --g-mom-damping 1.0, T1_s2rcMSBm10H05 for --held-momentum
    0.5): R66's d, H, m order would name the first T1_s2rcMSBd10m10, which this table does not take."""
    values = _r66_values(option, T1_R73_BASE_OPTION, T1_R73_ROW)
    tail = option.split()[len(T1_R73_BASE_OPTION.split()):]
    if (len(tail) != 4 or tail[0] not in T1_R73_STEPS or not re.fullmatch(r"\d\.\d", tail[1])
            or tail[2:] != ["--direction-momentum-damping", "1.0"]
            or name != f"{T1_R73_BASE_FILE}{T1_R73_STEPS[tail[0]]}{tail[1].replace('.', '')}"):
        raise ValueError(f"row {T1_R73_ROW}'s file {name} does not match its option string {option}")
    return values


T1_R73_FILES = {name: (option, *_r73_values(name, option)) for name, option in T1_R73_TABLE}
T1_R73_COMMAND = T1_CANDIDATE_COMMAND + " {option} --captain-ruling " + T1_R73_ROW
# each file's one further step is row R51's single-step file at the same value (T1_s2rcMSB plus the step's suffix): the
# launcher rebuilds that file, and the file's order repeats its change in row R51's words
T1_R73_SINGLE = {name: T1_R51_BASE_FILE + name.removeprefix(T1_R73_BASE_FILE) for name in T1_R73_FILES}
T1_R73_CHANGES = {name: T1_R51_CHANGES[single] for name, single in T1_R73_SINGLE.items()}
T1_R73_ORDER = ("experiments/harness/TRIES.md row R73, the captain's order of 2026-11-25 (\"All four\", the four files "
                "sent at once; \"ruling, no run\"; decision log 2026-11-25 (bbbb), T1's second refinement round) "
                "naming `{option}` ({file}) and `--captain-ruling R73`: file {file}, T1_s2rcMSBm10's command (row "
                "R51's file: candidate MS's step exactly as row R23 fixes it on candidate B's command, the atlas step "
                "at 1.5 x the calendar fraction, g_mom's within-type damping 0.5, MS's own damping 1.0 in place of "
                "harness.B0_DAMPING) with one further change, {change}, the change row R51's file {single} makes to "
                "T1_s2rcMSB's command; nothing else changed")


def _r73_row_names(name: str, path: Path | None = None) -> bool:
    """True when the ledger holds exactly one row R73 whose first cell starts "R73 (the captain's order", says "ruling,
    no run" and names the file's option string and `--captain-ruling R73` each as a whole backticked string, the
    string followed there by " (<file>)" as the row pairs each string with its file (row R67's check)."""
    lines = [ln for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines()
             if ln.startswith(f"| {T1_R73_ROW} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    option = T1_R73_FILES[name][0]
    named = set(re.findall(r"`([^`]*)`", first))
    return (first.startswith(f"{T1_R73_ROW} (the captain's order") and "ruling, no run" in first and option in named
            and f"--captain-ruling {T1_R73_ROW}" in named and f"`{option}` ({name})" in first)


def _r46_row_names(path: Path | None = None) -> bool:
    """Whether experiments/harness/TRIES.md (TRIES_FILE) holds exactly one row R46 whose first cell starts "R46 (the
    captain's ", says "ruling, no run" and names the file T1_s2rcB25MSH05 as a whole word (the row names the file and
    "if both beat", no option string). A name outside the first cell, or inside a longer name, does not pass."""
    lines = [ln for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines()
             if ln.startswith(f"| {T1_R46_ROW} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    return (first.startswith(f"{T1_R46_ROW} (the captain's ") and "ruling, no run" in first
            and re.search(rf"(?<![A-Za-z0-9_]){T1_R46_FILE}(?![A-Za-z0-9_])", first) is not None)


def _l12_hit(a) -> tuple[str, str] | None:
    """(row, file) when the command gives exactly one of row R66's thirty option sets, row R60's seven or row R46's
    (--direction-momentum all, and --feng-fraction-mult, --g-mom-damping, --held-momentum and
    --direction-momentum-damping each at the file's value or not given, as its option string has them), whatever
    --captain-ruling says; else None. The three rows' sets are disjoint, so the order they are tried in changes
    nothing. Reads nothing."""
    if getattr(a, "direction_momentum", None) not in DIRECTION_MOMENTUM_ARMS:
        return None
    given = (getattr(a, "feng_fraction_mult", None), getattr(a, "g_mom_damping", None),
             getattr(a, "held_momentum", None), getattr(a, "direction_momentum_damping", None))
    for row, files in ((T1_R73_ROW, T1_R73_FILES), (T1_R66_ROW, T1_R66_FILES), (T1_R60_ROW, T1_R60_FILES)):
        for name, (_, *values) in files.items():
            if given == tuple(values):
                return row, name
    return (T1_R46_ROW, T1_R46_FILE) if given == T1_R46_VALUES else None


def _l12_row_names(row: str, name: str) -> bool:
    """Row R60's or row R66's check for its file (_captain_order_names on the file's backticked option string) or row
    R46's (_r46_row_names)."""
    if row == T1_R60_ROW:
        return _captain_order_names(T1_R60_ROW, f"`{T1_R60_FILES[name][0]}`")
    if row == T1_R66_ROW:
        return _captain_order_names(T1_R66_ROW, f"`{T1_R66_FILES[name][0]}`")
    if row == T1_R73_ROW:
        return _r73_row_names(name)
    return _r46_row_names()


def _l12_ordered(a) -> bool:
    """The commands rows R60, R46 and R66 order, as _t1_l12 takes them: a file's option set exactly with its own row's
    --captain-ruling while that row names it, and no released bet's option."""
    ruling = getattr(a, "captain_ruling", None)
    if ruling not in (T1_R60_ROW, T1_R46_ROW, T1_R66_ROW, T1_R73_ROW):
        return False
    hit = _l12_hit(a)
    return (hit is not None and hit[0] == ruling and getattr(a, "released_bet", None) is None
            and getattr(a, "bet_strength", None) is None and _l12_row_names(*hit))


def _t1_l12(a) -> dict | None:
    """Row R60's seven files (T1's line 12, each one change to T1_s2rcMS20's command) and row R46's pair
    T1_s2rcB25MSH05 (row R44's Blood step on T1_s2rcB25MS's command), taken by _t1_r51 before row R51's own files: the
    record's values have _t1_r51's shape, so every step and the record go as row R51's files go. Returns None, reading
    nothing, unless the command gives exactly one of the eight option sets with its own row's --captain-ruling (R60
    for T1_R60_FILES's, R46 for T1_R46_OPTION): every other command, a file's set under any other ruling or none and
    any other set under R60 or R46 included, is routed as before, where each is refused. A file's set under its own
    row is refused before anything but the ledger is read unless the row names it (R60: its first cell the captain's
    order naming the option string and `--captain-ruling R60`, _captain_order_names; R46: its first cell starting "R46
    (the captain's ", saying "ruling, no run" and naming T1_s2rcB25MSH05, _r46_row_names), no --released-bet or
    --bet-strength is given, R1's hold stands, and the command is T1_CANDIDATE_COMMAND (whose own refusals then hold
    the rest) at the default seed on the test board (never a dry run).

    Row R66's thirty files (T1's third-round step-on slots, T1_R66_FILES) go the same way under --captain-ruling R66
    alone, each while row R66's first cell is the captain's order naming its option string and `--captain-ruling R66`
    (_captain_order_names), with R60's refusals; their record also names the slot and lists the file's changes."""
    ruling = getattr(a, "captain_ruling", None)
    if ruling not in (T1_R60_ROW, T1_R46_ROW, T1_R66_ROW, T1_R73_ROW):
        return None
    hit = _l12_hit(a)
    if hit is None or hit[0] != ruling:
        if ruling == T1_R73_ROW:  # row R73 orders its four sets and nothing else
            sys.exit(f"--captain-ruling {T1_R73_ROW} orders four files, each with exactly its own option string on "
                     f"{T1_CANDIDATE_COMMAND}: " + "; ".join(f"{f} `{v[0]}`" for f, v in T1_R73_FILES.items()))
        return None
    row, name = hit
    r60, r66, r73 = row == T1_R60_ROW, row == T1_R66_ROW, row == T1_R73_ROW
    if r73:
        option, mult, damping, held, ms_damping = T1_R73_FILES[name]
    elif r66:
        option, mult, damping, held, ms_damping = T1_R66_FILES[name]
    else:
        option, mult, damping, held, ms_damping = T1_R60_FILES[name] if r60 else (T1_R46_OPTION, *T1_R46_VALUES)
    base = T1_R73_BASE_FILE if r73 else T1_R66_BASE_FILE if r66 else T1_R60_BASE_FILE if r60 else T1_R46_BASE_FILE
    command = (T1_R73_COMMAND.format(option=option) if r73 else T1_R66_COMMAND.format(option=option) if r66
               else T1_R60_COMMAND.format(option=option) if r60 else T1_R46_COMMAND)
    what = f"row {row}'s file {name}"
    if not _l12_row_names(row, name):
        named = (f"the captain's order (\"ruling, no run\") naming `{option}` ({name}) and `--captain-ruling {row}`, "
                 "each as a whole backticked string" if r73 else
                 f"the captain's order (\"ruling, no run\") naming `{option}` and `--captain-ruling {row}`"
                 if r60 or r66 else
                 f"the captain's (starting \"{row} (the captain's \", \"ruling, no run\") naming {name}")
        sys.exit(f"{what} is refused: experiments/harness/TRIES.md has no single row {row} whose first cell is {named}")
    stacked = [flag for flag, attr in R51_STACKED if getattr(a, attr, None) is not None]
    if stacked:
        count = T1_R66_COUNT[len(T1_R66_CHANGES[name])] if r66 else "one change"
        sys.exit(f"{what} is refused with {', '.join(stacked)}: it is {base}'s command with {count}")
    if not getattr(a, "hold_unbracketed", True):
        sys.exit(f"{what} is refused with --no-hold-unbracketed: the cells without a group keep R1's hold")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"{what} changes {base} only: it needs the command {command} at the default seed on the test board")
    if r73:
        order = T1_R73_ORDER.format(option=option, file=name, change=T1_R73_CHANGES[name], single=T1_R73_SINGLE[name])
    elif r66:
        order = _r66_order(name)
    else:
        order = T1_R60_ORDER.format(option=option, file=name, change=T1_R60_CHANGES[name]) if r60 else T1_R46_ORDER
    out = {"direction_momentum": a.direction_momentum, "captain_ruling": row, "file": name, "base_file": base,
           "option": f"`{option}`", "feng_fraction_mult": float(mult),
           "g_mom_damping": G_MOM_DAMPINGS[0] if damping is None else float(damping),
           "held_momentum": None if held is None else float(held),
           "direction_momentum_damping": float(harness.B0_DAMPING) if ms_damping is None else float(ms_damping),
           "order": order, "command": command}
    if r66:
        out = {**out, "slot": T1_R66_SLOTS[name], "changes": list(T1_R66_CHANGES[name])}
    return out


def _l12_line(r51: dict) -> str:
    """g_mom's record text for a damping file of row R51's route: row R60's files name row R60, T1's line 12 and
    T1_s2rcMS20's command; row R66's name row R66, the file's slot and T1_s2rcMS20's command; every other file
    _r51_line's text."""
    if r51.get("captain_ruling") == T1_R60_ROW:
        return f"row {T1_R60_ROW}'s file {r51['file']}, T1's line 12, on {T1_R60_BASE_FILE}'s command"
    if r51.get("captain_ruling") == T1_R66_ROW:
        return (f"row {T1_R66_ROW}'s file {r51['file']}, T1's step-on {T1_R66_SLOTS[r51['file']]}, on "
                f"{T1_R66_BASE_FILE}'s command")
    if r51.get("captain_ruling") == T1_R73_ROW:
        return f"row {T1_R73_ROW}'s file {r51['file']}, T1's second refinement round, on {T1_R73_BASE_FILE}'s command"
    return _r51_line(r51)


def _l12_momentum(r51: dict, ms_record: dict | None) -> dict | None:
    """MS's record for a file of row R51's route: for row R60's files T1_s2rcMS20m10 and m20, and row R66's files with
    --direction-momentum-damping, the damping_option text names that file in place of row R51's (MS's step writes row
    R51's); every other record as it is."""
    row = r51.get("captain_ruling")
    if row not in (T1_R60_ROW, T1_R66_ROW, T1_R73_ROW) or not ms_record or "damping_option" not in ms_record:
        return ms_record
    rest = ("nothing else changed" if row == T1_R60_ROW else
            f"as row {T1_R51_ROW}'s file {T1_R73_BASE_FILE}, this file's base, has it; the file's one further change "
            "as its order names it" if row == T1_R73_ROW else
            f"the file's other changes as its order lists them ({T1_R66_COUNT[len(T1_R66_CHANGES[r51['file']])]} in "
            "all)")
    return {**ms_record, "damping_option": (f"--direction-momentum-damping {r51['direction_momentum_damping']:.1f} "
                                            f"(row {row}'s file {r51['file']}): MS's own damping in place of "
                                            f"harness.B0_DAMPING ({float(harness.B0_DAMPING):g}); {rest}")}


def _r51_row_names(option: str, path: Path | None = None) -> bool:
    """Whether experiments/harness/TRIES.md (TRIES_FILE) holds exactly one row R51 whose first cell is the captain's
    ruling (it starts "R51 (the captain's " and says "ruling, no run") naming `option` and `--captain-ruling R51`, each
    as a whole backticked string: a backticked string holding the option inside a longer one, or one the option holds,
    does not name it, and neither does a name outside the first cell."""
    lines = [ln for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines()
             if ln.startswith(f"| {T1_R51_ROW} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    named = set(re.findall(r"`([^`]*)`", first))
    return (first.startswith(f"{T1_R51_ROW} (the captain's ") and "ruling, no run" in first and option in named
            and f"--captain-ruling {T1_R51_ROW}" in named)


def _r51_file(a) -> str | None:
    """The row R51 file whose option set the command gives exactly (--direction-momentum all, and --feng-fraction-mult,
    --g-mom-damping, --held-momentum and --direction-momentum-damping each given at the file's value or not given, as
    the file's option string has them) with --captain-ruling R51, while row R51 names that option string
    (_r51_row_names); else None, and nothing is read unless the option set is one of the seven with R51."""
    if (getattr(a, "captain_ruling", None) != T1_R51_ROW
            or getattr(a, "direction_momentum", None) not in DIRECTION_MOMENTUM_ARMS):
        return None
    given = (getattr(a, "feng_fraction_mult", None), getattr(a, "g_mom_damping", None),
             getattr(a, "held_momentum", None), getattr(a, "direction_momentum_damping", None))
    for name, (option, *values) in T1_R51_FILES.items():
        if given == tuple(values):
            return name if _r51_row_names(option) else None
    return None


def _r51_ordered(a) -> bool:
    """The commands row R51 orders, as _t1_r51 holds them: one of its seven option sets exactly with --captain-ruling
    R51 while its row names that set, and no released bet's option."""
    return (_r51_file(a) is not None and getattr(a, "released_bet", None) is None
            and getattr(a, "bet_strength", None) is None)


def _t1_r51(a) -> dict | None:
    """Row R51's seven files (T1's line 10), each one change to T1_s2rcMSB's command, tested before every other route.
    Returns None when the command is not exactly one of them under its row (_r51_file; every such command is then
    routed as before), after refusing --direction-momentum-damping, which no other command takes. One of them is
    refused before anything but the ledger is read with --released-bet or --bet-strength, with --no-hold-unbracketed,
    and on any command but T1_CANDIDATE_COMMAND (whose own refusals then hold the rest) at the default seed on the test
    board (never a dry run). Returns the record's values: the file, the row, its option string, the base file, the four
    settings the build takes from here (the multiplier, g_mom's damping, Blood's damping or None, MS's damping), the
    order text and the command."""
    l12 = _t1_l12(a)  # rows R60 and R46 first, with their own refusals (None: every command as before)
    if l12 is not None:
        return l12
    name = _r51_file(a)
    if name is None:
        raw = getattr(a, "direction_momentum_damping", None)
        if raw is not None:
            m10 = T1_R51_FILES[T1_R51_M10_FILE][0]
            sys.exit(f"--direction-momentum-damping {raw} is refused: it takes {', '.join(T1_R51_DMD_CHOICES)} only, "
                     f"and only in row {T1_R51_ROW}'s file {T1_R51_M10_FILE}, the command "
                     f"{T1_R51_COMMAND.format(option=m10)} at the default seed on the test board, while "
                     f"experiments/harness/TRIES.md holds exactly one row {T1_R51_ROW} whose first cell is the "
                     f"captain's ruling (\"ruling, no run\") naming `{m10}` and `--captain-ruling {T1_R51_ROW}`")
        return None
    option, mult, damping, held, ms_damping = T1_R51_FILES[name]
    command = T1_R51_COMMAND.format(option=option)
    what = f"row {T1_R51_ROW}'s file {name}"
    stacked = [flag for flag, attr in R51_STACKED if getattr(a, attr, None) is not None]
    if stacked:
        sys.exit(f"{what} is refused with {', '.join(stacked)}: it is {T1_R51_BASE_FILE}'s command with one change")
    if not getattr(a, "hold_unbracketed", True):
        sys.exit(f"{what} is refused with --no-hold-unbracketed: the cells without a group keep R1's hold")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"{what} changes {T1_R51_BASE_FILE} only: it needs the command {command} at the default seed on the "
                 "test board")
    return {"direction_momentum": a.direction_momentum, "captain_ruling": T1_R51_ROW, "file": name,
            "base_file": T1_R51_BASE_FILE, "option": f"`{option}`", "feng_fraction_mult": float(mult),
            "g_mom_damping": G_MOM_DAMPINGS[0] if damping is None else float(damping),
            "held_momentum": None if held is None else float(held),
            "direction_momentum_damping": float(harness.B0_DAMPING) if ms_damping is None else float(ms_damping),
            "order": T1_R51_ORDER.format(option=option, file=name, change=T1_R51_CHANGES[name]), "command": command}


def _r51_line(r51: dict) -> str:
    """g_mom's record text for row R51's damping files: the row, the file and the base command."""
    return f"row {T1_R51_ROW}'s file {r51['file']}, T1's line 10, on {T1_R51_BASE_FILE}'s command"


def _r58_ordered(a) -> bool:
    """The command row R58 orders, as _t1_r58 holds it: its option set exactly with --captain-ruling R58, none of
    R58_STACKED's options, and row R58's first cell the captain's order naming the option string."""
    return (getattr(a, "captain_ruling", None) == T1_R58_ROW
            and getattr(a, "direction_momentum_damping", None) == T1_R58_DMD
            and getattr(a, "feng_fraction_mult", None) == T1_MSB_MULT
            and getattr(a, "direction_momentum", None) in DIRECTION_MOMENTUM_ARMS
            and all(getattr(a, attr, None) is None for _, attr in R58_STACKED)
            and _captain_order_names(T1_R58_ROW, T1_R58_OPTION))


def _t1_r58(a) -> dict | None:
    """Row R58's file T1_s2rcMSBm20 (T1's line 11), tested before row R51's route. Returns None unless the command
    names --captain-ruling R58: every other command is routed as before (--direction-momentum-damping 2.0 under any
    other ruling is refused by row R51's guard). Under R58 it is refused before anything but the ledger is read unless
    row R58's first cell is the captain's order naming T1_R58_OPTION and `--captain-ruling R58`, the option set is
    exactly T1_R58_OPTION_STRING's, none of R58_STACKED's options is given, R1's hold stands, and the command is
    T1_CANDIDATE_COMMAND (whose own refusals then hold the rest) at the default seed on the test board. Returns the
    record's values as row R51's route returns its files': the file, the row, its option string, the base file, the
    multiplier, g_mom's damping, Blood's (None), MS's damping 2.0, the order text and the command."""
    if getattr(a, "captain_ruling", None) != T1_R58_ROW:
        return None
    what = f"row {T1_R58_ROW}'s file {T1_R58_FILE}"
    if not _captain_order_names(T1_R58_ROW, T1_R58_OPTION):
        sys.exit(f"{what} is refused: experiments/harness/TRIES.md has no single row {T1_R58_ROW} whose first cell is "
                 f"the captain's order (\"ruling, no run\") naming {T1_R58_OPTION} and `--captain-ruling {T1_R58_ROW}`")
    given = (getattr(a, "feng_fraction_mult", None), getattr(a, "direction_momentum", None),
             getattr(a, "direction_momentum_damping", None))
    if given != (T1_MSB_MULT, DIRECTION_MOMENTUM_ARMS[0], T1_R58_DMD):
        sys.exit(f"--captain-ruling {T1_R58_ROW} orders one file, {T1_R58_FILE}: it needs exactly "
                 f"{T1_R58_OPTION_STRING} on {T1_CANDIDATE_COMMAND}")
    stacked = [flag for flag, attr in R58_STACKED if getattr(a, attr, None) is not None]
    if stacked:
        sys.exit(f"{what} is refused with {', '.join(stacked)}: it is {T1_R51_BASE_FILE}'s command with one change")
    if not getattr(a, "hold_unbracketed", True):
        sys.exit(f"{what} is refused with --no-hold-unbracketed: the cells without a group keep R1's hold")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"{what} changes {T1_R51_BASE_FILE} only: it needs the command {T1_R58_COMMAND} at the default seed "
                 "on the test board")
    return {"direction_momentum": a.direction_momentum, "captain_ruling": T1_R58_ROW, "file": T1_R58_FILE,
            "base_file": T1_R51_BASE_FILE, "option": T1_R58_OPTION, "feng_fraction_mult": float(T1_MSB_MULT),
            "g_mom_damping": G_MOM_DAMPINGS[0], "held_momentum": None,
            "direction_momentum_damping": float(T1_R58_DMD), "order": T1_R58_ORDER, "command": T1_R58_COMMAND}


def _r67_row_names(name: str, path: Path | None = None) -> bool:
    """Whether experiments/harness/TRIES.md (TRIES_FILE) holds exactly one row R67 whose first cell is the captain's
    order (it starts "R67 (the captain's order" and says "ruling, no run") naming row R67's file `name`'s option
    string and `--captain-ruling R67`, each as a whole backticked string, the option string followed there by
    " (<name>)" as the row pairs each string with its file. A backticked string holding the option inside a longer one,
    a string paired with another file, or a name outside the first cell does not pass."""
    lines = [ln for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines()
             if ln.startswith(f"| {T1_R67_ROW} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    option = T1_R67_FILES[name][0]
    named = set(re.findall(r"`([^`]*)`", first))
    return (first.startswith(f"{T1_R67_ROW} (the captain's order") and "ruling, no run" in first and option in named
            and f"--captain-ruling {T1_R67_ROW}" in named and f"`{option}` ({name})" in first)


def _r67_file(a) -> str | None:
    """The row R67 file whose option set the command gives exactly (--feng-fraction-mult 1.5, --direction-momentum all
    and --direction-momentum-damping at the file's literal value, T1_R67_FILES), whatever --captain-ruling says; else
    None. Reads nothing."""
    given = (getattr(a, "feng_fraction_mult", None), getattr(a, "direction_momentum", None),
             getattr(a, "direction_momentum_damping", None))
    for name, (_, value) in T1_R67_FILES.items():
        if given == (T1_MSB_MULT, DIRECTION_MOMENTUM_ARMS[0], value):
            return name
    return None


def _r67_ordered(a) -> bool:
    """The commands row R67 orders, as _t1_r67 holds them: one of its three option sets exactly with --captain-ruling
    R67, none of R67_STACKED's options, and row R67's first cell naming the set paired with its file. Reads nothing
    unless the command names --captain-ruling R67 and one of the three sets."""
    if getattr(a, "captain_ruling", None) != T1_R67_ROW:
        return False
    name = _r67_file(a)
    return (name is not None and all(getattr(a, attr, None) is None for _, attr in R67_STACKED)
            and _r67_row_names(name))


def _t1_r67(a) -> dict | None:
    """Row R67's three files (T1_R67_TABLE: T1_s2rcMSB's command with MS's own damping at 1.25, 1.5 or 0.75), on row
    R58's route, tested after row R58's and before row R51's. Returns None, reading nothing, unless the command names
    --captain-ruling R67: every other command is routed as before (the three values under any other ruling or none are
    refused by row R51's guard). Under R67 it is refused before anything but the ledger is read unless the option set
    is exactly one of the three (_r67_file), row R67's first cell names that option string paired with its file and
    `--captain-ruling R67` (_r67_row_names), none of R67_STACKED's options is given, R1's hold stands, and the command
    is T1_CANDIDATE_COMMAND (whose own refusals then hold the rest) at the default seed on the test board. Returns the
    record's values as row R58's route returns its file's: the file, the row, its option string, the base file, the
    multiplier, g_mom's damping, Blood's (None), MS's damping (the file's literal value), the order text and the
    command."""
    if getattr(a, "captain_ruling", None) != T1_R67_ROW:
        return None
    name = _r67_file(a)
    if name is None:
        sys.exit(f"--captain-ruling {T1_R67_ROW} orders three files, each with exactly its own option string on "
                 f"{T1_CANDIDATE_COMMAND}: " + "; ".join(f"{f} `{v[0]}`" for f, v in T1_R67_FILES.items()))
    option, value = T1_R67_FILES[name]
    command = T1_R67_COMMAND.format(option=option)
    what = f"row {T1_R67_ROW}'s file {name}"
    if not _r67_row_names(name):
        sys.exit(f"{what} is refused: experiments/harness/TRIES.md has no single row {T1_R67_ROW} whose first cell is "
                 f"the captain's order (\"ruling, no run\") naming `{option}` ({name}) and "
                 f"`--captain-ruling {T1_R67_ROW}`, each as a whole backticked string")
    stacked = [flag for flag, attr in R67_STACKED if getattr(a, attr, None) is not None]
    if stacked:
        sys.exit(f"{what} is refused with {', '.join(stacked)}: it is {T1_R51_BASE_FILE}'s command with one change")
    if not getattr(a, "hold_unbracketed", True):
        sys.exit(f"{what} is refused with --no-hold-unbracketed: the cells without a group keep R1's hold")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"{what} changes {T1_R51_BASE_FILE} only: it needs the command {command} at the default seed on the "
                 "test board")
    return {"direction_momentum": a.direction_momentum, "captain_ruling": T1_R67_ROW, "file": name,
            "base_file": T1_R51_BASE_FILE, "option": f"`{option}`", "feng_fraction_mult": float(T1_MSB_MULT),
            "g_mom_damping": G_MOM_DAMPINGS[0], "held_momentum": None, "direction_momentum_damping": float(value),
            "order": T1_R67_ORDER.format(option=option, file=name, value=value), "command": command}


def _r54_ordered(a) -> bool:
    """The command row R54 orders, as _t1_r54 holds it (and has already refused every other use of --late-mix by the
    time _bet_refusals asks): --late-mix with --captain-ruling R54, no released bet's option, and row R54's first cell
    the captain's order naming the option."""
    return (getattr(a, "late_mix", None) in LATE_MIX_CHOICES and getattr(a, "captain_ruling", None) == T1_R54_ROW
            and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None
            and _captain_order_names(T1_R54_ROW, T1_R54_OPTION))


def _late_mix_named_sha(path: Path | None = None) -> str | None:
    """The sha256 row N64 names for the shares file: the one `shares.json sha256 <hex>` string (LATE_MIX_SHA_NAMED) in
    the one row N64 of experiments/harness/TRIES.md (TRIES_FILE), in any of its cells (its outcome cell carries it once
    the step has run). None unless there is exactly one row N64 and exactly one distinct such string in it."""
    lines = [ln for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines()
             if ln.startswith(f"| {LATE_MIX_ROW} ")]
    if len(lines) != 1:
        return None
    named = set(LATE_MIX_SHA_NAMED.findall(lines[0]))
    return next(iter(named)) if len(named) == 1 else None


def _late_mix_bytes() -> tuple[bytes | None, str]:
    """The shares file's bytes as committed at HEAD (feng_relabel.committed_at_head, a `git show`), or None and why:
    not committed, or the tree's copy missing or not the same bytes."""
    import feng_relabel

    raw = feng_relabel.committed_at_head(LATE_MIX_SHARES)
    if raw is None:
        return None, f"the shares file {LATE_MIX_SHARES.name} of row {LATE_MIX_ROW} is not committed at HEAD"
    if not LATE_MIX_SHARES.is_file() or LATE_MIX_SHARES.read_bytes() != raw:
        return None, f"the tree's copy of the shares file {LATE_MIX_SHARES.name} is not the one committed at HEAD"
    return raw, ""


def _late_mix_values(raw: bytes, mix: str) -> dict:
    """The shares file's fields the build takes, checked: the row (N64), the anchor (--late-mix's value), the pool's
    own chamber mix it was computed from, the two moved shares (each strictly between 0 and 1) and the moved mix (every
    category of the chamber config, summing to 1). Refused otherwise."""
    try:
        rec = json.loads(raw)
        own = {str(k): float(v) for k, v in rec["own_shares"].items()}
        late = {str(k): float(v) for k, v in rec["shares"].items()}
        cm, blood = float(rec["cm_share_of_heart"]), float(rec["blood_share_of_chamber"])
        ok = (rec.get("row") == LATE_MIX_ROW and rec.get("late_mix") == mix and own and late
              and 0.0 < cm < 1.0 and 0.0 < blood < 1.0 and abs(sum(late.values()) - 1.0) <= 1e-9
              and abs(sum(own.values()) - 1.0) <= 1e-9)
    except (ValueError, KeyError, TypeError, AttributeError):
        ok = False
    if not ok:
        sys.exit(f"--late-mix {mix} is refused: the shares file is not row {LATE_MIX_ROW}'s step output for {mix} "
                 "(row, anchor, own mix, the two moved shares and the moved mix)")
    return {"own_shares": own, "cm_share_of_heart": cm, "blood_share_of_chamber": blood, "late_shares": late}


def _t1_r54(a) -> dict | None:
    """Row R54's file T1_s2rcMSB_LM16 (T1's late-mix line), tested before every other route. Returns None when
    --late-mix is not given (every command as before). Refused before anything but the ledger and the shares file is
    read: without --captain-ruling R54 or with any other row; unless row R54's first cell is the captain's order naming
    T1_R54_OPTION and `--captain-ruling R54`; on any command but T1_s2rcMSB's (--feng-fraction-mult 1.5
    --direction-momentum all on `t1 --g-mom --base lib3 --feng-groups recomputed`) with one of R54_STACKED's options
    or R54_SWITCHES, without R1's hold or R2, with --feng-none-ruling, at another seed or on a dry run; unless row N64
    names exactly one `shares.json sha256 <hex>`, the shares file committed at HEAD (the tree's copy the same bytes)
    has that sha256 and its fields are the step's. Returns the record's values: MSB's settings, the file, the row, the
    option, the base file, the shares file's path and sha256, the recorded own mix and the moved shares, the order
    text and the command."""
    mix = getattr(a, "late_mix", None)
    if mix is None:
        return None
    what = f"--late-mix {mix}"
    ruling = getattr(a, "captain_ruling", None)
    if ruling != T1_R54_ROW:
        sys.exit(f"{what} is row {T1_R54_ROW}'s file {T1_R54_FILE} alone and needs --captain-ruling {T1_R54_ROW}"
                 + ("" if ruling is None else f", not {ruling}"))
    if not _captain_order_names(T1_R54_ROW, T1_R54_OPTION):
        sys.exit(f"{what} is refused: experiments/harness/TRIES.md has no single row {T1_R54_ROW} whose first cell is "
                 f"the captain's ruling (\"ruling, no run\") naming {T1_R54_OPTION} and "
                 f"`--captain-ruling {T1_R54_ROW}`")
    stacked = [flag for flag, attr in R54_STACKED if getattr(a, attr, None) is not None]
    stacked += [flag for flag, attr in R54_SWITCHES if getattr(a, attr, False)]
    if stacked:
        sys.exit(f"{what} is refused with {', '.join(stacked)}: row {T1_R54_ROW}'s file is {T1_R51_BASE_FILE}'s "
                 "command with its chamber mix moved, nothing else changed")
    if (not getattr(a, "hold_unbracketed", True) or getattr(a, "feng_hygiene", None) is False
            or getattr(a, "feng_none_ruling", None) is not None):
        sys.exit(f"{what} is refused with --no-hold-unbracketed, --no-feng-hygiene or --feng-none-ruling: R1 and R2 "
                 f"stand as in {T1_R51_BASE_FILE}")
    if (getattr(a, "feng_fraction_mult", None) != T1_MSB_MULT
            or getattr(a, "direction_momentum", None) not in DIRECTION_MOMENTUM_ARMS or not getattr(a, "g_mom", False)
            or getattr(a, "base", None) != "lib3" or getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed"
            or a.seed != harness.BUILD_SEED_OFFSET + 1 or a.board is not None):
        sys.exit(f"{what} moves {T1_R51_BASE_FILE}'s chamber mix only: it needs the command {T1_R54_COMMAND} at the "
                 "default seed on the test board")
    sha = _late_mix_named_sha()
    if sha is None:
        sys.exit(f"{what} is refused: experiments/harness/TRIES.md has no single row {LATE_MIX_ROW} naming exactly one "
                 "`shares.json sha256 <hex>`")
    raw, why = _late_mix_bytes()
    if raw is None:
        sys.exit(f"{what} is refused: {why}")
    if hashlib.sha256(raw).hexdigest() != sha:
        sys.exit(f"{what} is refused: the committed shares file's sha256 is not the one row {LATE_MIX_ROW} names")
    try:
        shares_file = LATE_MIX_SHARES.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
    except ValueError:
        shares_file = str(LATE_MIX_SHARES)
    return {"direction_momentum": a.direction_momentum, "captain_ruling": T1_R54_ROW, "file": T1_R54_FILE,
            "base_file": T1_R51_BASE_FILE, "option": T1_R54_OPTION, "late_mix": mix,
            "feng_fraction_mult": float(T1_MSB_MULT), "g_mom_damping": G_MOM_DAMPINGS[0], "held_momentum": None,
            "direction_momentum_damping": float(harness.B0_DAMPING), "shares_file": shares_file,
            "shares_sha256": sha, "shares_row": LATE_MIX_ROW, **_late_mix_values(raw, mix), "order": T1_R54_ORDER,
            "command": T1_R54_COMMAND}


def _t1_late_mix_shares(own: dict, r54: dict) -> tuple[dict, dict]:
    """Row R54's shares in place of the pool's own: the pool's own chamber mix (`own`, harness.t1_chamber_shares as
    every lib3 file takes it) must be the one row N64's step recorded (LATE_MIX_TOL); it is then moved by
    t1_scope.with_share exactly as --cm-share and --blood-share move a mix (cardiomyocytes of heart cells, then blood
    of chamber cells, the other categories in proportion), and the result must be the step's recorded shares
    (LATE_MIX_TOL). Returns the shares and the record's late-mix part (shares and the file's hash; no count)."""
    rec_own = r54["own_shares"]
    if set(own) != set(rec_own) or any(abs(float(own[c]) - rec_own[c]) > LATE_MIX_TOL for c in own):
        sys.exit(f"--late-mix is refused: the pool's own chamber mix is not the one row {LATE_MIX_ROW}'s step recorded")
    moved = t1_scope.with_share(own, "cm", r54["cm_share_of_heart"])
    moved = t1_scope.with_share(moved, "blood", r54["blood_share_of_chamber"], within=list(moved))
    want = r54["late_shares"]
    if set(moved) != set(want) or any(abs(moved[c] - want[c]) > LATE_MIX_TOL for c in moved):
        sys.exit(f"--late-mix is refused: the moved shares are not row {LATE_MIX_ROW}'s recorded shares")
    heart = t1_scope.load_spec().heart
    return moved, {"row": LATE_MIX_ROW, "anchor": r54["late_mix"], "shares_file": r54["shares_file"],
                   "shares_sha256": r54["shares_sha256"], "own_shares": dict(own),
                   "own_cm_share_of_heart": t1_scope.share_of(own, "cm", heart),
                   "own_blood_share_of_chamber": float(own.get("blood", 0.0)),
                   "cm_share_of_heart": r54["cm_share_of_heart"],
                   "blood_share_of_chamber": r54["blood_share_of_chamber"],
                   "rule": "library 3's own chamber mix moved by vec.t1_scope.with_share, cardiomyocytes of heart "
                           "cells first, then blood of chamber cells, the other categories in proportion (the order "
                           "--cm-share and --blood-share take); quotas at the count rule, never weights"}


def _t1_ms_step(ref, library: np.ndarray, directions, fraction: float, groups, held, seed: int,
                damping: float | None = None, option_note: str | None = None):
    """Candidate MS's step on the drawn population, before g_mom: on the cells whose label's group has a direction
    (build_submission's covered test, `groups.get(label) in directions`), B0's per-category damped mean momentum
    (`t1_scope.momentum_shift` at harness.B0_DAMPING with its defaults: eps 0.5, min_cells 20, keep_depth; the held
    genes `held`; stream [seed, T1_MS_STREAM]; the chamber config's E9.5 and E10.5 stand-ins), fitted from E9.5 to the
    reference's portion libraries T1_MS_FIT_LIBRARIES (`library`: each reference row's library), then the atlas step on
    the same cells exactly as build_submission takes it (`t1_scope.external_shift`, the same directions, fraction,
    groups and held genes, stream [seed, T1_MS_DIRECTION_STREAM]); every other cell as drawn; stacked by
    `t1_scope._stack_rows`, covered cells first, as build_submission stacks them. The same calls in the same order as
    experiments/t1_m1/ms_read.py's `ms_population` (the practice arm m1rc_gmom_mom), whose fit runs on the practice
    task's own stages. Refused unless the reference's libraries other than the dissected one are exactly
    T1_MS_FIT_LIBRARIES. E9.5 is loaded here. Returns the step and its record (constants and stages; no cell number).
    damping (row R51's files only; None: harness.B0_DAMPING, the step and record exactly as before): the momentum's
    damping, passed to momentum_shift in place of harness.B0_DAMPING (file T1_s2rcMSBm10's 1.0; every other R51 file
    passes B0's value), nothing else changed; the record's damping field names it and, when it is not B0's, a
    damping_option field says so: `option_note` where it is given (row R58's file T1_s2rcMSBm20's 2.0, T1_R58_NOTE;
    row R67's three files, T1_R67_NOTES), else row R51's text as before."""
    portions = sorted(set(library.tolist()) - {T1_DISSECTED_LIBRARY})
    if portions != list(T1_MS_FIT_LIBRARIES):
        sys.exit(f"--direction-momentum is refused: E10.5's libraries other than {T1_DISSECTED_LIBRARY} are "
                 f"{portions}, not the fit's {list(T1_MS_FIT_LIBRARIES)} (row {T1_R23_ROW})")
    ms_damping = harness.B0_DAMPING if damping is None else float(damping)
    spec = t1_scope.load_spec()
    start = harness._shallow(load_stage(DATA / "T1/E9.5_RNA.h5ad", lineage=False))
    later = ref[np.flatnonzero(np.isin(library, T1_MS_FIT_LIBRARIES))].copy()
    t1_scope.classify(start, source="E9.5_RNA", spec=spec)
    t1_scope.classify(later, source="E10.5_RNA", spec=spec)
    kw = {"start_fallback": spec.fallbacks("E9.5_RNA"), "ref_fallback": spec.fallbacks("E10.5_RNA")}
    held = list(held)

    def step(pop):
        labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
        cov = np.array([groups.get(lab) in directions for lab in labels], dtype=bool)
        parts = []
        if cov.any():
            idx = np.flatnonzero(cov)
            mom, _ = t1_scope.momentum_shift(start, later, pop[idx].copy(), ms_damping, held,
                                             label_key=t1_scope.LABEL_KEY, seed=[seed, T1_MS_STREAM], **kw)
            moved, _ = t1_scope.external_shift(mom, directions, fraction, groups, held, label_key=t1_scope.LABEL_KEY,
                                               seed=[seed, T1_MS_DIRECTION_STREAM])
            parts.append((idx, moved))
        if not cov.all():
            idx = np.flatnonzero(~cov)
            parts.append((idx, pop[idx].copy()))
        return t1_scope._stack_rows(parts, pop.n_obs)

    record = {"rule": "B0's per-category damped mean momentum (vec.t1_scope.momentum_shift: log f = damping x "
                      "log((m_later + eps) / (m_start + eps)) per category over its own cells; f = 1 on the held genes "
                      "and for a category with fewer than min_cells own cells or drawn from stand-ins at either stage) "
                      "on the drawn cells whose label's group has a recomputed direction, then the atlas step on the "
                      "same cells, then g_mom; every other cell as drawn (R1's hold)",
              "damping": float(ms_damping), "eps": 0.5, "min_cells": 20, "keep_depth": True,
              "stream": [int(seed), T1_MS_STREAM], "direction_stream": [int(seed), T1_MS_DIRECTION_STREAM],
              "fit": {"start": "E9.5_RNA", "later": "E10.5_RNA", "later_libraries": list(T1_MS_FIT_LIBRARIES),
                      "library_3_in_fit": False},
              "held": "the build's held list: R2's union joined with the pre-registered list",
              "covered": "build_submission's covered cells: the label's group has a recomputed direction "
                         f"({sorted(directions)})",
              "practice_arm": "experiments/t1_m1/ms_read.py m1rc_gmom_mom (row T65), whose fit runs on the practice "
                              "task's own stages"}
    if damping is not None and ms_damping != harness.B0_DAMPING:  # R51's T1_s2rcMSBm10, R58's MSBm20, R67's three
        record["damping_option"] = (option_note if option_note is not None else
                                    f"--direction-momentum-damping {ms_damping:.1f} (row {T1_R51_ROW}'s file "
                                    f"{T1_R51_M10_FILE}): MS's own damping in place of harness.B0_DAMPING "
                                    f"({float(harness.B0_DAMPING):g}); nothing else changed")
    return step, record


def _t1_hm_step(ref, library: np.ndarray, damping: float, held, seed: int, base_file: str | None = None,
                row: str | None = None):
    """Blood's released momentum (--held-momentum; files T1_s2rcB25H05, H10 and H20), applied to the drawn population
    after the atlas step and before g_mom: on the drawn cells labelled T1_HM_LABEL only, B0's per-category damped mean
    momentum (`t1_scope.momentum_shift` with _t1_ms_step's constants: eps 0.5, min_cells 20, keep_depth; the held genes
    `held`, the build's held list, at f = 1; the chamber config's E9.5 and E10.5 stand-ins as _t1_ms_step passes them)
    fitted on Blood-labelled cells at both ends (within={"blood": ["Blood"]}) from E9.5 to the reference's portion
    libraries T1_MS_FIT_LIBRARIES (`library`: each reference row's library; library 3 never in the fit), at `damping`,
    on stream [seed, T1_HM_STREAM]; every other cell as it comes (the banked file's: the atlas step's output, R1's hold
    on Macrophage, Early VIC and the rest); stacked by `t1_scope._stack_rows`, Blood first, back in row order. Refused,
    before the draw (released stages only), unless the reference's libraries other than the dissected one are exactly
    T1_MS_FIT_LIBRARIES and Blood's own cells reach min_cells at both ends of the fit, so that the step is never idle.
    E9.5 is loaded here. Returns the step and its record (constants, stages, label and the fit's status; no cell
    number), which the step completes when it runs; the step stops the build unless it shifted the Blood cells.
    base_file (row R51's files T1_s2rcMSBH05 and H10; None: the record exactly as before): the same step on that file's
    command, after MS's step (whose momentum and atlas step move the cells with a direction, never Blood) and before
    g_mom; only the record's rule text names that file and MS's step in place of the banked file and the atlas step."""
    portions = sorted(set(library.tolist()) - {T1_DISSECTED_LIBRARY})
    if portions != list(T1_MS_FIT_LIBRARIES):
        sys.exit(f"--held-momentum is refused: E10.5's libraries other than {T1_DISSECTED_LIBRARY} are {portions}, not "
                 f"the fit's {list(T1_MS_FIT_LIBRARIES)}")
    spec = t1_scope.load_spec()
    start = harness._shallow(load_stage(DATA / "T1/E9.5_RNA.h5ad", lineage=False))
    later = ref[np.flatnonzero(np.isin(library, T1_MS_FIT_LIBRARIES))].copy()
    t1_scope.classify(start, source="E9.5_RNA", spec=spec)
    t1_scope.classify(later, source="E10.5_RNA", spec=spec)
    kw = {"start_fallback": spec.fallbacks("E9.5_RNA"), "ref_fallback": spec.fallbacks("E10.5_RNA")}
    held = list(held)
    within = {T1_HM_CATEGORY: [T1_HM_LABEL]}
    how = {}
    for end, stage, fallback in (("start", start, kw["start_fallback"]), ("later", later, kw["ref_fallback"])):
        rows, how[end] = t1_scope._type_rows(stage, T1_HM_CATEGORY, t1_scope.CATEGORY_KEY, t1_scope.LABEL_KEY, within,
                                             fallback, T1_HM_MIN_CELLS)
        if how[end] != "own, restricted" or len(rows) < T1_HM_MIN_CELLS:
            sys.exit(f"--held-momentum is refused: Blood's own cells at the fit's {end} stage do not reach min_cells "
                     "(the step would be idle)")
    record = {"rule": "Blood's own released change: B0's per-category damped mean momentum (vec.t1_scope."
                      "momentum_shift: log f = damping x log((m_later + eps) / (m_start + eps)) over the "
                      "Blood-labelled cells of category blood at each end; f = 1 on the held genes) on the drawn cells "
                      "labelled Blood, after the atlas step and before g_mom, which leaves Blood alone; every other "
                      "cell as the banked file has it (Macrophage and Early VIC keep R1's hold)",
              "label": T1_HM_LABEL, "category": T1_HM_CATEGORY, "damping": float(damping), "eps": T1_HM_EPS,
              "min_cells": T1_HM_MIN_CELLS, "keep_depth": True, "stream": [int(seed), T1_HM_STREAM],
              "fit": {"start": "E9.5_RNA", "later": "E10.5_RNA", "later_libraries": list(T1_MS_FIT_LIBRARIES),
                      "library_3_in_fit": False, "within": {T1_HM_CATEGORY: [T1_HM_LABEL]},
                      "start_cells": how["start"], "later_cells": how["later"]},
              "status": None,
              "held": "the build's held list: R2's union joined with the pre-registered list",
              "practice_arm": "experiments/t1_m1/heldmom_read.py m1rc_gmom_heldmom (row T64, candidate D) narrowed to "
                              "Blood, its fit on Blood-labelled cells at both ends; the fit's stages the test's"}
    if base_file is not None:  # row R51's H files: written only then, so every other file's record is as before
        record["rule"] = (record["rule"]
                          .replace("after the atlas step and before g_mom",
                                   "after MS's step (its momentum, then the atlas step, on the cells whose group has a "
                                   "direction, never Blood) and before g_mom")
                          .replace("as the banked file has it", f"as {base_file} has it"))
    if row is not None:  # row R73's H files (row and base_file given together): the row and the base file named
        record["rule"] = (f"{record['rule']}; Blood's step is row {row}'s one further change to {base_file}'s "
                          f"command (row {T1_R51_ROW}'s file, MS's own damping 1.0), as row {T1_R51_ROW}'s H files "
                          "make it to T1_s2rcMSB's")

    def step(pop):
        labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
        blood = labels == T1_HM_LABEL
        if not blood.any():
            sys.exit("--held-momentum: the draw holds no Blood cell, so the step would be idle")
        idx = np.flatnonzero(blood)
        moved, rep = t1_scope.momentum_shift(start, later, pop[idx].copy(), damping, held,
                                             label_key=t1_scope.LABEL_KEY, eps=T1_HM_EPS, min_cells=T1_HM_MIN_CELLS,
                                             keep_depth=True, within={T1_HM_CATEGORY: [T1_HM_LABEL]},
                                             seed=[seed, T1_HM_STREAM], **kw)
        rows = rep.rows
        if [r["type"] for r in rows] != [T1_HM_CATEGORY] or rows[0]["status"] != "shifted":
            sys.exit("--held-momentum: Blood's step did not shift the Blood cells (values not shown)")
        record["status"] = rows[0]["status"]
        parts = [(idx, moved)]
        if not blood.all():
            rest = np.flatnonzero(~blood)
            parts.append((rest, pop[rest].copy()))
        return t1_scope._stack_rows(parts, pop.n_obs)

    return step, record


def _feng_group_names(feng_groups: str = FENG_GROUPS_DEFAULT) -> list[str]:
    """The names of Feng's direction groups from Feng's statistics alone, before any T1 file is read. The groups
    feng_dirs.feng_directions returns do not depend on the gene list it is given, so it is called on one placeholder
    name (an empty list fails in dirlib.sym_index, whose empty index is a float array) and only the names are kept.
    The default makes exactly the call as before; recomputed asks for the recomputed groups' names."""
    if feng_groups == FENG_GROUPS_DEFAULT:
        return sorted(feng_dirs.feng_directions([""])[0])
    return sorted(_feng_directions([""], feng_groups)[0])


def _t1_union(hygiene_on: bool, feng_groups: str = FENG_GROUPS_DEFAULT) -> dict | None:
    """R2's union for --feng-hygiene, after R2's two refusals that need no T1 file, both made before anything else is
    read (R2: no E10.5 file is read before the build): hygiene.load_union refuses while configs/t1_hygiene.yaml pins no
    union_sha256 (until the build commit after N12) or pins another file's; then the Blood guard on the names of Feng's
    direction groups (on the groups --feng-groups names). Off: None, and nothing is read."""
    if not hygiene_on:
        return None
    try:
        union = hygiene.load_union()
    except hygiene.Refused as e:
        sys.exit(f"--feng-hygiene is refused: {e}")
    _no_blood(_feng_group_names() if feng_groups == FENG_GROUPS_DEFAULT else _feng_group_names(feng_groups))
    return union


def _t1_down_prior(on: bool, hold_unbracketed: bool, hygiene_on: bool) -> t1_scope.DownPrior | None:
    """The option of T12a and T12b (experiments/harness/TRIES.md; the captain's rulings of 2026-11-10 (a) and (d)): the
    strong-null down-tail prior after M1's steps. Off (the default): None, and nothing is read. On: refused unless the
    build is M1 with R1 and R2 (T12 departs from it) and unless configs/t1_strong_null.yaml says approved
    (`strong_null.check_approved`), both before anything else is read. The genes are ordered by the reference's level
    alone (`DownPriorTable.level_only`, the selection the backtests measured); the prior's cells are the ones M1 moves
    along Feng's direction (`covered_labels` None); `t1_scope.apply_down_prior` holds R2's and the pre-registered genes
    exactly as M1's own step does."""
    if not on:
        return None
    if not (hold_unbracketed and hygiene_on):
        sys.exit("--strong-null is refused without R1 and R2 (--no-hold-unbracketed, --no-feng-hygiene): T12 departs "
                 "from M1 with R1 and R2 (experiments/harness/TRIES.md T12a, T12b)")
    try:
        strong_null.check_approved()
    except strong_null.Refused as e:
        sys.exit(f"--strong-null is refused: {e}")
    return t1_scope.DownPrior(table=t1_scope.DownPriorTable.level_only())


def _t1_held(union: dict | None, direction_groups, build_genes) -> tuple[list[str], dict | None]:
    """The held-flat genes build_t1 passes to build_submission, and the build record's feng_hygiene field. Off (no
    union): the pre-registered list exactly as before R2, so the file is byte-identical, and None. On (R2): that list
    joined with the union (_t1_union), after the Blood guard again on the groups the build uses; the field records the
    union's file, sha256 and gene count, the genes per module, the commit and the sources' sha256 (docs/PLAN.md 4.1,
    item 4), and, by name, the union genes absent from the build's genes (N12 mapped the union on E9.5's names; a name
    E10.5 lacks cannot be held and changes nothing)."""
    held = t1_scope.held_flat_genes()
    if union is None:
        return held, None
    _no_blood(direction_groups)
    present = {str(g) for g in build_genes}
    absent = [g for g in union["genes"] if g not in present]
    commit, _ = provenance.git_state(provenance.REPO_ROOT)
    return sorted(set(held) | set(union["genes"])), {
        "file": union["file"], "sha256": union["sha256"], "commit": commit, "genes": len(union["genes"]),
        "union_genes_absent_from_build": {"count": len(absent), "genes": absent}, "modules": union["modules"],
        "sources": union["sources"]}


def _t1_g_mom(ref, categories: list[str], seed: int, covered_labels: set[str], feng_direction: bool = True,
              with_direction=None, damping: float | None = None, fraction_mult: float | None = None,
              line: str | None = None):
    """The test build's g_mom step (--g-mom, default off; T1 try T11 of experiments/harness/TRIES.md; only the captain's
    ruling puts it in an upload): E9.5 is loaded here and nowhere else, and the factors are fitted exactly as the
    candidate's (`vec.harness.g_mom_pair_factors`: start E9.5, reference E10.5, the chamber categories of E10.5's own
    mix, metacells of the E10.5 cells, damped by B0_DAMPING, zeroed and clipped). Returns the `post_cells` hook of
    `t1_scope.build_submission` and a record (categories only, no cell number) that the hook completes when it runs.
    The hook is the build's last count edit, after R1's hold and after Feng's direction (the order of the backtests:
    g_mom was judged on top of B0's per-type shift, and at test Feng's direction takes that place). It edits the cells
    of `vec.t1_state_transport.apply_state_factors`: a cell's source cell must sit in a fitted metacell; cells of R1's
    held labels (no Feng group: `r1_unbracketed_labels`, the candidate's own mask) keep the counts they have, so R1
    stands; R3b: the fit is the candidate's, over all `categories`, and only cells of G_MOM_EDIT_CATEGORIES are edited,
    the mask added to R1's (a category the fit skipped stops the build); the pre-registered held-flat genes are the
    candidate's own list (R2's union is held in Feng's direction only, by its ruling; if the captain wants R2 held
    in this step too, that is a new variant with its own row, since
    T11 was judged without it); each edited cell keeps its expected depth (`apply_factors`, keep_depth), the T1 depth
    lever stays closed; the cells, their number, shares and scope are untouched. The stream is [seed, 5], the
    candidate's. covered_labels: the labels the build moved along Feng's direction (`label_groups` mapped to a group
    with a direction); the hook checks that every label it edits is one of them, so R1's mask and the build's own mask
    of Feng-covered labels cannot disagree silently. feng_direction False (--feng-groups none): the same fit, stream
    and edit; the held labels are `_labels_without_group`'s, the readout config's own rule, so R1_FENG_GROUPS is not
    read (the same cells on every label configs/t1_chamber_classes.yaml places in those categories); the record
    says so (applied last, with no Feng direction in the file), and covered_labels is `_g_mom_covered`'s, so the
    check is a self-check. with_direction (--feng-groups recomputed: the groups with a recomputed direction after
    feng_relabel.yaml's fallback): the same fit, stream and edit; the held labels are `_labels_without_direction`'s,
    the labels whose group has no recomputed direction, so g_mom's mask and R1's hold in the build follow the same
    groups and R1_FENG_GROUPS (fixed under the authors' labels) is not read; the record says so. damping (candidate A,
    --g-mom-damping; None: today's factors and record): the same raw fit finished by
    `vec.t1_state_transport.finish_momentum` at that damping, as practice arm gmom_s1 (experiments/t1_m1/bets_read.py
    gmom_fits) finishes it; at B0_DAMPING that is the call g_mom_pair_factors makes itself. fraction_mult (the
    candidate's direction multiplier; None: not given) enters the record's text only, never the step: with damping
    1.5 or 2.0 it chooses the second line's text (multiplier 1.0), the fifth line's (multiplier 2.0) or the sixth
    line's (multiplier 2.5); at damping 1.0 the text is candidate A's whatever the multiplier (the third and seventh
    lines keep it). line (row R51's files only; None: as above): at a damping other than B0_DAMPING the text names that
    line (_r51_line) in place of every other's, so T1_s2rcMSBd10 and T1_s2rcMSBd15 name row R51, never candidate A or
    row R22's second line; the step is the same."""
    from vec import t1_state_transport as ts

    start = load_stage(DATA / "T1/E9.5_RNA.h5ad", lineage=False)
    factors, raw = harness.g_mom_pair_factors(ref, "E10.5_RNA", start, "E9.5_RNA", categories=categories)
    if damping is not None and damping != harness.B0_DAMPING:  # candidate A: the same fit, finished at `damping`
        factors = ts.finish_momentum(raw, damping)
    missing = sorted(set(G_MOM_EDIT_CATEGORIES) - set(factors.categories))
    if missing:
        raise ValueError(f"--g-mom has no factors for a category R3b edits: {missing}")
    held_flat = t1_scope.held_flat_genes()
    record = {"rule": "try T11's within-type factors at test: OT within each chamber category from E9.5 to E10.5, "
                      "metacells of the E10.5 cells, the per-category mean removed, damped by B0_DAMPING, zeroed and "
                      "clipped; applied last (after Feng's direction); R1's held labels unedited; keep_depth; R3b: "
                      "only the categories in edit_categories are edited; the factors are not scaled to the two-day "
                      "step",
              "edit_categories": list(G_MOM_EDIT_CATEGORIES),
              "categories_with_factors": sorted(factors.categories),
              "categories_without_factors": sorted(factors.skipped), "stream": [int(seed), 5],
              "held_flat": "the candidate's pre-registered list (R2's union is held in Feng's direction only)"}
    if not feng_direction:  # --feng-groups none: written only then, so every other file's record is as before
        record["rule"] = record["rule"].replace("applied last (after Feng's direction)",
                                                "applied last (no Feng direction in this file)")
        record["rule"] = record["rule"].replace("R1's held labels unedited",
                                                "labels the readout config's regexes give no group unedited")
        record["held_flat"] = "the candidate's pre-registered list"
        record["mask"] = ("the cells of edit_categories whose label the readout config's regexes give a group (the "
                          "readout config's own rule; the list of groups with a Feng direction is not read)")
    if with_direction is not None:  # --feng-groups recomputed: written only then
        record["rule"] = record["rule"].replace("applied last (after Feng's direction)",
                                                "applied last (after Feng's direction on the recomputed groups)")
        record["rule"] = record["rule"].replace("R1's held labels unedited",
                                                "labels of groups without a recomputed direction unedited")
        record["mask"] = ("the cells of edit_categories whose label's group has a recomputed direction after the "
                          f"fallback ({sorted(with_direction)}); R1_FENG_GROUPS is not read")
    if damping is not None:  # written only with --g-mom-damping, so every other file's record is as before
        record["damping"] = float(damping)
        if line is not None and damping != harness.B0_DAMPING:  # row R51's two g_mom files (T1's line 10)
            record["rule"] = record["rule"].replace("damped by B0_DAMPING",
                                                    f"damped by {float(damping):g} in place of B0_DAMPING ({line})")
        elif damping in T1_R22_LINE5_OPTIONS and fraction_mult == FENG_FRACTION_MULT_LINE5:  # row R22's fifth line
            record["rule"] = record["rule"].replace("damped by B0_DAMPING",
                                                    f"damped by {float(damping):g} in place of B0_DAMPING (row "
                                                    f"{T1_R22_ROW}'s fifth T1 search line, on the 2.0 file's command, "
                                                    "docs/PLAN.md decision log, T1 line 5)")
        elif damping in T1_R22_LINE6_OPTIONS and fraction_mult == FENG_FRACTION_MULT_LINE6:  # row R22's sixth line
            record["rule"] = record["rule"].replace("damped by B0_DAMPING",
                                                    f"damped by {float(damping):g} in place of B0_DAMPING (row "
                                                    f"{T1_R22_ROW}'s sixth T1 search line, on the 2.5 file's command, "
                                                    "docs/PLAN.md decision log 2026-11-24 (bb), T1 line 6)")
        elif damping in T1_R22_DAMPING_OPTIONS:  # row R22's second T1 search line; candidate A's text below as before
            record["rule"] = record["rule"].replace("damped by B0_DAMPING",
                                                    f"damped by {float(damping):g} in place of B0_DAMPING (row "
                                                    f"{T1_R22_ROW}'s second T1 search line, docs/PLAN.md decision log "
                                                    "2026-11-22)")
        elif damping != harness.B0_DAMPING:
            record["rule"] = record["rule"].replace("damped by B0_DAMPING",
                                                    f"damped by {float(damping):g} in place of B0_DAMPING (candidate "
                                                    "A, docs/PLAN.md decision log 2026-11-20 (i))")

    def hook(pop):
        labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
        cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
        if with_direction is not None:
            held_labels = _labels_without_direction(np.unique(labels), with_direction)
        else:
            held_labels = (ts.r1_unbracketed_labels(np.unique(labels)) if feng_direction
                           else _labels_without_group(np.unique(labels)))
        no_edit = np.isin(labels, sorted(held_labels)) | ~np.isin(cats, G_MOM_EDIT_CATEGORIES)
        edited = set(labels[~no_edit & np.isin(cats, list(factors.categories))])
        if not edited <= set(covered_labels):
            extra = sorted(edited - set(covered_labels))
            what = "Feng's direction does not cover" if feng_direction else "the readout config gives no group"
            raise ValueError(f"--g-mom would edit labels {what}: {extra}")
        moved, rec = ts.apply_state_factors(pop, factors, held_flat, [seed, 5], rescale=False, no_edit=no_edit)
        record["median_depth_scale"] = rec["median_depth_scale"]  # no cell number
        return moved

    return hook, record


def _t1_mid_rows(ref, in3: np.ndarray, seed: int) -> tuple[np.ndarray, dict]:
    """The mid base's pool (row R7; try T33's arm as a builder gets it): every library-3 cell of the reference, and in
    each chamber category a random fraction s / (1 + s) of the other libraries' cells, s being library 3's share of that
    category's cells in the file. Categories in sorted order, one generator, default_rng([seed, T1_MID_STREAM]). The
    other libraries' cells outside the chamber categories are left out (no draw takes them). Returns the kept rows in
    file order and a record of shares only (no cell number)."""
    spec = t1_scope.load_spec()
    b = harness._shallow(ref)
    t1_scope.classify(b, source="E10.5_RNA", spec=spec)
    cat = np.asarray(b.obs[t1_scope.CATEGORY_KEY]).astype(str)
    rng = np.random.default_rng([seed, T1_MID_STREAM])
    keep = [np.flatnonzero(in3)]
    record = {}
    for c in sorted(spec.in_scope("chamber")):
        rows3, rows_o = np.flatnonzero(in3 & (cat == c)), np.flatnonzero(~in3 & (cat == c))
        if len(rows3) + len(rows_o) == 0:
            continue
        if len(rows3) == 0:
            sys.exit(f"--base mid: library 3 holds no cell of the chamber category {c}, so the mid base is undefined")
        s = len(rows3) / (len(rows3) + len(rows_o))
        k = int(round(len(rows_o) * s / (1.0 + s)))
        kept = np.sort(rng.choice(rows_o, size=k, replace=False)) if k else np.array([], dtype=np.int64)
        keep.append(kept)
        record[c] = {"library3_share_in_file": s, "other_fraction_kept": k / len(rows_o) if len(rows_o) else 0.0,
                     "library3_share_in_pool": len(rows3) / (len(rows3) + k)}
    return np.sort(np.concatenate(keep)), record


def _t1_trend(ref, portions: np.ndarray, shares: dict, lam: float) -> tuple[dict, dict]:
    """--released-bet share_trend: the pool's own chamber shares moved along the like-for-like trend from E9.5 to
    E10.5's portion libraries (`portions`, a mask of the reference's rows), `t1_bets.trend_shares` at `lam`. E9.5 is
    loaded here for its labels' counts. Returns the shares and a record (shares and ratios, no count)."""
    start = load_stage(DATA / "T1/E9.5_RNA.h5ad", lineage=False)
    later = t1_bets.chamber_counts(ref[np.flatnonzero(portions)], "E10.5_RNA")
    new, rec = t1_bets.trend_shares(shares, later, t1_bets.chamber_counts(start, "E9.5_RNA"), lam)
    return t1_scope.load_shares(new), {"rule": "library 3's own chamber mix x (E10.5 portion libraries' share / E9.5's "
                                               "share) ** lam, renormalised; ratio 1 below the trend's min_cells",
                                       **rec}


def _t1_advance(ref, pool, portions: np.ndarray, shares: dict, n: int, seed: int, held: list[str], board, hook,
                tau: float) -> tuple:
    """--released-bet advance: the pool's cells drawn at `shares`, each category's quota split over z-quantile bins in
    proportion to exp(tau z) (`t1_bets.advance_reweight`), z on each category's momentum axis from E9.5 to E10.5's
    portion libraries (`portions`), the held genes at 0; then g_mom's hook and the board's genes
    (`t1_bets.assemble`), as `t1_scope.build_submission` does on the no-shift path (at tau 0 the same file, tested).
    Returns the submission, a report with the realised composition (stripped from the record on a test board), a
    record and the drawn population (before the hook), for its median library size."""
    if n < board.min_cells or (board.max_cells is not None and n > board.max_cells):
        raise ValueError("n_cells is outside the board's limits")
    spec = t1_scope.load_spec()
    start = harness._shallow(load_stage(DATA / "T1/E9.5_RNA.h5ad", lineage=False))
    later = ref[np.flatnonzero(portions)].copy()
    t1_scope.classify(start, source="E9.5_RNA", spec=spec)
    t1_scope.classify(later, source="E10.5_RNA", spec=spec)
    t1_scope.classify(pool, source="E10.5_RNA", spec=spec)  # in place, as build_submission classifies its ref
    axes, axes_rec = t1_bets.momentum_axes(start, later, [c for c, s in shares.items() if s > 0], held,
                                           start_fallback=spec.fallbacks("E9.5_RNA"),
                                           later_fallback=spec.fallbacks("E10.5_RNA"))
    z = t1_bets.advance_scores(pool, axes)
    pop, _ = t1_bets.advance_reweight(pool, shares, n, seed, z, tau, fallback=spec.fallbacks("E10.5_RNA"), spec=spec)
    sub = t1_bets.assemble(pop, board, hook)
    return sub, {"realised": t1_scope.composition(pop)}, {
        "rule": "library 3's cells drawn inside each chamber category by quotas over five z-quantile bins in "
                "proportion to each bin's sum of exp(tau z), z the cell's projection on its category's E9.5 -> "
                "E10.5-portion-libraries momentum axis, standardised within the category; category quotas at library "
                "3's own chamber mix; then g_mom", "tau": float(tau), "axes": axes_rec}, pop


def _bet_median_library(pool, base_shares: dict, bet_pop, n: int, seed: int) -> dict:
    """docs/DATA_POLICY.md 2026-11-04 (f)(3) for --released-bet: the bet's drawn cells' median library size over the
    none file's draw at the same seed (`reweight` on the same pool at the pool's own chamber mix, as build_submission
    draws it on that file's path). A ratio of two medians of released cells, no count."""
    spec = t1_scope.load_spec()
    t1_scope.classify(pool, source="E10.5_RNA", spec=spec)  # in place and idempotent, as build_submission does
    base, _ = t1_scope.reweight(pool, base_shares, n, seed, fallback=spec.fallbacks("E10.5_RNA"), spec=spec)
    return {"median_library_over_none_draw": t1_bets.median_library(bet_pop) / t1_bets.median_library(base),
            "against": "the none file's draw at the same seed (library 3's cells at their own chamber mix)"}


def _late_mix_median_library(pool, own_shares: dict, drawn: list, n: int, seed: int) -> dict:
    """docs/DATA_POLICY.md 2026-11-04 (f)(3) for row R54's file: the median library of the cells its build drew (taken
    in its post_cells, once) over that of T1_s2rcMSB's draw at the same seed (`reweight` on the same pool at the pool's
    own chamber mix, the call build_submission makes on that file's path). A ratio of two medians of released cells,
    no count; a shift beyond 5% goes to the captain by (f)(3), which row R54's named exception covers."""
    if len(drawn) != 1:
        raise ValueError("row R54's build drew its cells more or less than once")
    spec = t1_scope.load_spec()
    t1_scope.classify(pool, source="E10.5_RNA", spec=spec)  # in place and idempotent, as build_submission does
    base, _ = t1_scope.reweight(pool, own_shares, n, seed, fallback=spec.fallbacks("E10.5_RNA"), spec=spec)
    return {"median_library_over_msb_draw": float(drawn[0]) / t1_bets.median_library(base),
            "median_library_against": "T1_s2rcMSB's draw at the same seed (library 3's cells at their own chamber mix)"}


def build_t1(a) -> dict:
    feng_groups = getattr(a, "feng_groups", FENG_GROUPS_DEFAULT)
    # None unless --late-mix is given (row R54's file), tested before every route below, row R51's included, none of
    # which then sees it; then None unless the command names --captain-ruling R58 (row R58's file), tested before row
    # R51's route. Their refusals read nothing but the ledger (and, for R54, its committed shares file)
    r54 = _t1_r54(a)
    r58 = None if r54 is not None else _t1_r58(a)
    # then None unless the command names --captain-ruling R67 (row R67's three files, row R58's route over three
    # literals), tested after row R58's and before row R51's; its refusals read nothing but the ledger
    r67 = None if r54 is not None or r58 is not None else _t1_r67(a)
    # None unless the command is one of row R51's seven files (T1's line 10), tested before every route below, none of
    # which then sees it: its values stand in for theirs; its refusals read nothing but the ledger
    r51 = _t1_r51(a) if r54 is None and r58 is None and r67 is None else None
    ruled = r58 if r58 is not None else r67  # row R58's or row R67's values: one literal damping on MSB's command
    line = r51 if ruled is None else ruled  # row R51's, R58's or R67's values, which stand in for the routes below
    if r54 is not None:  # T1_s2rcMSB's steps exactly: MS's step at B0's damping, the multiplier 1.5, g_mom at 0.5
        hm, ms, candidate = None, r54, r54
    elif line is None:
        hm = _t1_hm(a)  # None unless --held-momentum is given (Blood's momentum); its refusals read only the ledger
        ms = _t1_ms(a)  # None unless --direction-momentum is given (candidate MS, row R23); its refusals read nothing
        candidate = _t1_candidate(a)  # None unless --g-mom-damping or --feng-fraction-mult is given; reads nothing
    else:  # MS's step always; Blood's in the two H files; the multiplier and g_mom's damping from the row's values
        hm, ms, candidate = (None if line["held_momentum"] is None else line), line, line
    pooled = _t1_pooled(a)  # None unless --base pooled --feng-groups recomputed --captain-ruling R37; likewise
    _bet_refusals(a)  # nothing at all unless --released-bet (or one of its two companions) is given
    if feng_groups == "none":
        _feng_none_refusals(a)
    elif feng_groups == "recomputed":
        _feng_recomputed_refusals(a)
    elif getattr(a, "feng_none_ruling", None) is not None:
        sys.exit("--feng-none-ruling is refused without --feng-groups none")
    if a.g_mom and a.strong_null:
        sys.exit("--g-mom is refused with --strong-null: the two count edits have no ruled order and were never tested "
                 "together; a new row in experiments/harness/TRIES.md and the captain's ruling fix it first")
    if a.g_mom and not a.hold_unbracketed:
        sys.exit("--g-mom is refused with --no-hold-unbracketed: g_mom's factors at test add to M1 with R1 and R2 as "
                 "they are (experiments/harness/TRIES.md, try T11)")
    if a.lib3_pair and (a.base != "lib3" or a.cm_share is None or a.blood_share is not None):
        sys.exit("--lib3-pair names the second file of the pair try T32 priced (library 3's cells with the "
                 "cardiomyocyte share moved): it needs --base lib3 with --cm-share, and no --blood-share")
    if a.base == "lib3" and (a.cm_share is not None or a.blood_share is not None) and not a.lib3_pair:
        sys.exit("--base lib3 is refused with --cm-share or --blood-share: the pair priced in try T30 keeps the "
                 "second upload on the pooled base (experiments/harness/TRIES.md, row R6). --lib3-pair builds try "
                 "T32's second file instead; it is no entry, and it is uploaded only by the captain's own ruling")
    if a.base == "lib3" and (not a.hold_unbracketed or a.strong_null or a.herd):
        sys.exit("--base lib3 is refused with --no-hold-unbracketed, --strong-null or --herd: each would fit or pull a "
                 "step on library 3's cells alone, which no try has read (experiments/harness/TRIES.md, row R6)")
    if a.base == "mid" and (a.cm_share != T1_MID_CM_SHARE or a.blood_share is not None or a.lib3_pair
                            or not a.hold_unbracketed or a.strong_null or a.herd):
        sys.exit("--base mid builds one file, T1's second upload by the captain's ruling of 2026-11-13 (row R7 of "
                 f"experiments/harness/TRIES.md; try T33): it needs --cm-share {T1_MID_CM_SHARE:.2f}, and it is "
                 "refused with another share, --blood-share, --lib3-pair, --no-hold-unbracketed, --strong-null or "
                 "--herd, which no try has read on that base")
    # R2's refusals come before anything is read. Under none R2 is off: it holds genes flat in Feng's direction, and
    # there is none (an explicit --feng-hygiene, the one case _feng_hygiene_on could refuse, was refused above)
    hygiene_on = False if feng_groups == "none" else _feng_hygiene_on(a.feng_hygiene, a.hold_unbracketed)
    prior = _t1_down_prior(a.strong_null, a.hold_unbracketed, hygiene_on)  # T12's, likewise
    union = _t1_union(hygiene_on) if feng_groups == FENG_GROUPS_DEFAULT else _t1_union(hygiene_on, feng_groups)
    board, dry = _board("t1", a.board)
    ref = load_stage(DATA / "T1/E10.5_RNA.h5ad", lineage=False)
    momentum = _t1_momentum(a.hold_unbracketed)
    shares = harness.t1_chamber_shares(harness._shallow(ref), "E10.5_RNA")
    fit_categories = [c for c, s in shares.items() if s > 0]  # g_mom's fit, as the candidate's: E10.5's own mix
    pool = ref
    mid_record = None
    if a.base in ("lib3", "mid"):
        names = [str(x) for x in ref.obs_names]
        library = np.array([x.rsplit("_", 1)[-1] if "_" in x else "" for x in names])
        kinds = set(library.tolist())
        if not all(k.isdigit() for k in kinds) or T1_DISSECTED_LIBRARY not in kinds or len(kinds) < 2:
            sys.exit(f"--base {a.base}: the reference's cell names do not separate library 3 from the others (every "
                     "name must end in _<library number>, and library 3 must be one of at least two)")
        in3 = library == T1_DISSECTED_LIBRARY
        if a.base == "lib3":
            pool = ref[np.flatnonzero(in3)].copy()
        else:
            rows, mid_record = _t1_mid_rows(ref, in3, a.seed)
            pool = ref[rows].copy()
        shares = harness.t1_chamber_shares(harness._shallow(pool), "E10.5_RNA")  # the pool's own chamber mix
    late_own = late_record = None
    if r54 is not None:  # row R54's file: the pool's own mix moved to row N64's shares (its route refused --cm-share
        # and --blood-share, so row R6's refusal of those on --base lib3 is not met and stays as it is)
        late_own = shares
        shares, late_record = _t1_late_mix_shares(late_own, r54)
    if a.cm_share is not None:
        shares = t1_scope.with_share(shares, "cm", a.cm_share)
    if a.blood_share is not None:  # blood's share of all chamber cells, the other categories in proportion
        shares = t1_scope.with_share(shares, "blood", a.blood_share, within=list(shares))
    directions, ext = _feng_directions(list(ref.var_names), feng_groups)
    held, feng_hygiene = _t1_held(union, directions, ref.var_names)
    groups = feng_dirs.label_groups(ref.obs["celltype"].astype(str).unique())
    calendar = feng_dirs.bracket_fraction(10.5, 12.5)
    # candidate B (--feng-fraction-mult): the multiplier times the calendar fraction, as double_step_read.py forms its
    # doubled arm (MULTIPLIER x fraction); without the option the calendar fraction itself, today's value
    fraction = calendar if candidate is None else candidate["feng_fraction_mult"] * calendar
    n = harness.count_for(board, ref.n_obs)
    extra = {} if prior is None else {"down_prior": prior}
    within_type, g_mom_record = ({}, None)
    if a.g_mom:
        covered = _g_mom_covered(groups, directions, feng_groups)  # the build's Feng-covered labels (none: self-check)
        if feng_groups == "recomputed":  # g_mom's mask follows the groups with a recomputed direction
            hook, g_mom_record = _t1_g_mom(ref, fit_categories, a.seed, covered, with_direction=sorted(directions),
                                           damping=None if candidate is None else candidate["g_mom_damping"],
                                           fraction_mult=None if candidate is None else candidate["feng_fraction_mult"],
                                           **({} if r51 is None else {"line": _l12_line(r51)}))
        else:
            hook, g_mom_record = _t1_g_mom(ref, fit_categories, a.seed, covered,
                                           feng_direction=feng_groups != "none")
        within_type = {"post_cells": hook}
    bet = getattr(a, "released_bet", None)
    bet_record = bet_pop = None
    own_shares = shares
    if bet == "share_trend":  # refused above unless --feng-groups none, so the pool is library 3's (row R8's command)
        shares, bet_record = _t1_trend(ref, library != T1_DISSECTED_LIBRARY, shares, a.bet_strength)
    hm_record = None
    if hm is not None:  # Blood's momentum (refused above unless the banked file's command, so --g-mom and library 3):
        # the banked file's draw and atlas step in build_submission, then Blood's step and g_mom's hook in post_cells.
        # Row R51's H files: the same step after MS's step below (it wraps this hook), on T1_s2rcMSB's command
        hm_step, hm_record = _t1_hm_step(ref, library, hm["held_momentum"], held, a.seed,
                                         **({} if r51 is None else {"base_file": r51["base_file"]}),
                                         **({"row": T1_R73_ROW} if r51 is not None
                                            and r51["captain_ruling"] == T1_R73_ROW else {}))
        hm_hook = within_type["post_cells"]
        within_type = {"post_cells": lambda pop: hm_hook(hm_step(pop))}
    ms_record = None
    if ms is not None:  # candidate MS (refused above unless the recomputed command with --g-mom): the same draw, then
        # the momentum and the atlas step on the covered cells, then g_mom's hook, all inside post_cells. Under row
        # R26 (the pair) the same calls, with `fraction` and g_mom's damping candidate B's from _t1_candidate above;
        # under row R51, R58 or R67 the same calls with its values, and MS's own damping passed to the step (row R58's
        # and each row R67 file's with its own record text); under row R54 R26's calls exactly, on the moved shares
        note = T1_R58_NOTE if r58 is not None else None if r67 is None else T1_R67_NOTES[r67["file"]]
        ms_step, ms_record = _t1_ms_step(ref, library, directions, fraction, groups, held, a.seed,
                                         **({} if line is None else {"damping": line["direction_momentum_damping"]}),
                                         **({} if note is None else {"option_note": note}))
        g_mom_hook = within_type["post_cells"]
        late_drawn = []

        def ms_post(pop):
            if r54 is not None:  # row R54's file: the drawn cells' median library, for its (f)(3) record
                late_drawn.append(t1_bets.median_library(pop))
            return g_mom_hook(ms_step(pop))

        sub, rep = t1_scope.build_submission(pool, board, shares, n, a.seed, held_flat=held, ref_source="E10.5_RNA",
                                             herd=a.herd, post_cells=ms_post)
        if r54 is not None:  # docs/DATA_POLICY.md 2026-11-04 (f)(3): the draw's median library over T1_s2rcMSB's draw
            late_record = {**late_record, **_late_mix_median_library(pool, late_own, late_drawn, n, a.seed)}
    elif bet == "advance":
        sub, rep, bet_record, bet_pop = _t1_advance(ref, pool, library != T1_DISSECTED_LIBRARY, shares, n, a.seed,
                                                    held, board, within_type["post_cells"], a.bet_strength)
    else:
        sub, rep = t1_scope.build_submission(pool, board, shares, n, a.seed, held_flat=held,
                                             ref_source="E10.5_RNA", directions=directions,
                                             direction_fraction=fraction, direction_groups=groups, herd=a.herd,
                                             **momentum, **extra, **within_type)
    if bet == "share_trend":  # the draw build_submission made at the trended shares (the same call), for (f)(3)
        spec = t1_scope.load_spec()
        bet_pop = t1_scope.reweight(pool, shares, n, a.seed, fallback=spec.fallbacks("E10.5_RNA"), spec=spec)[0]
    if bet_record is not None:
        bet_record = {**bet_record, **_bet_median_library(pool, own_shares, bet_pop, n, a.seed)}
    strong = None
    if prior is not None:  # T12's part (d): the file is written only if the target-free checks pass
        strong = {"selection": "the reference's level alone (vec.t1_scope.DownPriorTable.level_only), no table",
                  "constants": {"budget": prior.budget, "band": prior.band, "f_min": prior.f_min, "tol": prior.tol,
                                "mmd_bound": prior.mmd_bound},
                  **rep["down_prior"]}
        if not strong["checks"]["passed"]:
            sys.exit(f"the strong-null prior's target-free checks failed: {strong['checks']}")
    check = t1_scope.write_checked(sub, a.out, board)
    rec = {"entry": "M1", "board": board.key, "dry_run": dry, "n_cells": n, "seed": a.seed, "fraction": fraction,
           "herd": bool(a.herd), "hold_unbracketed": bool(a.hold_unbracketed), "feng_hygiene": feng_hygiene,
           "base": a.base,
           "cm_share_of_heart": t1_scope.share_of(shares, "cm", t1_scope.load_spec().heart),
           "blood_share_of_chamber": float(shares.get("blood", 0.0)),
           "shares": shares, "external": ext, "check_ok": check.ok, "realised": rep["realised"]}
    if a.lib3_pair:  # written only for that file, so every other file's record is as before
        rec["lib3_pair"] = True
    if feng_groups == "none":  # likewise: no direction step ran, so it has no fraction; no Feng file was read
        rec["fraction"] = None
        rec["feng_groups"] = {"groups": feng_groups, "file": ext["file"], "sha256": ext["sha256"],
                              "direction_step": "removed: no cell is moved along any direction; \"entry\": \"M1\" "
                                                "names the design this file comes from"}
        if getattr(a, "feng_none_ruling", None) is not None:  # only after a relabelling routed to recomputed
            rec["feng_groups"]["ruling_over_the_relabelling_route"] = a.feng_none_ruling
    if feng_groups == "recomputed":  # likewise: written only then; "external" names the same file
        rec["feng_groups"] = {"groups": feng_groups, "file": ext["file"], "sha256": ext["sha256"],
                              "relabel_config_sha256": ext["relabel_config_sha256"], "route": ext["route"],
                              "with_direction": ext["with_direction"], "fallbacks": ext["fallbacks"],
                              "relabel_gate": ext["relabel_gate"],
                              "groups_without_direction": "no cell of theirs moves: R1 holds them at their E10.5 "
                                                          "counts and g_mom's mask holds them back"}
    if candidate is not None:  # written only with --g-mom-damping or --feng-fraction-mult ("fraction" is the applied)
        rec["t1_candidate"] = {**candidate, "calendar_fraction": calendar, "applied_fraction": fraction,
                               "command": T1_CANDIDATE_COMMAND}
    if ms is not None:  # written only with --direction-momentum, so every other file's record is as before
        # row R26's file (the pair; candidate B's values came through _t1_candidate above) names its own command; the
        # banked command's MS file carries its command (with its row) in its values
        rec["t1_candidate"] = {**ms, "momentum": ms_record, "calendar_fraction": calendar, "applied_fraction": fraction,
                               "command": ms["command"] if "command" in ms else
                               T1_MSB_COMMAND if ms["captain_ruling"] == T1_R26_ROW else T1_MS_COMMAND}
    if hm is not None:  # written only with --held-momentum, so every other file's record is as before
        rec["t1_candidate"] = {**hm, "momentum": hm_record, "calendar_fraction": calendar, "applied_fraction": fraction}
    if r51 is not None:  # row R51's files, in place of the fields above: MS's record as in every MS file, Blood's too
        rec["t1_candidate"] = {**r51, "momentum": _l12_momentum(r51, ms_record), "calendar_fraction": calendar,
                               "applied_fraction": fraction}
        if hm_record is not None:
            rec["t1_candidate"]["held_momentum_step"] = hm_record
    if ruled is not None:  # R58's file or an R67 file, in place of the fields above: MS's record as in every MS file
        rec["t1_candidate"] = {**ruled, "momentum": ms_record, "calendar_fraction": calendar,
                               "applied_fraction": fraction}
    if r54 is not None:  # row R54's file, in place of the fields above; its late mix in a field of its own
        rec["t1_candidate"] = {**{k: v for k, v in r54.items() if k not in ("own_shares", "late_shares")},
                               "momentum": ms_record, "calendar_fraction": calendar, "applied_fraction": fraction}
        rec["late_mix"] = late_record
    if pooled is not None:  # written only for row R37's file, so every other file's record is as before
        rec["base_exception"] = pooled
    if mid_record is not None:  # likewise: shares only, no cell number
        rec["mid_base"] = {"stream": [int(a.seed), T1_MID_STREAM], "categories": mid_record}
    if a.g_mom:
        rec["g_mom"] = g_mom_record
    if bet_record is not None:  # written only with --released-bet, so every other file's record is as before
        rec["released_bet"] = {"arm": _bet_arm(bet, a.bet_strength), "captain_ruling": a.captain_ruling,
                               **bet_record}
    return rec if strong is None else {**rec, "strong_null": strong}


def _heart_line_ordered() -> bool:
    """Whether experiments/harness/TRIES.md holds exactly one row R28 whose first cell is the captain's order naming
    `--trend-strength 0.5`, `--trend-strength 1.0` and `--captain-ruling R28` (_captain_order_names)."""
    return all(_captain_order_names(HEART_TS_ROW, opt) for opt in HEART_TS_OPTIONS.values())


def _heart_d_line_ordered() -> bool:
    """Whether experiments/harness/TRIES.md holds exactly one row R29 whose first cell is the captain's order ("ruling,
    no run") naming `--entry organ_scope`, `--trend-strength 0.5`, `--trend-strength 1.0` and `--captain-ruling R29`
    (_captain_order_names)."""
    return all(_captain_order_names(HEART_OTS_ROW, opt) for opt in (HEART_OTS_ENTRY, *HEART_TS_OPTIONS.values()))


def _heart_further_ordered(row: str) -> bool:
    """Whether `row` is a ruling row's id (HEART_FURTHER_ROW_FORM) and experiments/harness/TRIES.md holds exactly one
    row `row` whose first cell is the captain's order ("ruling, no run") naming `--entry organ_scope`,
    `--trend-strength 0.75`, `--trend-strength 1.25` and `--captain-ruling <row>` (_captain_order_names)."""
    return (HEART_FURTHER_ROW_FORM.fullmatch(row) is not None
            and all(_captain_order_names(row, opt) for opt in (HEART_OTS_ENTRY, *HEART_FURTHER_OPTIONS.values())))


def _heart_ots15_ordered() -> bool:
    """Whether experiments/harness/TRIES.md holds exactly one row R52 whose first cell is the captain's order ("ruling,
    no run") naming `--entry organ_scope`, `--trend-strength 1.5` and `--captain-ruling R52` (_captain_order_names)."""
    return all(_captain_order_names(HEART_OTS15_ROW, opt) for opt in (HEART_OTS_ENTRY, HEART_OTS15_OPTION))


def _ruling_first_cell(row: str) -> str:
    """The first cell of TRIES.md's one row `row`, as _captain_order_names reads it (for the build record; read once
    that check has passed)."""
    lines = [ln for ln in TRIES_FILE.read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {row} ")]
    if len(lines) != 1:
        sys.exit(f"experiments/harness/TRIES.md no longer holds exactly one row {row}")
    return lines[0][2:].split("|", 1)[0].strip()


def _further_refusals(a, raw: str, ruling: str | None) -> float:
    """--trend-strength 0.75 or 1.25 (file D's two further strengths): its refusals, before anything is read. Refused
    except as heart --entry organ_scope --trend-strength <value> --z-mode raw --captain-ruling <ROW> on the test board,
    with <ROW> a ruling row's id, and unless TRIES.md holds exactly one row <ROW> whose first cell is the captain's
    order naming `--entry organ_scope`, the two options and `--captain-ruling <ROW>`; then, since row R52's change,
    unless <ROW> is R49 (HEART_FURTHER_ROW), the pin coming last so that every refusal before it is word for word as it
    was. Returns the strength."""
    command = HEART_FURTHER_COMMANDS[float(raw)].format(row="<ROW>")
    if (a.entry != "organ_scope" or a.z_mode != "raw" or a.expression != "b0" or a.trend or a.depth or a.herd
            or a.cfm):
        sys.exit(f"--trend-strength {raw} is refused except as {command}: --entry organ_scope and --z-mode raw only, "
                 "and with --expression copy, --trend, --depth, --herd and --cfm")
    if a.board not in (None, TEST_BOARDS["heart"]):
        sys.exit(f"--trend-strength {raw} is refused on {a.board}: the heart search line builds on the test board "
                 "only (no dry run)")
    if ruling is None:
        sys.exit(f"--entry organ_scope --trend-strength {raw} needs --captain-ruling <ROW>, the row of "
                 "experiments/harness/TRIES.md that is the captain's order for file D's two further strengths")
    if HEART_FURTHER_ROW_FORM.fullmatch(ruling) is None:
        sys.exit(f"--entry organ_scope --trend-strength {raw} is refused with --captain-ruling {ruling}: <ROW> is a "
                 "ruling row's id, R and its number (e.g. R44)")
    if not _heart_further_ordered(ruling):
        sys.exit(f"--entry organ_scope --trend-strength {raw} is refused: experiments/harness/TRIES.md must hold "
                 f"exactly one row {ruling} whose first cell is the captain's order (\"ruling, no run\") naming "
                 f"{HEART_OTS_ENTRY}, {', '.join(HEART_FURTHER_OPTIONS.values())} and `--captain-ruling {ruling}`")
    if ruling != HEART_FURTHER_ROW:
        sys.exit(f"--entry organ_scope --trend-strength {raw} is refused with --captain-ruling {ruling}: only row "
                 f"{HEART_FURTHER_ROW} orders it")
    return float(raw)


def _ots15_refusals(a, raw: str, ruling: str | None) -> float:
    """--trend-strength 1.5 (file D's strength 1.5, row R52): its refusals, before anything is read. Refused except as
    heart --entry organ_scope --trend-strength 1.5 --z-mode raw --captain-ruling R52 on the test board (no dry run),
    and unless TRIES.md holds exactly one row R52 whose first cell is the captain's order ("ruling, no run") naming
    `--entry organ_scope`, `--trend-strength 1.5` and `--captain-ruling R52` (_heart_ots15_ordered). Returns 1.5."""
    if (a.entry != "organ_scope" or a.z_mode != "raw" or a.expression != "b0" or a.trend or a.depth or a.herd
            or a.cfm):
        sys.exit(f"--trend-strength {raw} is refused except as {HEART_OTS15_COMMAND}: --entry organ_scope and --z-mode "
                 "raw only, and with --expression copy, --trend, --depth, --herd and --cfm")
    if a.board not in (None, TEST_BOARDS["heart"]):
        sys.exit(f"--trend-strength {raw} is refused on {a.board}: the heart search line builds on the test board "
                 "only (no dry run)")
    if ruling is None:
        sys.exit(f"--entry organ_scope --trend-strength {raw} is the captain's order of row {HEART_OTS15_ROW} alone "
                 f"and needs --captain-ruling {HEART_OTS15_ROW}")
    if ruling != HEART_OTS15_ROW:
        sys.exit(f"--entry organ_scope --trend-strength {raw} is refused with --captain-ruling {ruling}: only row "
                 f"{HEART_OTS15_ROW} orders it")
    if not _heart_ots15_ordered():
        sys.exit(f"--entry organ_scope --trend-strength {raw} is refused: experiments/harness/TRIES.md must hold "
                 f"exactly one row {HEART_OTS15_ROW} whose first cell is the captain's order (\"ruling, no run\") "
                 f"naming {HEART_OTS_ENTRY}, {HEART_OTS15_OPTION} and `--captain-ruling {HEART_OTS15_ROW}`")
    return HEART_OTS15_STRENGTH


def _trend_strength_refusals(a) -> float | None:
    """--trend-strength (the heart search line, row R28 on the full trim; row R29 on the organ cut; file D's two
    further strengths under row R49, the ruling row given; file D's strength 1.5 under row R52): its refusals, before
    anything is read. Returns the strength, or None when the option is not given (today's commands exactly;
    --captain-ruling alone is refused). Every value but 0.75, 1.25 and 1.5 (0.25 included) meets the refusals it met
    before those were added, word for word, and 0.75 and 1.25 meet theirs of before row R52's change. With
    --trend-form (row R55's trend-form line) _trend_form_refusals decides alone, first."""
    form = getattr(a, "trend_form", None)
    if form is not None:
        if getattr(a, "captain_ruling", None) == HEART_R74_ROW:  # row R74's four files, from its one table
            return _r74_refusals(a, form)
        return _trend_form_refusals(a, form)
    raw = getattr(a, "trend_strength", None)
    ruling = getattr(a, "captain_ruling", None)
    commands = (f"{HEART_TS_COMMANDS[0.5]}, {HEART_TS_COMMANDS[1.0]}, {HEART_OTS_COMMANDS[0.5]} or "
                f"{HEART_OTS_COMMANDS[1.0]}")
    if raw is None:
        if ruling is not None:
            sys.exit(f"heart --captain-ruling goes with --trend-strength only ({commands})")
        return None
    if raw in HEART_FURTHER_CHOICES:  # file D's two further strengths: the ruling row is the one given, R49
        return _further_refusals(a, raw, ruling)
    if raw == HEART_OTS15_CHOICE:  # file D's strength 1.5: row R52's alone
        return _ots15_refusals(a, raw, ruling)
    if raw not in HEART_TS_CHOICES:
        sys.exit(f"--trend-strength {raw} is refused: it takes {' or '.join(HEART_TS_CHOICES)} only")
    s = float(raw)
    if (a.entry not in ("core_trim", "organ_scope") or a.z_mode != "raw" or a.expression != "b0" or a.trend or a.depth
            or a.herd or a.cfm):
        sys.exit(f"--trend-strength is refused except as {commands}: --entry core_trim or organ_scope and --z-mode raw "
                 "only, and with --expression copy, --trend, --depth, --herd and --cfm")
    if a.board not in (None, TEST_BOARDS["heart"]):
        sys.exit(f"--trend-strength is refused on {a.board}: the heart search line builds on the test board only "
                 "(no dry run)")
    if a.entry == "organ_scope":  # file D's line (row R29); row R28 never orders it
        if ruling is None:
            sys.exit(f"--entry organ_scope --trend-strength {raw} is the captain's order of row {HEART_OTS_ROW} alone "
                     f"and needs --captain-ruling {HEART_OTS_ROW}")
        if ruling != HEART_OTS_ROW:
            sys.exit(f"--entry organ_scope --trend-strength {raw} is refused with --captain-ruling {ruling}: only row "
                     f"{HEART_OTS_ROW} orders it")
        if not _heart_d_line_ordered():
            sys.exit(f"--entry organ_scope --trend-strength {raw} is refused: experiments/harness/TRIES.md must hold "
                     f"exactly one row {HEART_OTS_ROW} whose first cell is the captain's order (\"ruling, no run\") "
                     f"naming {HEART_OTS_ENTRY}, {HEART_TS_OPTIONS[0.5]}, {HEART_TS_OPTIONS[1.0]} and "
                     f"`--captain-ruling {HEART_OTS_ROW}`")
        return s
    if ruling is None:
        sys.exit(f"--trend-strength {raw} is the captain's order of row {HEART_TS_ROW} alone and needs "
                 f"--captain-ruling {HEART_TS_ROW}")
    if ruling != HEART_TS_ROW:
        sys.exit(f"--trend-strength {raw} is refused with --captain-ruling {ruling}: only row {HEART_TS_ROW} orders "
                 "it")
    if not _heart_line_ordered():
        sys.exit(f"--trend-strength {raw} is refused: experiments/harness/TRIES.md must hold exactly one row "
                 f"{HEART_TS_ROW} whose first cell is the captain's order (\"ruling, no run\") naming "
                 f"{HEART_TS_OPTIONS[0.5]}, {HEART_TS_OPTIONS[1.0]} and `--captain-ruling {HEART_TS_ROW}`")
    return s


def _trend_strength_step(ctx: harness.BuildContext, strength: float, form: str | None = None, whole=None) -> tuple:
    """Row H28's strength-scaled b0_trend step (experiments/harness/h28_heart_trend_strength.py trend_expression,
    placebo off): `vec.harness._b0_t2(ctx, trend=True)` line by line with the b0_trend_factors mapping times strength
    and momentum_shift's damping (B0's own factor for a lineage the mapping leaves out) times strength. At strength 1
    it is `_b0_t2(ctx, trend=True)` bit for bit. Returns the stage and b0_trend_factors' record. With `form` (row R55's
    trend-form line; `whole` the whole E10.5) the step is that form of it (_trend_form_step); without, as before."""
    if form is not None:
        return _trend_form_step(ctx, strength, form, whole)
    harness._b0_trend_task_check(ctx.task)
    s = float(strength)
    pool, rows, _, coords = harness.draw_copy(ctx)
    start, end = ctx.inputs["start"], pool

    def typed(x, tag: str):
        m = np.asarray(x.obs["merged"]).astype(str)
        return harness._shallow(x, b0_type=np.where(m == "Untracked", f"Untracked@{tag}", m))

    like = typed(pool, "end")[rows]
    logf, rec = harness.b0_trend_factors(ctx, list(like.var_names),
                                         sorted(set(np.asarray(like.obs["b0_type"]).astype(str))))
    mapping = {k: s * v for k, v in logf.items()}
    shifted, _ = t1_scope.momentum_shift(typed(start, "start"), typed(end, "end"), like, harness.B0_DAMPING * s, (),
                                         key="b0_type", seed=[ctx.seed, 3], logf_by_type=mapping)
    return harness.finish(ctx, pool, dense(shifted[:, ctx.board.genes].X), coords), rec


def _trend_strength_record(s: float, rec: dict, entry: str = "core_trim", ruling: str | None = None,
                           order_row: str | None = None, form: str | None = None, r74: str | None = None) -> dict:
    """What a heart search line file's record names: the option, the strength, the file, the ruling row and the order
    text, the command, the step, its budgets (B0's and the test budget, 0.5 x the strength day) and b0_trend's record
    without any cell number (_trend_record). On the organ cut (row R29) the file, row, order, command and step are
    file D's line's; on the full trim they are row R28's as before. At file D's two further strengths (0.75, 1.25)
    they are that line's, with `ruling` the row given on the command line, and the record also holds that row's first
    cell as TRIES.md has it (order_row). At file D's strength 1.5 they are row R52's, with R52's first cell
    (order_row). With `form` they are row R55's trend-form file's (_trend_form_record)."""
    if r74 is not None:  # row R74's file: its own record, the table's strength string verbatim
        return _r74_record(r74, s, rec, order_row)
    if form is not None:
        return _trend_form_record(form, s, rec, order_row)
    if s in HEART_FURTHER_FILES:
        return {"option": "--trend-strength", "strength": s, "file": HEART_FURTHER_FILES[s], "captain_ruling": ruling,
                "order": HEART_FURTHER_ORDER.format(row=ruling), "order_row": order_row,
                "command": HEART_FURTHER_COMMANDS[s].format(row=ruling), "step": HEART_FURTHER_STEP,
                "budget_days": rec["budget_days"], "test_budget_days": s * rec["budget_days"],
                "trend": _trend_record(rec)}
    if s == HEART_OTS15_STRENGTH:
        return {"option": "--trend-strength", "strength": s, "file": HEART_OTS15_FILE,
                "captain_ruling": HEART_OTS15_ROW, "order": HEART_OTS15_ORDER, "order_row": order_row,
                "command": HEART_OTS15_COMMAND, "step": HEART_OTS15_STEP, "budget_days": rec["budget_days"],
                "test_budget_days": s * rec["budget_days"], "trend": _trend_record(rec)}
    if entry == "organ_scope":
        line = (HEART_OTS_FILES[s], HEART_OTS_ROW, HEART_OTS_ORDER, HEART_OTS_COMMANDS[s], HEART_OTS_STEP)
    else:
        line = (HEART_TS_FILES[s], HEART_TS_ROW, HEART_TS_ORDER, HEART_TS_COMMANDS[s], HEART_TS_STEP)
    file, row, order, command, step = line
    return {"option": "--trend-strength", "strength": s, "file": file, "captain_ruling": row,
            "order": order, "command": command, "step": step,
            "budget_days": rec["budget_days"], "test_budget_days": s * rec["budget_days"],
            "trend": _trend_record(rec)}


# The heart's trend-form line (row R55 of experiments/harness/TRIES.md, the captain's order of 2026-11-24, ruling, no
# run: "Yes to OTK125, OTX125, OTW125 and OTP125: my named exception on the heart's scope and size for these four files,
# bound to their screen"; decision log 2026-11-24 (www); the HEART PORTFOLIO design of the step-change round as its
# cold challenge fixes it, OKS125 struck; the screen is row H35, experiments/harness/h35_heart_trend_forms.py):
# OTS125's command (file D at strength 1.25, row R49's) with one change, the direction of the trend step
# (vec.trend_forms): --trend-form category (OTK125), contrast (OTX125), whole (OTW125) or portion (OTP125). Off by
# default; taken only as heart --entry organ_scope --trend-strength 1.25 --trend-form <F> --z-mode raw --captain-ruling
# R55 at the default seed on the test board (no dry run), and only while TRIES.md holds exactly one row R55 whose first
# cell is the captain's order ("ruling, no run") naming `--entry organ_scope`, `--trend-strength 1.25`,
# `--trend-form <F>` and `--captain-ruling R55`, each as a whole backticked string (_heart_tf_ordered). The rows,
# coordinates, cell number, size field and half's rows are OTS125's (its draw, organ-cut pool, count rule, size
# procedure and seed); the counts of every moved cell are drawn afresh. Every command without --trend-form is routed,
# refused and built as before. The files: T2heart_organ_ts125_<tag>.h5ad, then halve.py at seed 1000001 as for OTS125.
HEART_TF_ROW = "R55"
HEART_TF_FORMS = ("category", "contrast", "whole", "portion")  # vec.trend_forms.FORMS, in the design's order
HEART_TF_FILES = {"category": "OTK125", "contrast": "OTX125", "whole": "OTW125", "portion": "OTP125"}
HEART_TF_TAGS = {"category": "tk", "contrast": "tx", "whole": "tw", "portion": "tp"}  # the design's file tags
HEART_TF_STRENGTH_CHOICE = "1.25"  # OTS125's strength, the option's literal value as row R55 names it
HEART_TF_STRENGTH = float(HEART_TF_STRENGTH_CHOICE)
HEART_TF_BASE_FILE = "OTS125"
HEART_TF_COMMANDS = {f: (f"heart --entry organ_scope --trend-strength {HEART_TF_STRENGTH_CHOICE} --trend-form {f} "
                         f"--z-mode raw --captain-ruling {HEART_TF_ROW}") for f in HEART_TF_FORMS}
# the whole backticked strings row R55's first cell must hold for a form's file, besides `--trend-form <F>`
HEART_TF_NAMED = ("--entry organ_scope", f"--trend-strength {HEART_TF_STRENGTH_CHOICE}",
                  f"--captain-ruling {HEART_TF_ROW}")
HEART_TF_ORDER = ("experiments/harness/TRIES.md row R55, the captain's order of 2026-11-24 (ruling, no run: \"Yes to "
                  "OTK125, OTX125, OTW125 and OTP125: my named exception on the heart's scope and size for these four "
                  "files, bound to their screen\"; decision log 2026-11-24 (www); the screen is row H35): file {file}, "
                  "OTS125's command (file D at strength 1.25, row R49's) with one change, the trend step's {form} "
                  "form; nothing else changed")
HEART_TF_STEP = ("row H28's strength-scaled b0_trend step at strength 1.25 (OTS125's) with each lineage's log factor "
                 "formed by the form (vec.trend_forms), before the size step, on OTS125's draw, organ-cut pool, count "
                 "rule, size procedure and seed; the damping and the stream [seed, 3] are OTS125's, and the counts of "
                 "every moved cell are drawn afresh")
HEART_TF_RECORD = {  # each form's own record fields, beside the trend (no cell number in any)
    "category": ("split_lineage", "split_types", "collapsed"),
    "contrast": ("shared_lineages", "shared_norm", "lineages"),
    "whole": ("lineages",),
    "portion": ("lineages", "banked_factor_kept"),
}


def _heart_tf_ordered(form: str, path: Path | None = None) -> bool:
    """Whether experiments/harness/TRIES.md (TRIES_FILE) holds exactly one row R55 whose first cell is the captain's
    order (it starts "R55 (the captain's order" and says "ruling, no run") naming `--entry organ_scope`,
    `--trend-strength 1.25`, `--trend-form <form>` and `--captain-ruling R55`, each as a whole backticked string: a
    backticked string holding one of them inside a longer one, or one it holds, does not name it, and neither does a
    name outside the first cell."""
    lines = [ln for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines()
             if ln.startswith(f"| {HEART_TF_ROW} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    named = set(re.findall(r"`([^`]*)`", first))
    return (first.startswith(f"{HEART_TF_ROW} (the captain's order") and "ruling, no run" in first
            and f"--trend-form {form}" in named and all(s in named for s in HEART_TF_NAMED))


def _trend_form_refusals(a, form: str) -> float:
    """--trend-form (the trend-form line, row R55): its refusals, before anything is read, and before every refusal of
    --trend-strength (whose values and rows stay as they were for every command without --trend-form). Refused except
    as heart --entry organ_scope --trend-strength 1.25 --trend-form <F> --z-mode raw --captain-ruling R55 at the default
    seed on the test board, F one of the four forms, and unless TRIES.md holds exactly one row R55 whose first cell is
    the captain's order naming the form's file (_heart_tf_ordered). Returns 1.25, OTS125's strength."""
    if form not in HEART_TF_FORMS:
        sys.exit(f"--trend-form {form} is refused: it takes {', '.join(HEART_TF_FORMS)} only")
    command = HEART_TF_COMMANDS[form]
    if (getattr(a, "trend_strength", None) != HEART_TF_STRENGTH_CHOICE or a.entry != "organ_scope" or a.z_mode != "raw"
            or a.expression != "b0" or a.trend or a.depth or a.herd or a.cfm):
        sys.exit(f"--trend-form {form} is refused except as {command}: OTS125's line only (--entry organ_scope, "
                 "--trend-strength 1.25, --z-mode raw), and with --expression copy, --trend, --depth, --herd and --cfm")
    if a.board not in (None, TEST_BOARDS["heart"]) or a.seed != harness.BUILD_SEED_OFFSET + 1:
        sys.exit(f"--trend-form {form} is refused on another board or seed: it builds {command} at the default seed "
                 "on the test board only (no dry run)")
    ruling = getattr(a, "captain_ruling", None)
    if ruling is None:
        sys.exit(f"--trend-form {form} is the captain's order of row {HEART_TF_ROW} alone and needs --captain-ruling "
                 f"{HEART_TF_ROW}")
    if ruling != HEART_TF_ROW:
        sys.exit(f"--trend-form {form} is refused with --captain-ruling {ruling}: only row {HEART_TF_ROW} orders it")
    if not _heart_tf_ordered(form):
        sys.exit(f"--trend-form {form} is refused: experiments/harness/TRIES.md must hold exactly one row "
                 f"{HEART_TF_ROW} whose first cell is the captain's order (\"ruling, no run\") naming "
                 + ", ".join(f"`{s}`" for s in (*HEART_TF_NAMED[:2], f"--trend-form {form}"))
                 + f" and `{HEART_TF_NAMED[2]}`, each as a whole backticked string")
    return HEART_TF_STRENGTH


def _trend_form_step(ctx: harness.BuildContext, strength: float, form: str, whole) -> tuple:
    """Row R55's trend-form step: vec.trend_forms.form_step on file D's context (the reference input D's pool, `whole`
    the whole E10.5 as loaded for the organ cut), OTS125's step with the form's one change. An unmapped label refuses,
    as the organ cut's does. Returns the stage and the form's record."""
    from vec import trend_forms

    try:
        stage, _, _, rec = trend_forms.form_step(ctx, strength, form, whole)
    except t2_build.UnmappedT2LabelError as e:
        sys.exit(str(e))
    return stage, rec


def _trend_form_record(form: str, s: float, rec: dict, order_row: str | None) -> dict:
    """What a trend-form file's record names: the option, the form, the strength, the file and its tag, the base file,
    row R55 with its order text and its first cell as TRIES.md holds it (order_row), the command, the step, the form's
    rule, the budgets, b0_trend's record (_trend_record: no cell number), the form's own fields and each type's
    status in momentum_shift (shifted or skipped)."""
    out = {"option": "--trend-form", "form": form, "strength": s, "file": HEART_TF_FILES[form],
           "tag": HEART_TF_TAGS[form], "base_file": HEART_TF_BASE_FILE, "captain_ruling": HEART_TF_ROW,
           "order": HEART_TF_ORDER.format(file=HEART_TF_FILES[form], form=form), "order_row": order_row,
           "command": HEART_TF_COMMANDS[form], "step": HEART_TF_STEP, "rule": rec["rule"],
           "budget_days": rec["budget_days"], "test_budget_days": s * rec["budget_days"],
           "trend": _trend_record(rec["trend"]), **{k: rec[k] for k in HEART_TF_RECORD[form]},
           "types": rec["types"], "shifted": rec["shifted"]}
    if form == "whole":
        out["trend_whole"] = _trend_record(rec["trend_whole"])
    return out


# Row R74 (the captain's order of 2026-11-25, ruling, no run; decision log 2026-11-25 (jjjj)): the heart's refinement
# round on OTK125's command (row R55's category file), four files, each that command with one change, the category
# form's strength, taken from this one literal table in the row's order: the file and its option string as the row
# names them (each a whole backticked string followed by " (FILE)"), each only with `--entry organ_scope`,
# `--z-mode raw` and `--captain-ruling R74`, at the default seed on the test board. The strength is the table's string,
# matched whole ("1.0", "1.5", "0.75", "1.75"; never R55's constant 1.25, HEART_TF_STRENGTH); the record is the file's
# own (its name, strength string, command and order). The launcher experiments/builds/heart_r74_line.sh holds its own
# literal case table, which tests/test_heart_r74_line.py holds equal to this one.
HEART_R74_ROW = "R74"
HEART_R74_TABLE = (
    ("OTK10", "--trend-strength 1.0 --trend-form category"),
    ("OTK15", "--trend-strength 1.5 --trend-form category"),
    ("OTK075", "--trend-strength 0.75 --trend-form category"),
    ("OTK175", "--trend-strength 1.75 --trend-form category"),
)
HEART_R74_FORM = "category"
HEART_R74_BASE_FILE = "OTK125"
HEART_R74_NAMED = ("--entry organ_scope", "--z-mode raw", f"--captain-ruling {HEART_R74_ROW}")


def _r74_strength(name: str, option: str) -> str:
    """The strength string of row R74's file `name`, read from its option string: refuses (ValueError) unless the
    string is exactly `--trend-strength <S> --trend-form category`, S a decimal d.d or d.dd other than R55's 1.25, and
    the name is OTK with S's digits (OTK10 for 1.0, OTK075 for 0.75)."""
    words = option.split(" ")
    if (len(words) != 4 or words[0] != "--trend-strength" or words[2:] != ["--trend-form", HEART_R74_FORM]
            or re.fullmatch(r"[0-9]\.[0-9][0-9]?", words[1]) is None or words[1] == HEART_TF_STRENGTH_CHOICE
            or name != "OTK" + words[1].replace(".", "")):
        raise ValueError(f"row {HEART_R74_ROW}'s file {name} does not match its option string {option}")
    return words[1]


HEART_R74_FILES = {name: (option, _r74_strength(name, option)) for name, option in HEART_R74_TABLE}
HEART_R74_BY_OPTION = {option: name for name, option in HEART_R74_TABLE}
HEART_R74_STEMS = {name: f"T2heart_organ_ts{s.replace('.', '')}_tk" for name, (_, s) in HEART_R74_FILES.items()}
HEART_R74_COMMAND = "heart --entry organ_scope {option} --z-mode raw --captain-ruling " + HEART_R74_ROW
HEART_R74_ORDER = ("experiments/harness/TRIES.md row R74, the captain's order of 2026-11-25 (ruling, no run; decision "
                   "log 2026-11-25 (jjjj); his named exception on the heart's scope and size for each of the four "
                   "files) naming `--entry organ_scope`, `--z-mode raw`, `{option}` ({file}) and `--captain-ruling "
                   "R74`: file {file}, OTK125's command (row R55's category file) with one change, the category "
                   "form's strength {strength} in place of 1.25, then halved by experiments/builds/halve.py at seed "
                   "1000001; nothing else changed")
HEART_R74_STEP = ("row H28's strength-scaled b0_trend step with the category form (vec.trend_forms), OTK125's step, at "
                  "strength {strength} (row R74's table) in place of 1.25: the form's log factors and momentum_shift's "
                  "fallback damping (vec.harness.B0_DAMPING) both multiplied by it, before the size step, on OTK125's "
                  "draw, organ-cut pool, count rule, size procedure, seed and stream [seed, 3]; the counts of every "
                  "moved cell are drawn afresh")


def _heart_r74_ordered(name: str, path: Path | None = None) -> bool:
    """Whether TRIES.md holds exactly one row R74 whose first cell is the captain's order ("ruling, no run") naming
    `--entry organ_scope`, `--z-mode raw`, file `name`'s option string and `--captain-ruling R74`, each as a whole
    backticked string, the option string followed there by " (name)"."""
    lines = [ln for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines()
             if ln.startswith(f"| {HEART_R74_ROW} ")]
    if len(lines) != 1 or name not in HEART_R74_FILES:
        return False
    first = lines[0][2:].split("|", 1)[0]
    named = set(re.findall(r"`([^`]*)`", first))
    option = HEART_R74_FILES[name][0]
    return (first.startswith(f"{HEART_R74_ROW} (the captain's order") and "ruling, no run" in first
            and all(s in named for s in (*HEART_R74_NAMED, option)) and f"`{option}` ({name})" in first)


def _r74_refusals(a, form: str) -> float:
    """--trend-form with --captain-ruling R74: exactly one of row R74's four files, or refused before anything is read.
    Returns that file's strength, the table's (never R55's HEART_TF_STRENGTH)."""
    option = f"--trend-strength {getattr(a, 'trend_strength', None)} --trend-form {form}"
    name = HEART_R74_BY_OPTION.get(option)
    if name is None:
        sys.exit(f"--captain-ruling {HEART_R74_ROW} orders four files, each with exactly its own strings on "
                 f"{HEART_R74_COMMAND.format(option='<OPTION>')}: "
                 + "; ".join(f"{f} `{v[0]}`" for f, v in HEART_R74_FILES.items()))
    command = HEART_R74_COMMAND.format(option=option)
    what = f"row {HEART_R74_ROW}'s file {name}"
    if (a.entry != "organ_scope" or a.z_mode != "raw" or a.expression != "b0" or a.trend or a.depth or a.herd
            or a.cfm):
        sys.exit(f"{what} is refused except as {command}: --entry organ_scope and --z-mode raw only, and with "
                 "--expression copy, --trend, --depth, --herd and --cfm")
    if a.board not in (None, TEST_BOARDS["heart"]) or a.seed != harness.BUILD_SEED_OFFSET + 1:
        sys.exit(f"{what} is refused on another board or seed: it builds {command} at the default seed on the test "
                 "board only (no dry run)")
    if not _heart_r74_ordered(name):
        sys.exit(f"{what} is refused: experiments/harness/TRIES.md must hold exactly one row {HEART_R74_ROW} whose "
                 "first cell is the captain's order (\"ruling, no run\") naming "
                 + ", ".join(f"`{s}`" for s in (*HEART_R74_NAMED[:2], option))
                 + f" and `{HEART_R74_NAMED[2]}`, each as a whole backticked string, `{option}` followed there by "
                 f"\" ({name})\"")
    return float(HEART_R74_FILES[name][1])


def _r74_file(a) -> str | None:
    """Row R74's file this call builds (after _r74_refusals), else None (every other heart command)."""
    if getattr(a, "trend_form", None) is None or getattr(a, "captain_ruling", None) != HEART_R74_ROW:
        return None
    return HEART_R74_BY_OPTION[f"--trend-strength {a.trend_strength} --trend-form {a.trend_form}"]


def _r74_record(name: str, s: float, rec: dict, order_row: str | None) -> dict:
    """Row R74's file's record: the category form's fields (R55's) under the file's own name, strength string,
    command, order and step, and the row's first cell; no cell number."""
    option, raw = HEART_R74_FILES[name]
    return {"option": f"`{option}`", "form": HEART_R74_FORM, "strength": s, "strength_option": raw, "file": name,
            "file_stem": HEART_R74_STEMS[name], "base_file": HEART_R74_BASE_FILE, "captain_ruling": HEART_R74_ROW,
            "order": HEART_R74_ORDER.format(option=option, file=name, strength=raw), "order_row": order_row,
            "command": HEART_R74_COMMAND.format(option=option), "step": HEART_R74_STEP.format(strength=raw),
            "rule": rec["rule"], "budget_days": rec["budget_days"], "test_budget_days": s * rec["budget_days"],
            "trend": _trend_record(rec["trend"]), **{k: rec[k] for k in HEART_TF_RECORD[HEART_R74_FORM]},
            "types": rec["types"], "shifted": rec["shifted"]}


def build_heart(a) -> dict:
    if getattr(a, "insitu_chords", None) is not None or getattr(a, "trend_form", None) == HEART_OIC_FORM:
        return _build_heart_oic(a)  # file OIC125 (row R70): its own path, refused unless exactly its command
    if getattr(a, "run_profile", None) is not None:  # file OXP125 (row R56): its own path; no other command goes there
        return _build_heart_xpc(a)
    strength = _trend_strength_refusals(a)
    form = getattr(a, "trend_form", None)  # row R55's trend-form line (its refusals passed); None on every other
    r74 = _r74_file(a)  # row R74's file (its refusals passed), else None
    # file D's two further strengths and its strength 1.5: the ruling row's first cell for the record, read once its
    # check has passed (and row R55's, whose files take strength 1.25)
    order_row = (_ruling_first_cell(a.captain_ruling)
                 if r74 is not None or strength in HEART_FURTHER_FILES or strength == HEART_OTS15_STRENGTH
                 else None)
    board, dry = _board("heart", a.board)
    genes = board.genes
    if a.trend and (a.entry == "s2" or a.depth or a.herd or a.cfm):
        sys.exit("--trend is refused with --entry s2, --depth, --herd and --cfm (no try read those combinations)")
    copy = a.expression == "copy"
    # The organ cut (HEART_ORGAN_SCOPE; off unless asked for): --entry organ_scope builds that one combination only
    organ = a.entry == "organ_scope"
    # (with --trend-strength, _trend_strength_refusals has already held it to file D's line of row R29, to its two
    # further strengths under row R49 or to its strength 1.5 under row R52)
    if organ and strength is None and (not (copy and a.z_mode == "raw") or a.trend or a.depth or a.herd or a.cfm):
        sys.exit(f"--entry organ_scope is refused except as {HEART_ORGAN_SCOPE} (the copy draw from the organ cut, "
                 "sized in the raw mode to the core forecast times rms(organ pool) / rms(core) on E10.5), and with "
                 "--trend, --depth, --herd and --cfm")
    # The scope bet (HEART_ISO_TRIM; off unless asked for): --entry s2_trim with --expression copy --z-mode raw, the one
    # trimmed copy file. Every other use of --expression copy outside s1 is refused, as before.
    iso_trim = copy and a.entry == "s2_trim" and a.z_mode == "raw"
    if copy and ((a.entry not in ("s1", "core_trim", "organ_scope") and not iso_trim) or a.trend or a.depth or a.herd
                 or a.cfm):
        sys.exit("--expression copy is refused with --entry s2, --trend, --depth, --herd and --cfm, and with --entry "
                 "s2_trim unless --z-mode raw is given (it is the copy baseline's draw from --entry s1's pool, the "
                 f"full E10.5, with no step; or {HEART_ISO_TRIM}: the same draw from the half trim, sized in the raw "
                 "mode for H2)")
    # The full trim (HEART_CORE_TRIM; off unless asked for): --entry core_trim builds that one combination only
    core_trim = a.entry == "core_trim"
    if core_trim and strength is None and (not (copy and a.z_mode == "raw") or a.trend or a.depth or a.herd or a.cfm):
        sys.exit(f"--entry core_trim is refused except as {HEART_CORE_TRIM} (the copy draw from the full trim, sized "
                 "in the raw mode to the core's own forecast), and with --trend, --depth, --herd and --cfm")
    scope_option = HEART_SCOPE_OPTIONS[a.entry] if iso_trim or core_trim else None
    ref = load_stage(DATA / "T2/E10.5.h5ad", genes=genes, drop_duplicates=True)
    start = load_stage(DATA / "T2/E9.5.h5ad", genes=genes, drop_duplicates=True)
    pool, trim = ref, None
    if a.entry == "s2_trim":  # docs/S2_RULES.md (amended 2026-11-02): H2 built as a trim to half E10.5's scope
        core = np.asarray(ref.obs["core"], dtype=bool)
        C = np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3]
        keep = t2_build.radial_trim_mask(C, core, 0.5)
        pool = ref[keep].copy()
        trim = {"rule": "surrounding cells removed farthest in plane from the core first, to half the in-plane scope",
                "kept_share": float(keep.mean()), "xy_scope_before": t2_build.xy_scope(C, core),
                "xy_scope_after": t2_build.xy_scope(C[keep], core[keep])}
    elif core_trim:  # the full trim (no ruling): radial_trim_mask at fraction 0
        core = np.asarray(ref.obs["core"], dtype=bool)
        C = np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3]
        keep = t2_build.radial_trim_mask(C, core, HEART_SCOPE_OPTIONS["core_trim"]["frac"])
        pool = ref[keep].copy()
        trim = {"rule": HEART_SCOPE_OPTIONS["core_trim"]["rule"], "frac": HEART_SCOPE_OPTIONS["core_trim"]["frac"],
                "kept_share": float(keep.mean()), "xy_scope_before": t2_build.xy_scope(C, core),
                "xy_scope_after": t2_build.xy_scope(C[keep], core[keep])}
    elif organ:  # file D: the organ cut of E10.5 (row R21)
        core = np.asarray(ref.obs["core"], dtype=bool)
        C = np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3]
        if t2_build.T2_LABEL_KEY not in ref.obs:
            sys.exit(f"E10.5 has no obs[{t2_build.T2_LABEL_KEY!r}] column: the organ cut needs the released labels")
        try:
            keep, trim = t2_build.organ_scope(C, ref.obs[t2_build.T2_LABEL_KEY].to_numpy(), core)
        except t2_build.UnmappedT2LabelError as e:
            sys.exit(str(e))
        organ_ratio = t2_build.organ_core_rms_ratio(C, keep, core)
        pool = ref[keep].copy()
        # the pool's rows are E10.5's (released data), never a test board's count
        trim.update(pool_rows=int(keep.sum()), pool_rows_of="the released E10.5 after duplicate records are dropped",
                    organ_core_rms_ratio=organ_ratio)
    roles ={"ref": ("T2/E10.5.h5ad", 10.5), "start": ("T2/E9.5.h5ad", 9.5)}
    inputs = {"ref": pool, "start": start}
    if a.trend or strength is not None:  # b0_trend (row H11; row R28's strength): the earlier stages on the board's
        # panel, loaded as the two above, never trimmed
        for role, (f, day) in HEART_TREND_STAGES.items():
            inputs[role] = load_stage(DATA / f, genes=genes, drop_duplicates=True)
            roles[role] = (f, day)
    view = _view("test:T2h", "T2", "extrap", "heart", board.key, roles, "ref", 12.5)
    count = _t2_count(board, pool.n_obs)
    ctx = harness.BuildContext(task=view, board=harness._candidate_board(board), inputs=inputs, seed=a.seed)
    trend_rec: dict = {}
    if copy:  # the copy baseline: B0's own draw (draw_copy, same context, seed and count) through finish, no step
        stage = harness._copy(ctx)
    elif strength is not None:  # the heart search line (row R28; rows R29, R49 and R52 on the organ cut): row H28's
        # step at the strength, before the size step, on this entry's own context (C's or D's); with --trend-form, row
        # R55's form of OTS125's step (the whole E10.5, `ref`, for the forms that read it)
        stage, ts_rec = _trend_strength_step(ctx, strength, form, ref)
    elif a.trend:
        stage = harness._b0_t2(ctx, trend=True, trend_record=trend_rec)  # b0_trend, the harness's B0 otherwise
    else:
        stage = harness._b0_t2(ctx, depth=a.depth, herd=a.herd, cfm=a.cfm)  # the harness's B0 with its guard and checks
    hyp = {"s1": "H1", "s2": "H2", "s2_trim": "H2", "core_trim": t2_build.SIZE_CORE_SCOPE,
           "organ_scope": t2_build.SIZE_ORGAN_SCOPE}[a.entry]
    z_mode = None if a.z_mode == "recorded" else a.z_mode  # None: the mode the size JSON recorded
    # the size step checks its reference against the full E10.5 cloud it measured; it rescales the (trimmed) cells
    if organ:
        coords, size = t2_build.size_organ_scope(stage.coords, SIZES["heart"], organ_ratio, source_stage="E10.5",
                                                 reference_cloud=ref, expect_board="heart")
    elif z_mode is None:
        coords, size = t2_build.size_or_hold(stage.coords, SIZES["heart"], hyp, source_stage="E10.5",
                                             reference_cloud=ref, expect_board="heart")
    else:
        coords, size = t2_build.size_or_hold(stage.coords, SIZES["heart"], hyp, source_stage="E10.5",
                                             reference_cloud=ref, z_mode=z_mode, expect_board="heart")
    _, report = t2_build.build_submission(board, stage.X, coords, genes=list(genes), out=a.out)
    if not report.ok:
        raise SystemExit(f"the built file fails its check: {report.errors}")
    entry = f"copy {a.entry} ({hyp})" if copy else f"B0 {a.entry} ({hyp})"
    if strength is not None:
        entry = f"trend strength {strength:g} {a.entry} ({hyp})"
    if form is not None:
        entry = f"trend strength {strength:g} form {form} {a.entry} ({hyp})"
    return {"entry": entry, "board": board.key, "dry_run": dry, "n_cells": int(stage.n),
            "seed": a.seed, "count_decision": count, "trim": trim, "depth_step": bool(a.depth),
            **({"cfm_step": t2_build.CFM_RULE} if a.cfm else {}),
            **({"trend": _trend_record(trend_rec)} if a.trend else {}),
            **({"trend_strength": _trend_strength_record(strength, ts_rec, a.entry, a.captain_ruling, order_row,
                                                         form, r74)}
               if strength is not None else {}),
            **({"expression": "copy"} if copy else {}),
            **({"z_mode": z_mode} if z_mode is not None else {}),
            # on the organ cut (row R29) scope_option stays file D's own record (_organ_option_record, below)
            **({"scope_option": {**_scope_option_record(a.entry, scope_option), "command": HEART_TS_COMMANDS[strength],
                                 "ruling": HEART_TS_ORDER}}
               if strength is not None and not organ else
               {"scope_option": _scope_option_record(a.entry, scope_option)} if scope_option is not None else {}),
            **({"scope_option": _organ_option_record(trim)} if organ else {}),
            "size": {k: size[k] for k in size if not isinstance(size[k], (dict, list))}}


def _organ_option_record(trim: dict) -> dict:
    """What the organ cut's build record names: the command, the order it was built under, the scope rule, the label
    config and its sha256, L and its percentile, the record fields of the kept share and of the pool's rows (E10.5's,
    released data; no test board's count), the size field and its source, and the halving that follows."""
    return {"command": HEART_ORGAN_SCOPE, "ruling": HEART_ORGAN_RULING, "scope": "organ_scope",
            "scope_rule": trim["rule"], "label_config": trim["label_config"],
            "label_config_sha256": trim["label_config_sha256"], "link_percentile": trim["link_percentile"],
            "link_length": trim["link_length"], "kept_share_field": "trim.kept_share",
            "pool_rows_field": "trim.pool_rows", "size_hypothesis": t2_build.SIZE_ORGAN_SCOPE,
            "size_source_field": t2_build.SIZE_ORGAN_SOURCE, "halving": HEART_HALVING}


def _scope_option_record(entry: str, opt: dict) -> dict:
    """What a scope option's build record names (no ruling; no cell number): its scope rule and trim fraction, the
    record field that holds its kept share (a share of the released E10.5, never a count), the size JSON field its size
    comes from, and the halving that follows as a separate step."""
    return {"command": HEART_ISO_TRIM if entry == "s2_trim" else HEART_CORE_TRIM, "ruling": None,
            "scope_rule": opt["rule"], "trim_fraction": opt["frac"], "kept_share_field": "trim.kept_share",
            "size_hypothesis": opt["hypothesis"], "size_source_field": opt["size_source"], "halving": HEART_HALVING}


TREND_LINEAGE_FIELDS = ("lineage", "stages", "days", "own_factor", "weights", "cosine", "norm_ratio",
                        "mean_abs_logf_trend", "mean_abs_logf_b0")


def _trend_record(rec: dict) -> dict:
    """b0_trend's record (vec.harness.b0_trend_factors) for the build record: the rule, the stage roles and days, the
    transition weights, and per lineage the usable stages, cosine and size ratio (L2 norm, trend / B0) against B0's own
    factor and the two mean absolute log factors, plus the medians. Built from named fields only, so no cell number of
    any stage or pool enters it (the harness's per-stage cell numbers and min_cells are left out), on a dry run as on
    the test board."""
    return {"rule": "vec.harness.B0_TREND_RULE", "stages": list(rec["stages"]), "days": list(rec["days"]),
            "weights": [float(w) for w in harness.trend_transition_weights(rec["days"])],
            "budget_days": rec["budget_days"],
            "lineages": [{k: row[k] for k in TREND_LINEAGE_FIELDS} for row in rec["lineages"]],
            "median_cosine": rec["median_cosine"], "median_norm_ratio": rec["median_norm_ratio"]}


# XPC, file OXP125 (row R56 of experiments/harness/TRIES.md, the captain's order of 2026-11-24, ruling, no run: "... and
# yes to XPC: my named exception on the heart's scope and size for this file and my A15 exception, bound to its
# read."; decision log 2026-11-24 (www)): OTS125's command (file D at trend strength 1.25, row R49) with E10.5's own run
# profile removed. OTS125's draw, organ-cut pool, count rule, trend step (row H28's at 1.25, stream [seed, 3]), size
# procedure and seed; between the trend step and the size step, vec.t2_build.apply_run_profile with the profile row
# H36's test profile wrote (HEART_XPC_PROFILE, its sha256 given on the command line; stream [seed, 37]); then halve.py
# at seed 1000001 as for OTS125. Off by default: build_heart hands this path only a command that gives --run-profile,
# so every other command is exactly as before. Taken only as HEART_XPC_COMMAND on the test board at the default seed,
# every other heart option at its default, and only while TRIES.md holds exactly one row R56 whose first cell is the
# captain's order ("ruling, no run") naming `--entry organ_scope`, `--trend-strength 1.25`, `--run-profile` and
# `--captain-ruling R56` (_captain_order_names), and the profile is the file of that sha256: a passed test profile on
# the board's panel in its order (experiments/heart_profile/xpc_table.py load_profile).
HEART_XPC_ROW = "R56"
HEART_XPC_CHOICE = "1.25"  # the trend strength's literal value, OTS125's
HEART_XPC_STRENGTH = float(HEART_XPC_CHOICE)
HEART_XPC_NAMES = (HEART_OTS_ENTRY, f"`--trend-strength {HEART_XPC_CHOICE}`", "`--run-profile`")
HEART_XPC_FILE = "OXP125"
HEART_XPC_SHA256 = re.compile(r"[0-9a-f]{64}")
HEART_XPC_PROFILE = REPO_ROOT / "experiments" / "heart_profile" / "results" / "xpc_profile.json"
HEART_XPC_COMMAND = (f"heart --entry organ_scope --trend-strength {HEART_XPC_CHOICE} --run-profile {{sha256}} "
                     f"--z-mode raw --captain-ruling {HEART_XPC_ROW}")  # .format(sha256=...)
HEART_XPC_GIVEN = ("which", "out", "board", "entry", "trend_strength", "z_mode", "captain_ruling", "run_profile")
HEART_XPC_ORDER = ("experiments/harness/TRIES.md row R56, the captain's order of 2026-11-24 (ruling, no run: \"... "
                   "and yes to XPC: my named exception on the heart's scope and size for this file and my A15 "
                   "exception, bound to its read.\"; decision log 2026-11-24 (www)): file OXP125, OTS125's command "
                   "with E10.5's run profile removed, bound to row H36's read")
HEART_XPC_STEP = ("row H28's strength-scaled b0_trend step at 1.25 as OTS125 builds it, then the run-profile step "
                  "(vec.t2_build.apply_run_profile: exp(-delta) per gene on every drawn cell, keep_depth, stream "
                  "[seed, 37]), before the size step, on D's draw, organ-cut pool, count rule, size procedure and seed "
                  "(OTS125's)")


def _xpc_refusals(a) -> float:
    """--run-profile (file OXP125, row R56): its refusals, before anything is read. Refused unless the value is a
    sha256 (64 lowercase hexadecimal digits), the command is HEART_XPC_COMMAND on the test board at the default seed
    with every other heart option at its default (any option added later included), --captain-ruling is R56 and
    TRIES.md holds exactly one row R56 whose first cell is the captain's order naming HEART_XPC_NAMES and
    `--captain-ruling R56`. Returns the strength, 1.25."""
    raw = a.run_profile
    command = HEART_XPC_COMMAND.format(sha256="<sha256>")
    if not isinstance(raw, str) or HEART_XPC_SHA256.fullmatch(raw) is None:
        sys.exit(f"--run-profile {raw} is refused: it takes the sha256 (64 lowercase hexadecimal digits) of "
                 "experiments/heart_profile/results/xpc_profile.json, row H36's test profile")
    defaults = vars(build_parser().parse_args(["heart", "--entry", "organ_scope", "--out", "x"]))
    other = sorted(k for k, v in vars(a).items()
                   if k not in HEART_XPC_GIVEN and (k not in defaults or v != defaults[k]))
    if (a.entry != "organ_scope" or a.trend_strength != HEART_XPC_CHOICE or a.z_mode != "raw" or other
            or a.board not in (None, TEST_BOARDS["heart"])):
        sys.exit(f"--run-profile is refused except as {command}, on the test board at the default seed with every "
                 "other heart option at its default (file OXP125, row R56)")
    if a.captain_ruling is None:
        sys.exit(f"--run-profile is the captain's order of row {HEART_XPC_ROW} alone and needs --captain-ruling "
                 f"{HEART_XPC_ROW}")
    if a.captain_ruling != HEART_XPC_ROW:
        sys.exit(f"--run-profile is refused with --captain-ruling {a.captain_ruling}: only row {HEART_XPC_ROW} orders "
                 "it")
    if not all(_captain_order_names(HEART_XPC_ROW, s) for s in HEART_XPC_NAMES):
        sys.exit(f"--run-profile is refused: experiments/harness/TRIES.md must hold exactly one row {HEART_XPC_ROW} "
                 f"whose first cell is the captain's order (\"ruling, no run\") naming {', '.join(HEART_XPC_NAMES)} "
                 f"and `--captain-ruling {HEART_XPC_ROW}`")
    return HEART_XPC_STRENGTH


def _xpc_profile(sha256: str, genes) -> np.ndarray:
    """The run profile OXP125 applies: HEART_XPC_PROFILE, the file of `sha256`, a passed test profile naming the
    board's genes in their order (xpc_table.load_profile); refused otherwise, before any stage is read."""
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "heart_profile"))
    import xpc_table  # noqa: E402 (only this path reads it)

    if not HEART_XPC_PROFILE.is_file():
        sys.exit(f"--run-profile {sha256} is refused: {HEART_XPC_PROFILE} does not exist (row H36's test profile "
                 "writes it)")
    try:
        return xpc_table.load_profile(HEART_XPC_PROFILE, sha256, list(genes))
    except ValueError as e:
        sys.exit(f"--run-profile {sha256} is refused: {e}")


def _run_profile_step(ctx: harness.BuildContext, strength: float, delta: np.ndarray) -> tuple:
    """OXP125's expression: _trend_strength_step's own lines at `strength` (OTS125's step: the same draw, mapping,
    damping and stream [seed, 3]), then vec.t2_build.apply_run_profile on the shifted counts in the board's gene order
    at stream [seed, 37], then harness.finish on the step's own draw and coordinates. A zero profile gives
    _trend_strength_step's stage bit for bit (tested). Returns the stage and b0_trend_factors' record."""
    harness._b0_trend_task_check(ctx.task)
    s = float(strength)
    pool, rows, _, coords = harness.draw_copy(ctx)
    start, end = ctx.inputs["start"], pool

    def typed(x, tag: str):
        m = np.asarray(x.obs["merged"]).astype(str)
        return harness._shallow(x, b0_type=np.where(m == "Untracked", f"Untracked@{tag}", m))

    like = typed(pool, "end")[rows]
    logf, rec = harness.b0_trend_factors(ctx, list(like.var_names),
                                         sorted(set(np.asarray(like.obs["b0_type"]).astype(str))))
    mapping = {k: s * v for k, v in logf.items()}
    shifted, _ = t1_scope.momentum_shift(typed(start, "start"), typed(end, "end"), like, harness.B0_DAMPING * s, (),
                                         key="b0_type", seed=[ctx.seed, 3], logf_by_type=mapping)
    _, X, _ = t2_build.apply_run_profile(t1_scope.counts_of(shifted, list(ctx.board.genes)), delta,
                                         [ctx.seed, t2_build.RUN_PROFILE_STREAM])
    return harness.finish(ctx, pool, X, coords), rec


def _build_heart_xpc(a) -> dict:
    """File OXP125 (row R56): build_heart's organ-cut path at strength 1.25, line by line as OTS125 builds it, with
    the run-profile step (_run_profile_step) in place of the trend step alone. Its record is OTS125's, with its
    trend_strength block naming OXP125, row R56 and its first cell, and a run_profile block (the option, the sha256,
    the file, the rule; no cell number)."""
    strength = _xpc_refusals(a)
    order_row = _ruling_first_cell(HEART_XPC_ROW)
    board, dry = _board("heart", a.board)
    genes = board.genes
    delta = _xpc_profile(a.run_profile, genes)
    ref = load_stage(DATA / "T2/E10.5.h5ad", genes=genes, drop_duplicates=True)
    start = load_stage(DATA / "T2/E9.5.h5ad", genes=genes, drop_duplicates=True)
    core = np.asarray(ref.obs["core"], dtype=bool)
    C = np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3]
    if t2_build.T2_LABEL_KEY not in ref.obs:
        sys.exit(f"E10.5 has no obs[{t2_build.T2_LABEL_KEY!r}] column: the organ cut needs the released labels")
    try:
        keep, trim = t2_build.organ_scope(C, ref.obs[t2_build.T2_LABEL_KEY].to_numpy(), core)
    except t2_build.UnmappedT2LabelError as e:
        sys.exit(str(e))
    organ_ratio = t2_build.organ_core_rms_ratio(C, keep, core)
    pool = ref[keep].copy()
    trim.update(pool_rows=int(keep.sum()), pool_rows_of="the released E10.5 after duplicate records are dropped",
                organ_core_rms_ratio=organ_ratio)
    roles = {"ref": ("T2/E10.5.h5ad", 10.5), "start": ("T2/E9.5.h5ad", 9.5)}
    inputs = {"ref": pool, "start": start}
    for role, (f, day) in HEART_TREND_STAGES.items():
        inputs[role] = load_stage(DATA / f, genes=genes, drop_duplicates=True)
        roles[role] = (f, day)
    view = _view("test:T2h", "T2", "extrap", "heart", board.key, roles, "ref", 12.5)
    count = _t2_count(board, pool.n_obs)
    ctx = harness.BuildContext(task=view, board=harness._candidate_board(board), inputs=inputs, seed=a.seed)
    stage, ts_rec = _run_profile_step(ctx, strength, delta)
    coords, size = t2_build.size_organ_scope(stage.coords, SIZES["heart"], organ_ratio, source_stage="E10.5",
                                             reference_cloud=ref, expect_board="heart")
    _, report = t2_build.build_submission(board, stage.X, coords, genes=list(genes), out=a.out)
    if not report.ok:
        raise SystemExit(f"the built file fails its check: {report.errors}")
    command = HEART_XPC_COMMAND.format(sha256=a.run_profile)
    return {"entry": f"trend strength {strength:g} organ_scope ({t2_build.SIZE_ORGAN_SCOPE}) run profile",
            "board": board.key, "dry_run": dry, "n_cells": int(stage.n), "seed": a.seed, "count_decision": count,
            "trim": trim, "depth_step": False,
            "trend_strength": {"option": "--trend-strength", "strength": strength, "file": HEART_XPC_FILE,
                               "captain_ruling": HEART_XPC_ROW, "order": HEART_XPC_ORDER, "order_row": order_row,
                               "command": command, "step": HEART_XPC_STEP, "budget_days": ts_rec["budget_days"],
                               "test_budget_days": strength * ts_rec["budget_days"], "trend": _trend_record(ts_rec)},
            "run_profile": {"option": "--run-profile", "sha256": a.run_profile, "file": HEART_XPC_FILE,
                            "profile": "experiments/heart_profile/results/xpc_profile.json", "read": "row H36",
                            "rule": t2_build.RUN_PROFILE_RULE, "stream": t2_build.RUN_PROFILE_STREAM,
                            "captain_ruling": HEART_XPC_ROW, "command": command},
            "z_mode": "raw", "scope_option": _organ_option_record(trim),
            "size": {k: size[k] for k in size if not isinstance(size[k], (dict, list))}}


# File OIC125 (row R70 of experiments/harness/TRIES.md, the captain's order of 2026-11-25, ruling, no run: his A15
# exception for OIC125's outside sources and his named exception on the heart's scope and size for OIC125; read by row
# H38). OTK125's command (the category form at strength 1.25 on file D's organ cut) with each chamber category's in-situ
# chord added where it exists (vec.trend_forms.insitu_step; HEART_MOON1_design.md P1 with HEART_MOON1_challenge.md's
# fixes). Expression only: draw, rows, coordinates, count, size step and half are OTK125's. FIX 6's faults (a P-carrying
# type skipped, an applied log factor off P's entries not OTK125's) refuse the build.
HEART_OIC_ROW = "R70"
HEART_OIC_FORM = "insitu"
HEART_OIC_CHOICE = "1.25"
HEART_OIC_STRENGTH = float(HEART_OIC_CHOICE)
HEART_OIC_FILE = "OIC125"
HEART_OIC_CHORDS = REPO_ROOT / "experiments" / "heart_ext" / "results" / "chords_E10.5.json"
HEART_OIC_COMMAND = (f"heart --entry organ_scope --trend-strength {HEART_OIC_CHOICE} --trend-form {HEART_OIC_FORM} "
                     f"--insitu-chords {{sha256}} --z-mode raw --captain-ruling {HEART_OIC_ROW}")
HEART_OIC_GIVEN = ("which", "out", "board", "entry", "trend_strength", "z_mode", "captain_ruling", "insitu_chords",
                   "trend_form")
HEART_OIC_ORDER = ("experiments/harness/TRIES.md row R70, the captain's order of 2026-11-25 (ruling, no run: \"my A15 "
                   "exception for OIC125's outside sources and my named exception on the heart's scope and size for "
                   "OIC125\"; decision log 2026-11-25 (qqq)): file OIC125, OTK125's command with the in-situ chords "
                   "of row H38's G0 added; nothing else changed")


def _oic_read_passed(sha256: str, path: Path | None = None) -> list[str]:
    """Why row H38's outcome cell does not allow OIC125 with this chord file (empty when it does)."""
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "harness"))
    import h38_row

    try:
        text = (path or TRIES_FILE).read_text(encoding="utf-8")
    except OSError:
        return ["TRIES.md cannot be read"]
    got, why = h38_row.passed_sha256(text)
    if not why and got != sha256:
        why = ["row H38's G0 statement records another sha256 for chords_E10.5.json"]
    return why


def _oic_refusals(a) -> float:
    """The refusals of OIC125's command: exactly HEART_OIC_COMMAND, row R70's order, row H38's PASS."""
    raw = getattr(a, "insitu_chords", None)
    command = HEART_OIC_COMMAND.format(sha256="<sha256>")
    if not isinstance(raw, str) or HEART_XPC_SHA256.fullmatch(raw) is None:
        sys.exit(f"--insitu-chords {raw} is refused: it takes the sha256 (64 lowercase hexadecimal digits) of "
                 "experiments/heart_ext/results/chords_E10.5.json, row H38's G0 chords at the test")
    defaults = vars(build_parser().parse_args(["heart", "--entry", "organ_scope", "--out", "x"]))
    other = sorted(k for k, v in vars(a).items()
                   if k not in HEART_OIC_GIVEN and (k not in defaults or v != defaults[k]))
    if (getattr(a, "trend_form", None) != HEART_OIC_FORM or a.entry != "organ_scope"
            or a.trend_strength != HEART_OIC_CHOICE or a.z_mode != "raw" or other
            or a.board not in (None, TEST_BOARDS["heart"])):
        sys.exit(f"--insitu-chords and --trend-form {HEART_OIC_FORM} are refused except as {command}, on the test "
                 "board at the default seed with every other heart option at its default (file OIC125, row R70)")
    if a.captain_ruling is None:
        sys.exit(f"file {HEART_OIC_FILE} is the captain's order of row {HEART_OIC_ROW} alone and needs "
                 f"--captain-ruling {HEART_OIC_ROW}")
    if a.captain_ruling != HEART_OIC_ROW:
        sys.exit(f"file {HEART_OIC_FILE} is refused with --captain-ruling {a.captain_ruling}: only row "
                 f"{HEART_OIC_ROW} orders it")
    if not _captain_order_names(HEART_OIC_ROW, HEART_OIC_FILE):
        sys.exit(f"file {HEART_OIC_FILE} is refused: experiments/harness/TRIES.md must hold exactly one row "
                 f"{HEART_OIC_ROW} whose first cell is the captain's order (\"ruling, no run\") naming "
                 f"{HEART_OIC_FILE} and `--captain-ruling {HEART_OIC_ROW}`")
    why = _oic_read_passed(raw)
    if why:
        sys.exit(f"file {HEART_OIC_FILE} is refused: its gates have not read PASS ({'; '.join(why)})")
    return HEART_OIC_STRENGTH


def _oic_chords(sha256: str, genes) -> dict:
    """chords_E10.5.json as category -> P over `genes` (NaN where none); refused unless it is the file of `sha256`
    and G0's test chords for OIC125."""
    if not HEART_OIC_CHORDS.is_file():
        sys.exit(f"--insitu-chords {sha256} is refused: {HEART_OIC_CHORDS} does not exist (row H38's G0 writes it)")
    raw = HEART_OIC_CHORDS.read_bytes()
    if hashlib.sha256(raw).hexdigest() != sha256:
        sys.exit(f"--insitu-chords {sha256} is refused: {HEART_OIC_CHORDS} is not the file of that sha256")
    doc = json.loads(raw)
    if doc.get("row") != "H38" or doc.get("ref") != "E10.5" or doc.get("file") != HEART_OIC_FILE:
        sys.exit(f"--insitu-chords {sha256} is refused: the file is not row H38's test chords for {HEART_OIC_FILE}")
    idx = {g: i for i, g in enumerate(doc["genes"])}
    out = {}
    for k, v in doc["P"].items():
        if len(v) != len(idx):
            sys.exit(f"--insitu-chords {sha256} is refused: the chord of {k} does not cover the file's genes")
        out[k] = np.array([np.nan if g not in idx or v[idx[g]] is None else float(v[idx[g]]) for g in genes])
    return out


def _build_heart_oic(a) -> dict:
    """File OIC125: file D's organ cut of E10.5, OTK125's draw and category step with the in-situ chords added, then
    file D's size step (the half is halve.py's, as every heart file)."""
    from vec import trend_forms

    strength = _oic_refusals(a)
    order_row = _ruling_first_cell(HEART_OIC_ROW)
    board, dry = _board("heart", a.board)
    genes = board.genes
    ref = load_stage(DATA / "T2/E10.5.h5ad", genes=genes, drop_duplicates=True)
    start = load_stage(DATA / "T2/E9.5.h5ad", genes=genes, drop_duplicates=True)
    core = np.asarray(ref.obs["core"], dtype=bool)
    C = np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3]
    if t2_build.T2_LABEL_KEY not in ref.obs:
        sys.exit(f"E10.5 has no obs[{t2_build.T2_LABEL_KEY!r}] column: the organ cut needs the released labels")
    try:
        keep, trim = t2_build.organ_scope(C, ref.obs[t2_build.T2_LABEL_KEY].to_numpy(), core)
    except t2_build.UnmappedT2LabelError as e:
        sys.exit(str(e))
    organ_ratio = t2_build.organ_core_rms_ratio(C, keep, core)
    pool = ref[keep].copy()
    trim.update(pool_rows=int(keep.sum()), pool_rows_of="the released E10.5 after duplicate records are dropped",
                organ_core_rms_ratio=organ_ratio)
    roles = {"ref": ("T2/E10.5.h5ad", 10.5), "start": ("T2/E9.5.h5ad", 9.5)}
    inputs = {"ref": pool, "start": start}
    for role, (f, day) in HEART_TREND_STAGES.items():
        inputs[role] = load_stage(DATA / f, genes=genes, drop_duplicates=True)
        roles[role] = (f, day)
    view = _view("test:T2h", "T2", "extrap", "heart", board.key, roles, "ref", 12.5)
    count = _t2_count(board, pool.n_obs)
    ctx = harness.BuildContext(task=view, board=harness._candidate_board(board), inputs=inputs, seed=a.seed)
    chords = _oic_chords(a.insitu_chords, list(pool.var_names))
    try:
        stage, _, _, rec = trend_forms.insitu_step(ctx, strength, chords)
    except (ValueError, t2_build.UnmappedT2LabelError) as e:
        sys.exit(f"file {HEART_OIC_FILE} is refused: {e}")
    coords, size = t2_build.size_organ_scope(stage.coords, SIZES["heart"], organ_ratio, source_stage="E10.5",
                                             reference_cloud=ref, expect_board="heart")
    _, report = t2_build.build_submission(board, stage.X, coords, genes=list(genes), out=a.out)
    if not report.ok:
        raise SystemExit(f"the built file fails its check: {report.errors}")
    command = HEART_OIC_COMMAND.format(sha256=a.insitu_chords)
    return {"entry": f"trend strength {strength:g} form {HEART_OIC_FORM} organ_scope ({t2_build.SIZE_ORGAN_SCOPE})",
            "board": board.key, "dry_run": dry, "n_cells": int(stage.n), "seed": a.seed, "count_decision": count,
            "trim": trim, "depth_step": False,
            "trend_strength": {"option": "--trend-form", "form": HEART_OIC_FORM, "strength": strength,
                               "file": HEART_OIC_FILE, "base_file": "OTK125", "captain_ruling": HEART_OIC_ROW,
                               "order": HEART_OIC_ORDER, "order_row": order_row, "command": command,
                               "step": trend_forms.INSITU_RULE, "budget_days": rec["budget_days"],
                               "test_budget_days": strength * rec["budget_days"], "trend": _trend_record(rec["trend"]),
                               "types": rec["types"], "shifted": rec["shifted"], "carried": rec["carried"]},
            "insitu_chords": {"option": "--insitu-chords", "sha256": a.insitu_chords, "file": HEART_OIC_FILE,
                              "chords": "experiments/heart_ext/results/chords_E10.5.json", "read": "row H38",
                              "captain_ruling": HEART_OIC_ROW, "command": command},
            "z_mode": "raw", "scope_option": _organ_option_record(trim),
            "size": {k: size[k] for k in size if not isinstance(size[k], (dict, list))}}


def floor_copy_coords(stage_file: str) -> np.ndarray:
    """The floor copy's coordinates as the emulator's bake protocol draws it on a board whose floor is a copy of
    `stage_file` (a T2 path under data/raw): vec.evaluate.subsample of the stage as released (no duplicate drop), at
    SUBSAMPLE_FRAC and seed t2_build.COPY_PARITY_COPY_SEED, rows in released order, float32 as a file holds them."""
    a = read_stage(DATA / stage_file)
    C = np.asarray(a.obsm[COORD_KEY], dtype=np.float32)[:, :3]
    return subsample(Stage(X=np.zeros((len(C), 0), dtype=np.float32), coords=C), SUBSAMPLE_FRAC,
                     t2_build.COPY_PARITY_COPY_SEED).coords


def _series_back_ordered() -> bool:
    """Whether experiments/harness/TRIES.md holds exactly one row R48 whose first cell is the captain's order ("ruling,
    no run") naming `--expression series_back` and `--captain-ruling R48` (_captain_order_names)."""
    return _captain_order_names(EMBRYO_OVSB_ROW, EMBRYO_OVSB_OPTION)


def _series_back_refusals(a) -> bool:
    """--expression series_back's refusals (row R48), before any board or stage is read. Returns whether the option is
    given; without it every command is what it was (embryo --captain-ruling alone is refused)."""
    ruling = getattr(a, "captain_ruling", None)
    if a.expression != "series_back":
        # --neighbourhood run_sorted (row E57's file) takes --captain-ruling too, with its own refusals; without it,
        # every command and its refusal are what they were
        if ruling is not None and getattr(a, "neighbourhood", NEIGHBOURHOODS[0]) == NEIGHBOURHOODS[0]:
            sys.exit(f"embryo --captain-ruling goes with --expression series_back only ({EMBRYO_OVSB_COMMAND})")
        return False
    o_command = (a.method == "fgw_shape0" and a.positions == t2_build.POSITIONS[1] and a.shape_interp == "overlay"
                 and a.shape_interp_test_ok and not a.herd and a.draw == t2_build.DRAWS[0]
                 and a.move_coupling == t2_build.MOVE_COUPLINGS[0] and a.row_order == t2_build.ROW_ORDERS[0]
                 and not a.row_order_test_ok and not getattr(a, "copy_rows", False)
                 and not getattr(a, "copy_rows_test_ok", False) and not getattr(a, "copy_rows_size", False)
                 and not getattr(a, "copy_rows_size_test_ok", False) and a.seed == harness.BUILD_SEED_OFFSET + 1
                 and a.board in (None, TEST_BOARDS["embryo"]))
    if not o_command:
        sys.exit(f"--expression series_back is refused except on O's command with it, {EMBRYO_OVSB_COMMAND}: the "
                 "test board, the default seed and nothing else changed (row R48)")
    if ruling is None:
        sys.exit(f"--expression series_back is the captain's order of row {EMBRYO_OVSB_ROW} alone and needs "
                 f"--captain-ruling {EMBRYO_OVSB_ROW}")
    if ruling != EMBRYO_OVSB_ROW:
        sys.exit(f"--expression series_back is refused with --captain-ruling {ruling}: only row {EMBRYO_OVSB_ROW} "
                 "orders it")
    if not _series_back_ordered():
        sys.exit(f"--expression series_back is refused: experiments/harness/TRIES.md must hold exactly one row "
                 f"{EMBRYO_OVSB_ROW} whose first cell is the captain's order (\"ruling, no run\") naming "
                 f"{EMBRYO_OVSB_OPTION} and `--captain-ruling {EMBRYO_OVSB_ROW}`")
    return True


def _mix_lean_ordered() -> bool:
    """Whether experiments/harness/TRIES.md holds exactly one row R59 whose first cell is the captain's order ("ruling,
    no run") naming `--mix-lean later` and `--captain-ruling R59` (_captain_order_names)."""
    return _captain_order_names(EMBRYO_OT23_ROW, EMBRYO_OT23_OPTION)


def _is_ot23_command(a) -> bool:
    """Whether the parsed command is EMBRYO_OT23_COMMAND (O's command with --mix-lean later --captain-ruling R59) and
    nothing else: every option the embryo parser knows holds the value that command gives it (the test board, named or
    left out, and the default seed), --out aside. An option added to the parser later is held at its default too."""
    want = vars(build_parser().parse_args([*EMBRYO_OT23_COMMAND.split(), "--out", str(a.out)]))
    got = dict(vars(a))
    if got.get("board") == TEST_BOARDS["embryo"]:
        got["board"] = None
    return got == want


def _mix_lean_refusals(a) -> bool:
    """--mix-lean's refusals (row R59), before any board or stage is read. Returns whether the option is given; without
    it every command is what it was, and --captain-ruling R59 is refused (it goes with --mix-lean later only)."""
    lean, ruling = getattr(a, "mix_lean", None), getattr(a, "captain_ruling", None)
    if lean is None:
        if ruling == EMBRYO_OT23_ROW:
            sys.exit(f"embryo --captain-ruling {EMBRYO_OT23_ROW} goes with --mix-lean later only "
                     f"({EMBRYO_OT23_COMMAND})")
        return False
    if lean not in MIX_LEANS:  # the parser's choices refuse it first; kept for a caller that skips the parser
        sys.exit(f"--mix-lean takes {', '.join(MIX_LEANS)} only (row {EMBRYO_OT23_ROW})")
    if ruling is None:
        sys.exit(f"--mix-lean later is the captain's order of row {EMBRYO_OT23_ROW} alone and needs --captain-ruling "
                 f"{EMBRYO_OT23_ROW}")
    if ruling != EMBRYO_OT23_ROW:
        sys.exit(f"--mix-lean later is refused with --captain-ruling {ruling}: only row {EMBRYO_OT23_ROW} orders it")
    if not _is_ot23_command(a):
        sys.exit(f"--mix-lean later is refused except on O's command with it, {EMBRYO_OT23_COMMAND}: the test board, "
                 f"the default seed and nothing else changed (row {EMBRYO_OT23_ROW})")
    if not _mix_lean_ordered():
        sys.exit(f"--mix-lean later is refused: experiments/harness/TRIES.md must hold exactly one row "
                 f"{EMBRYO_OT23_ROW} whose first cell is the captain's order (\"ruling, no run\") naming "
                 f"{EMBRYO_OT23_OPTION} and `--captain-ruling {EMBRYO_OT23_ROW}`")
    return True


def _r75_row_names(name: str, path: Path | None = None) -> bool:
    """True when the ledger holds exactly one row R75 whose first cell starts "R75 (the captain's order", says "ruling,
    no run" and names the file's option string and `--captain-ruling R75` each as a whole backticked string, the string
    followed there by " (<file>)" as the row pairs each string with its file (row R73's check)."""
    lines = [ln for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines()
             if ln.startswith(f"| {EMBRYO_R75_ROW} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    option = EMBRYO_R75_FILES[name][0]
    named = set(re.findall(r"`([^`]*)`", first))
    return (first.startswith(f"{EMBRYO_R75_ROW} (the captain's order") and "ruling, no run" in first
            and option in named and f"--captain-ruling {EMBRYO_R75_ROW}" in named and f"`{option}` ({name})" in first)


def _r75_file(a) -> str | None:
    """The file of row R75's table whose command (EMBRYO_R75_COMMAND: O's command with the file's option string and
    --captain-ruling R75) is the parsed command and nothing else: every option the embryo parser knows holds the value
    that command gives it (the test board, named or left out, and the default seed), --out aside, as _is_ot23_command
    holds row R59's. None when no file's is."""
    got = dict(vars(a))
    if got.get("board") == TEST_BOARDS["embryo"]:
        got["board"] = None
    for name, (option, *_) in EMBRYO_R75_FILES.items():
        want = vars(build_parser().parse_args([*EMBRYO_R75_COMMAND.format(option=option).split(), "--out", str(a.out)]))
        if got == want:
            return name
    return None


def _r75_refusals(a) -> dict | None:
    """Row R75's refusals, the first of build_embryo's, before any board or stage is read and before the older rows'
    refusals read the ruling. With --captain-ruling R75 the command must be one file's of EMBRYO_R75_TABLE exactly
    (_r75_file) and the ledger's row R75 must name that file's string (_r75_row_names); a lean only row R75 has
    (EMBRYO_R75_NEW_LEANS) is refused under any other ruling or none, and --expression depth with --shape-interp overlay
    without a ruling. Returns the file's values ({file, option, expression, mix_lean, t_mix}, t_mix the exact Fraction
    or None) or None: every other command goes on to the older refusals as it did."""
    ruling, lean = getattr(a, "captain_ruling", None), getattr(a, "mix_lean", None)
    if ruling != EMBRYO_R75_ROW:
        if lean in EMBRYO_R75_NEW_LEANS:
            if ruling is None:
                sys.exit(f"--mix-lean {lean} is the captain's order of row {EMBRYO_R75_ROW} alone and needs "
                         f"--captain-ruling {EMBRYO_R75_ROW}")
            sys.exit(f"--mix-lean {lean} is refused with --captain-ruling {ruling}: only row {EMBRYO_R75_ROW} "
                     "orders it")
        if (ruling is None and getattr(a, "expression", None) == EMBRYO_R75_DEPTH
                and getattr(a, "shape_interp", None) == "overlay"):
            sys.exit(f"--expression depth with --shape-interp overlay is row {EMBRYO_R75_ROW}'s file alone and needs "
                     f"--captain-ruling {EMBRYO_R75_ROW}: "
                     + EMBRYO_R75_COMMAND.format(option=EMBRYO_R75_FILES['T2embryo_s2x_ovd.h5ad'][0]))
        return None
    name = _r75_file(a)
    if name is None:
        sys.exit(f"--captain-ruling {EMBRYO_R75_ROW} orders five files, each only as O's command with exactly its own "
                 f"option string, {EMBRYO_R75_COMMAND.format(option='<OPTION>')}, on the test board at the default "
                 "seed and nothing else changed: "
                 + "; ".join(f"{f} `{v[0]}`" for f, v in EMBRYO_R75_FILES.items()))
    option, expression, mix_lean, t_mix = EMBRYO_R75_FILES[name]
    if not _r75_row_names(name):
        sys.exit(f"row {EMBRYO_R75_ROW}'s file {name} is refused: experiments/harness/TRIES.md must hold exactly one "
                 f"row {EMBRYO_R75_ROW} whose first cell is the captain's order (\"ruling, no run\") naming `{option}` "
                 f"({name}) and `--captain-ruling {EMBRYO_R75_ROW}`, each as a whole backticked string")
    return {"file": name, "option": option, "expression": expression, "mix_lean": mix_lean, "t_mix": t_mix}


def _r78_row_names(name: str, path: Path | None = None) -> bool:
    """Whether experiments/harness/TRIES.md (or `path`) holds exactly one row R78 whose first cell is the captain's
    order ("ruling, no run") naming the file's option string and `--captain-ruling R78`, each as a whole backticked
    string, the option string followed there by " (<file>)", as the row pairs each string with its file (row R75's
    check, on row R78)."""
    lines = [ln for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines()
             if ln.startswith(f"| {EMBRYO_R78_ROW} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    option = EMBRYO_R78_FILES[name][0]
    named = set(re.findall(r"`([^`]*)`", first))
    return (first.startswith(f"{EMBRYO_R78_ROW} (the captain's order") and "ruling, no run" in first
            and option in named and f"--captain-ruling {EMBRYO_R78_ROW}" in named and f"`{option}` ({name})" in first)


def _r78_file(a) -> str | None:
    """The row R78 file whose command the parsed arguments are, exactly: O's command with the file's one option string
    and --captain-ruling R78, nothing else changed (the test board named or not, the default seed); else None."""
    got = dict(vars(a))
    if got.get("board") == TEST_BOARDS["embryo"]:
        got["board"] = None
    for name, (option, *_) in EMBRYO_R78_FILES.items():
        want = vars(build_parser().parse_args([*EMBRYO_R78_COMMAND.format(option=option).split(), "--out", str(a.out)]))
        if got == want:
            return name
    return None


def _r78_refusals(a) -> dict | None:
    """Row R78's check, the first of build_embryo's. None for a command that is not row R78's (no --captain-ruling R78
    and no lean only row R78 has); for exactly one of its two commands, while row R78's first cell names it, the file's
    values (file, option, expression, lean and exact t_mix); anything else under --captain-ruling R78, and a lean only
    row R78 has under any other ruling or none, exits here, before anything is read."""
    ruling, lean = getattr(a, "captain_ruling", None), getattr(a, "mix_lean", None)
    if ruling != EMBRYO_R78_ROW:
        if lean in EMBRYO_R78_NEW_LEANS:
            if ruling is None:
                sys.exit(f"--mix-lean {lean} is the captain's order of row {EMBRYO_R78_ROW} alone and needs "
                         f"--captain-ruling {EMBRYO_R78_ROW}")
            sys.exit(f"--mix-lean {lean} is refused with --captain-ruling {ruling}: only row {EMBRYO_R78_ROW} "
                     "orders it")
        return None
    name = _r78_file(a)
    if name is None:
        sys.exit(f"--captain-ruling {EMBRYO_R78_ROW} orders two files, each only as O's command with exactly its own "
                 f"option string, {EMBRYO_R78_COMMAND.format(option='<OPTION>')}, on the test board at the default "
                 "seed and nothing else changed: "
                 + "; ".join(f"{f} `{v[0]}`" for f, v in EMBRYO_R78_FILES.items()))
    option, expression, mix_lean, t_mix = EMBRYO_R78_FILES[name]
    if not _r78_row_names(name):
        sys.exit(f"row {EMBRYO_R78_ROW}'s file {name} is refused: experiments/harness/TRIES.md must hold exactly one "
                 f"row {EMBRYO_R78_ROW} whose first cell is the captain's order (\"ruling, no run\") naming `{option}` "
                 f"({name}) and `--captain-ruling {EMBRYO_R78_ROW}`, each as a whole backticked string")
    return {"file": name, "option": option, "expression": expression, "mix_lean": mix_lean, "t_mix": t_mix}


def _bayes_ordered() -> bool:
    """Row R71's first cell is the captain's order ("ruling, no run") naming BF on O and on SO only and
    `--captain-ruling R71`, in exactly one row."""
    return _captain_order_names(EMBRYO_BF_ROW, EMBRYO_BF_WORDS)


def _bf_file(a) -> str | None:
    """The file of row R71 whose command the arguments are exactly (the test board may be named), else None."""
    got = dict(vars(a))
    if got.get("board") == TEST_BOARDS["embryo"]:
        got["board"] = None
    for name, command in EMBRYO_BF_COMMANDS.items():
        if got == vars(build_parser().parse_args([*command.split(), "--out", str(a.out)])):
            return name
    return None


def _bayes_refusals(a) -> bool:
    """Row R71's option --shape-interp overlay_bayes: taken only as one of the two files' commands, with
    --captain-ruling R71, while row R71's first cell orders it, and written under that file's name; every refusal comes
    before anything is read. --captain-ruling R71 without the option is refused. True when the option is taken."""
    shape, ruling = getattr(a, "shape_interp", SHAPE_INTERP_OFF), getattr(a, "captain_ruling", None)
    if shape not in shape_interp.BAYES_FORMS:
        if ruling == EMBRYO_BF_ROW:
            sys.exit(f"embryo --captain-ruling {EMBRYO_BF_ROW} goes with --shape-interp overlay_bayes only (row "
                     f"{EMBRYO_BF_ROW}'s files: " + "; ".join(EMBRYO_BF_COMMANDS.values()) + ")")
        return False
    if ruling is None:
        sys.exit(f"--shape-interp overlay_bayes is the captain's order of row {EMBRYO_BF_ROW} alone and needs "
                 f"--captain-ruling {EMBRYO_BF_ROW}")
    if ruling != EMBRYO_BF_ROW:
        sys.exit(f"--shape-interp overlay_bayes is refused with --captain-ruling {ruling}: only row {EMBRYO_BF_ROW} "
                 "orders it")
    name = _bf_file(a)
    if name is None:
        sys.exit("--shape-interp overlay_bayes is refused except on O's or SO's command with it, "
                 + " or ".join(EMBRYO_BF_COMMANDS.values()) + f": the test board, the default seed and nothing else "
                 f"changed (row {EMBRYO_BF_ROW}: BF on O and on SO only)")
    if Path(a.out).name != name:
        sys.exit(f"this command writes row {EMBRYO_BF_ROW}'s file {name} only: --out must name it")
    if not _bayes_ordered():
        sys.exit(f"--shape-interp overlay_bayes is refused: experiments/harness/TRIES.md must hold exactly one row "
                 f"{EMBRYO_BF_ROW} whose first cell is the captain's order (\"ruling, no run\") naming "
                 f"{EMBRYO_BF_WORDS} and `--captain-ruling {EMBRYO_BF_ROW}`")
    return True


def _bf_prior():
    """Row R71's prior at test: the released embryo specimens of vec.shape_interp.BAYES_PRIOR_STAGES, each loaded as a
    prediction input is (synthetic records dropped), coordinates only."""
    coords = {}
    for stage in shape_interp.BAYES_PRIOR_STAGES:
        a = load_stage(DATA / "T2" / f"{stage}.h5ad", lineage=False, drop_duplicates=True)
        coords[stage] = np.asarray(a.obsm[COORD_KEY], dtype=float)[:, :3]
    return shape_interp.bayes_prior(coords)


def _run_sorted_refusals(a) -> str | None:
    """--neighbourhood run_sorted's refusals (row E57's file), before any board or stage is read. Returns the ruling
    row when the option is given, else None: without it every command is what it was. Refused on any command but O's
    with it (so with --expression other than mixture, --herd, --draw midpoint, --move-coupling aligned, --row-order
    copy_parity, --copy-rows and its flags, another seed or a dry run), without --captain-ruling, and unless
    experiments/harness/TRIES.md holds exactly one row of the given ID whose first cell is the captain's ruling naming
    `--neighbourhood run_sorted` and `--captain-ruling <ROW>` (_captain_ruling_row)."""
    mode = getattr(a, "neighbourhood", NEIGHBOURHOODS[0])
    if mode == NEIGHBOURHOODS[0]:
        return None
    if mode not in NEIGHBOURHOODS:  # argparse's choices hold this already; the placebo is never built
        sys.exit(f"--neighbourhood {mode} is never built (harness only)")
    o_command = (a.method == "fgw_shape0" and a.positions == t2_build.POSITIONS[1] and a.expression == "mixture"
                 and a.shape_interp == "overlay" and a.shape_interp_test_ok and not a.herd
                 and a.draw == t2_build.DRAWS[0] and a.move_coupling == t2_build.MOVE_COUPLINGS[0]
                 and a.row_order == t2_build.ROW_ORDERS[0] and not a.row_order_test_ok
                 and not getattr(a, "copy_rows", False) and not getattr(a, "copy_rows_test_ok", False)
                 and not getattr(a, "copy_rows_size", False) and not getattr(a, "copy_rows_size_test_ok", False)
                 and a.seed == harness.BUILD_SEED_OFFSET + 1 and a.board in (None, TEST_BOARDS["embryo"]))
    if not o_command:
        sys.exit(f"--neighbourhood run_sorted is refused except on O's command with it, "
                 f"{EMBRYO_RS_COMMAND.format(row='<ROW>')}: the test board, the default seed and nothing else changed "
                 "(row E57's file)")
    ruling = getattr(a, "captain_ruling", None)
    if ruling is None:
        sys.exit("--neighbourhood run_sorted is the captain's ruling alone: it needs --captain-ruling <ROW>, the row "
                 f"of experiments/harness/TRIES.md whose first cell is his ruling naming {EMBRYO_RS_OPTION} and "
                 "`--captain-ruling <ROW>`")
    if not _captain_ruling_row(ruling, EMBRYO_RS_OPTION):
        sys.exit(f"--neighbourhood run_sorted is refused: experiments/harness/TRIES.md has no single row {ruling} "
                 f"whose first cell is the captain's ruling (\"ruling, no run\") naming {EMBRYO_RS_OPTION} and "
                 f"`--captain-ruling {ruling}`")
    return ruling


def _copy_rows_refusals(a) -> None:
    """--copy-rows' refusals, before any board or stage is read (off by default: nothing here acts without it)."""
    copy_rows, ok = getattr(a, "copy_rows", False), getattr(a, "copy_rows_test_ok", False)
    sized, sized_ok = getattr(a, "copy_rows_size", False), getattr(a, "copy_rows_size_test_ok", False)
    if sized_ok and not sized:  # row E52's option (off by default)
        raise SystemExit("--copy-rows-size-test-ok is for --copy-rows-size only")
    if sized and not copy_rows:
        raise SystemExit("--copy-rows-size is for --copy-rows only")
    if sized and a.expression != COPY_ROWS_SIZE_EXPRESSION:
        raise SystemExit("--copy-rows-size is for --expression mixture only")
    if ok and not copy_rows:
        raise SystemExit("--copy-rows-test-ok is for --copy-rows only")
    if not copy_rows:
        return
    if a.method != "fgw_shape0" or a.positions != t2_build.POSITIONS[1]:
        raise SystemExit("--copy-rows is for --method fgw_shape0 --positions floor_exact only")
    if (a.row_order != t2_build.ROW_ORDERS[0] or a.shape_interp != SHAPE_INTERP_OFF or a.herd
            or a.draw != t2_build.DRAWS[0] or a.move_coupling != t2_build.MOVE_COUPLINGS[0]
            or a.expression not in COPY_ROWS_EXPRESSIONS):
        raise SystemExit("--copy-rows is refused with --row-order copy_parity, --shape-interp, --herd, --draw "
                         "midpoint, --move-coupling aligned or an --expression other than mixture or series")


def _same_bits(C, ref32: np.ndarray) -> bool:
    """Whether coordinates are ref32 exactly: shape, row order and float32 bits."""
    C32 = np.ascontiguousarray(np.asarray(C)[:, :3], dtype=np.float32)
    R = np.ascontiguousarray(ref32, dtype=np.float32)
    return C32.shape == R.shape and np.array_equal(C32.view(np.uint32), R.view(np.uint32))


def _copy_rows_size_target() -> tuple[float, str, str]:
    """--copy-rows-size: the embryo size JSON's pre-registered full-cloud RMS, read by code with build_embryo's checks
    and message (valid, raw mode), and the file's path and sha256. Refusals show no value."""
    path = SIZES["embryo"]
    raw = Path(path).read_bytes()
    size = json.loads(raw)
    if not size.get("valid") or size["preregistered"]["z_mode"] != "raw":
        raise SystemExit("the embryo size JSON is not a valid raw-mode result")
    target = float(size["preregistered"]["pred_rms_full"])  # exactly as build_embryo reads it for the positions file
    if not np.isfinite(target) or target <= 0:
        raise SystemExit("the embryo size JSON's target is not finite and positive (values not shown)")
    try:
        shown = str(Path(path).resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        shown = str(path)
    return target, shown, hashlib.sha256(raw).hexdigest()


def _copy_rows_rms_ok(S32: np.ndarray, target: float) -> bool:
    """The RMS self-check: the stored float32 coordinates' RMS radius within COPY_ROWS_SIZE_RMS_TOL of the target."""
    return bool(abs(t2_build.rms_radius(np.asarray(S32, dtype=np.float64)) / target - 1.0) <= COPY_ROWS_SIZE_RMS_TOL)


def _copy_rows_parity_ok(S32: np.ndarray, floor32: np.ndarray) -> bool:
    """The prefix parity self-check, count-free: for every count n from COPY_ROWS_PARITY_FROM to the copy's own, the
    SVD parity of the rows the scorer's _match_n takes at metric seed COPY_ROWS_PARITY_SEED from the stored
    coordinates equals the floor copy's (the rows depend on the counts only)."""
    from vec.evaluate import match_rows, svd_parity

    S, F = np.asarray(S32, dtype=np.float32), np.asarray(floor32, dtype=np.float32)
    if S.shape != F.shape:
        return False
    for n in range(COPY_ROWS_PARITY_FROM, len(S) + 1):
        zeros = np.zeros((n, 3))
        if (svd_parity(match_rows(S, zeros, COPY_ROWS_PARITY_SEED)[0])
                != svd_parity(match_rows(F, zeros, COPY_ROWS_PARITY_SEED)[0])):
            return False
    return True


def _build_copy_rows(a, board, dry: bool, genes, files: dict, inputs: dict, count: dict) -> dict:
    """--copy-rows (row E50 as specified; the captain's order of 2026-11-21): the floor copy's rows, order, count and
    float32 coordinates (floor_copy_coords of the earlier bracket as released), the expression of the handed-over
    form's rule (mixture or series). No size rule, centring, jitter or row-order step. The limits are checked here and
    every refusal on this path shows no value; the coordinates are checked bit for bit before writing and on the file
    read back."""
    import anndata as ad

    stage_file = files["before"][0]
    before, after = inputs["before"], inputs["after"]
    released = read_stage(DATA / stage_file, list(genes))
    floor = floor_copy_coords(stage_file)
    Eb, Ea = dense(before[:, genes].X).astype(np.float32), dense(after[:, genes].X).astype(np.float32)
    try:
        dropped = before.uns.get("qc", {}).get("dropped_rows", [])
        res = t2_build.copy_rows_cloud(np.asarray(released.obsm[COORD_KEY])[:, :3], dropped, Eb, Ea,
                                       dense(released.X).astype(np.float32), 0.5, t2_build.COPY_PARITY_COPY_SEED,
                                       a.seed)
    except ValueError as e:
        raise SystemExit(f"--copy-rows refused: {e}") from None
    if not _same_bits(res.coords, floor):
        raise SystemExit("--copy-rows: the coordinates are not the floor copy's, bit for bit (values not shown)")
    sized = getattr(a, "copy_rows_size", False)  # row E52 (off by default): the rescale, then only writing
    coords, written = res.coords, floor
    if sized:
        target, size_file, size_sha = _copy_rows_size_target()
        try:
            coords = t2_build.copy_rows_sized(res.coords, target)
        except ValueError:
            raise SystemExit("--copy-rows-size refused: the copy cannot be rescaled (values not shown)") from None
        written = coords.astype(np.float32)  # what the file stores: round to nearest
        if not _copy_rows_rms_ok(written, target):
            raise SystemExit("--copy-rows-size: the stored coordinates' RMS is not the target's (values not shown)")
        if not _copy_rows_parity_ok(written, floor):
            raise SystemExit("--copy-rows-size: a prefix parity differs from the floor copy's (values not shown)")
    expr_rec = None
    if a.expression == "mixture":
        X = np.empty((len(res.rows), len(genes)), dtype=np.float32)
        X[res.stage == 0], X[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
    else:  # series, the handed-over exact-position file's own rule
        X, expr_rec = t2_build.series_rendered_expression(before, after, res, genes, 0.5, [a.seed, 7])
    if len(X) < board.min_cells or (board.max_cells is not None and len(X) > board.max_cells):
        raise SystemExit("the copy-rows count is outside the board's limits (values not shown)")
    try:
        _, report = t2_build.build_submission(board, X, coords, genes=list(genes), out=a.out,
                                              allow_duplicate_coords=True, allow_below_truth=True)
    except ValueError:
        raise SystemExit("build_submission refused the copy-rows file (values not shown)") from None
    if not report.ok:
        raise SystemExit("the built copy-rows file fails its check (values not shown)")
    back = ad.read_h5ad(a.out)
    if not _same_bits(back.obsm[COORD_KEY], written):
        Path(a.out).unlink(missing_ok=True)
        provenance.stamp_path(Path(a.out)).unlink(missing_ok=True)
        if sized:
            raise SystemExit("--copy-rows-size: the written file's coordinates are not the rescaled copy's, bit for "
                             "bit; the file was removed (values not shown)")
        raise SystemExit("--copy-rows: the written file's coordinates are not the floor copy's, bit for bit; the file "
                         "was removed (values not shown)")
    rec = dict(res.record)
    if expr_rec is not None:
        rec.update({f"series_{k}": v for k, v in expr_rec.items() if not isinstance(v, (dict, list))})
    rec["expression"] = a.expression
    out = {"entry": a.method, "board": board.key, "dry_run": dry, "n_cells": int(len(X)), "seed": a.seed,
           "count_decision": count, "target_rms": None,
           "record": {k: v for k, v in rec.items() if not isinstance(v, (dict, list))},
           "copy_rows": True,
           "copy_rows_copy": f"vec.evaluate.subsample({stage_file} as released, {SUBSAMPLE_FRAC}, "
                             f"{t2_build.COPY_PARITY_COPY_SEED})",
           "count_rule": COPY_ROWS_COUNT_RULE, "size_rule": "not applied", "allow_below_truth": True,
           COPY_ROWS_N: int(len(floor))}
    if sized:  # row E52 only: the record says the size step and where its target came from
        out.update({"copy_rows_size": True, "size_rule": COPY_ROWS_SIZE_RULE, "target_rms": target,
                    "size_file": size_file, "size_sha256": size_sha})
    return out


def build_embryo(a) -> dict:
    r78 = _r78_refusals(a)  # row R78's check first: when it returns a file, R75's, R59's and series_back's are skipped
    r75 = (_r75_refusals(a)  # row R75's five files (off by default): first, before the older refusals read the ruling
           if r78 is None else None)
    mix_lean = _mix_lean_refusals(a) if (r75 is None  # row R59's option (off by default): before any read
                                         and r78 is None) else False
    # row R48's option (off by default): before anything is read; with --mix-lean the command is O's, mixture included
    bayes = _bayes_refusals(a)  # row R71's option, before the series_back check refuses its ruling
    series_back = _series_back_refusals(a) if not (mix_lean or bayes or r75 or r78) else False
    rs_ruling = _run_sorted_refusals(a)  # row E57's file (off by default): before anything is read
    _copy_rows_refusals(a)
    if a.positions != t2_build.POSITIONS[0] and (a.method != "fgw_shape0" or a.herd
                                                 or a.draw != t2_build.DRAWS[0]):
        raise SystemExit(f"--positions {a.positions} is for --method fgw_shape0, without --herd and --draw midpoint")
    if a.row_order_test_ok and a.row_order != t2_build.ROW_ORDERS[1]:
        raise SystemExit("--row-order-test-ok is for --row-order copy_parity only")
    shape = None if a.shape_interp == SHAPE_INTERP_OFF else a.shape_interp
    if shape is not None and shape not in BUILDABLE_SHAPES:  # by name: no control
        raise SystemExit(f"--shape-interp {shape}: a harness-only control or an unknown form is never built")
    if a.shape_interp_test_ok and shape is None:
        raise SystemExit("--shape-interp-test-ok is for --shape-interp with a form only")
    if shape is not None and (a.method != "fgw_shape0" or a.positions != t2_build.POSITIONS[1]):
        raise SystemExit("--shape-interp is for --method fgw_shape0 --positions floor_exact only")
    board, dry = _board("embryo", a.board)
    if (a.row_order != t2_build.ROW_ORDERS[0] and (not dry or str(board.split).startswith("test"))
            and not a.row_order_test_ok):  # before anything loads: nothing about the test count is computed
        raise SystemExit("--row-order copy_parity on a test board needs the captain's yes (--row-order-test-ok); a "
                         "refusal there shows one bit of the test count")
    if (shape is not None and (not dry or str(board.split).startswith("test"))
            and not a.shape_interp_test_ok):  # before anything loads
        raise SystemExit("--shape-interp on a test board needs the captain's yes (--shape-interp-test-ok)")
    if (getattr(a, "copy_rows", False) and (not dry or str(board.split).startswith("test"))
            and not a.copy_rows_test_ok):  # before anything loads: a refusal or a pass there shows one bit
        raise SystemExit("--copy-rows on a test board needs the captain's yes (--copy-rows-test-ok); a refusal or a "
                         "pass there shows one bit about the board's limits")
    if (getattr(a, "copy_rows_size", False) and (not dry or str(board.split).startswith("test"))
            and not a.copy_rows_size_test_ok):  # row E52's option: before anything loads
        raise SystemExit("--copy-rows-size on a test board needs the captain's yes (--copy-rows-size-test-ok)")
    genes = board.genes
    files = {"before": ("T2/E7.5.h5ad", 7.5), "after": ("T2/E8.0.h5ad", 8.0)}
    if a.method == "pseudobulk_shift":  # only it uses the stage before the bracket (E7.25 lacks two panel genes)
        files["earlier"] = ("T2/E7.25.h5ad", 7.25)
    inputs = {r: load_stage(DATA / f, genes=genes, drop_duplicates=True) for r, (f, _) in files.items()}
    count = _t2_count(board, inputs["before"].n_obs)  # provenance; under --copy-rows it does not set the count
    if getattr(a, "copy_rows", False):  # off by default
        return _build_copy_rows(a, board, dry, genes, files, inputs, count)
    view = _view("test:T2e", "T2", "interp", "embryo", board.key, files, "before", 7.75)
    ctx = harness.BuildContext(task=view, board=harness._candidate_board(board), inputs=inputs, seed=a.seed)
    size = json.loads(SIZES["embryo"].read_text())
    if not size.get("valid") or size["preregistered"]["z_mode"] != "raw":
        raise SystemExit("the embryo size JSON is not a valid raw-mode result")
    target = float(size["preregistered"]["pred_rms_full"])  # what apply_size_procedure applies in raw mode
    if a.method == "pseudobulk_shift":
        stage = harness._pseudobulk_shift(ctx)
        coords, rec = t2_build.apply_size_procedure(stage.coords, SIZES["embryo"], expect_board="embryo")
        X = stage.X
    else:
        before, after = inputs["before"], inputs["after"]
        bf_prior = _bf_prior() if bayes else None  # row R71: the released specimens, after every refusal
        Eb, Ea = dense(before[:, genes].X).astype(np.float32), dense(after[:, genes].X).astype(np.float32)
        res = t2_build.interpolate_clouds(np.asarray(before.obsm[COORD_KEY], float)[:, :3],
                                          np.asarray(after.obsm[COORD_KEY], float)[:, :3], Eb, Ea, 0.5,
                                          harness.count_for(board, before.n_obs), a.seed, target_rms=target,
                                          t_pos=EMBRYO_T_POS[a.method],
                                          with_partners=a.expression in ("blend", "latent"), herd=a.herd,
                                          move_coupling=a.move_coupling, draw=a.draw,
                                          log_depth=(None if a.draw == t2_build.DRAWS[0]
                                                     else (t2_build.log_panel_counts(before, genes),
                                                           t2_build.log_panel_counts(after, genes))),
                                          with_coupling=a.expression == "latent",
                                          **({} if shape is None else {"shape_form": shape}),
                                          **({} if bf_prior is None else {"shape_prior": bf_prior}),
                                          **({"t_mix": EMBRYO_OT23_T_MIX} if mix_lean else {}),  # row R59 only
                                          **({} if r75 is None or r75["t_mix"] is None  # row R75's leans only
                                             else {"t_mix": float(r75["t_mix"])}),
                                          **({} if r78 is None else {"t_mix": float(r78["t_mix"])}),
                                          positions=a.positions)
        if (shape in shape_interp.COMPOSED_FORMS and (not dry or str(board.split).startswith("test"))
                and res.record.get("shape_interp_blocks") != OSR_TEST_BLOCKS):  # row E55, cell 3: "must be"
            raise SystemExit(f"--shape-interp {shape} on a test board needs the brackets' blocks {OSR_TEST_BLOCKS}; "
                             f"they are {res.record.get('shape_interp_blocks')}: refused, nothing written")
        expr_rec = None
        if a.expression == "mixture":
            X = np.empty((len(res.rows), len(genes)), dtype=np.float32)
            X[res.stage == 0], X[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
        elif a.expression == "depth":  # the harness's fgw_depth (try E7, experiments/harness/TRIES.md)
            X, expr_rec = t2_build.depth_matched_expression(before, after, res, genes, 0.5, [a.seed, 7])
        elif a.expression == "latent":  # the harness's fgw_latent (try E14): torch and scvi-tools load only here
            from vec import t2_latent

            X, expr_rec = t2_latent.latent_midpoint_counts(t2_latent.task_model(before, after, genes), before, after,
                                                           res, genes, 0.5, [a.seed, t2_latent.EXPRESSION_STREAM])
        elif a.expression == "series":  # the harness's series rendering (tries E24a-E24c, E25a-E25b)
            X, expr_rec = t2_build.series_rendered_expression(before, after, res, genes, 0.5, [a.seed, 7])
        elif a.expression == "render_hi":  # the harness's depth top-up to the deeper bracket (tries E19c and E24b)
            X, expr_rec = t2_build.depth_rendered_expression(before, after, res, genes, 0.5, [a.seed, 7], mode="hi")
        elif a.expression == "series_back":  # row R48's file: the E8.0 rows in E7.5's run profile and depth
            X, expr_rec = t2_build.series_back_expression(before, after, res, genes, 0.5, [a.seed, 7])
        else:  # the harness's fgw_blend (try E8)
            X, expr_rec = t2_build.blend_expression(before, after, res, genes, 0.5, [a.seed, 7])
        guard = t2_build.cloud_guard(res.coords, res.pool_coords, seed=a.seed, source_weights=res.pool_weights)
        if not guard.ok:
            raise SystemExit(f"cloud guard flags {guard.flags}")
        coords, rec = res.coords, res.record
        if expr_rec is not None:
            prefix = a.expression if a.expression in ("latent", "series", "render_hi", "series_back") else "depth"
            rec = {**rec, "expression": a.expression,
                   **{f"{prefix}_{k}": v for k, v in expr_rec.items() if not isinstance(v, (dict, list))
                      and (prefix, k) != ("latent", "retried")}}  # a test build's record holds no cell number
        if rs_ruling is not None:  # row E57's file: the expression rows re-assigned on the final coordinates
            step = harness.rs_step(before, after, res, a.neighbourhood)
            X = np.asarray(X)[step.src]
    if a.row_order != t2_build.ROW_ORDERS[0]:  # off by default: the rows re-ordered, nothing else changed
        if not board.truth_cells:
            raise SystemExit("--row-order copy_parity needs the board's truth_cells (read by code)")
        try:
            order, order_rec = t2_build.copy_parity_order(coords, floor_copy_coords(files["before"][0]),
                                                          board.truth_cells, a.seed)
        except ValueError as e:
            raise SystemExit(f"--row-order copy_parity refused: {e}") from None
        X, coords = np.asarray(X)[order], np.asarray(coords)[order]
        rec = {**rec, **order_rec, "row_order_copy": f"vec.evaluate.subsample({files['before'][0]} as released, "
               f"{SUBSAMPLE_FRAC}, {t2_build.COPY_PARITY_COPY_SEED})"}
    _, report = t2_build.build_submission(board, X, coords, genes=list(genes), out=a.out)
    if not report.ok:
        raise SystemExit(f"the built file fails its check: {report.errors}")
    out = {"entry": a.method, "board": board.key, "dry_run": dry, "n_cells": int(len(X)), "seed": a.seed,
           "count_decision": count, "target_rms": target,
           "record": {k: v for k, v in rec.items() if not isinstance(v, (dict, list))}}
    if shape is not None:  # written only when the option is given: flat values, no cell number
        out["shape_interp"] = {"form": shape, "module": "vec.shape_interp", "ruling": None,
                               "row": ("E55 (as specified)" if shape in shape_interp.COMPOSED_FORMS
                                       else "E49 (as drafted)"),
                               "order": SHAPE_INTERP_ORDER,
                               **{k.removeprefix("shape_interp_"): v for k, v in rec.items()
                                  if k.startswith("shape_interp_") and not isinstance(v, (dict, list))}}
    if series_back:  # row R48's file only: what it is and who ordered it, flat values, no cell number
        out["series_back"] = {"option": "--expression series_back", "file": EMBRYO_OVSB_FILE,
                              "captain_ruling": EMBRYO_OVSB_ROW, "order": EMBRYO_OVSB_ORDER,
                              "command": EMBRYO_OVSB_COMMAND, "o_command": EMBRYO_O_COMMAND,
                              "rule": t2_build.SERIES_BACK_RULE, "step0": EMBRYO_OVSB_STEP0}
    if mix_lean:  # row R59's file only: what it is and who ordered it, flat values, no cell number
        out["mix_lean"] = {"option": "--mix-lean later", "file": EMBRYO_OT23_FILE, "captain_ruling": EMBRYO_OT23_ROW,
                           "order": EMBRYO_OT23_ORDER, "command": EMBRYO_OT23_COMMAND, "o_command": EMBRYO_O_COMMAND,
                           "t_mix": EMBRYO_OT23_T_MIX, "t_mix_exact": "2/3", "rule": EMBRYO_OT23_RULE,
                           "gate": EMBRYO_OT23_GATE}
    if r75 is not None:  # row R75's file: what it is and who ordered it, flat values, no cell number
        f75, opt75, t75 = r75["file"], r75["option"], r75["t_mix"]
        out["r75"] = {"option": opt75, "file": f75, "captain_ruling": EMBRYO_R75_ROW,
                      "order": EMBRYO_R75_ORDER.format(option=opt75, file=f75, what=EMBRYO_R75_WHAT[f75]),
                      "command": EMBRYO_R75_COMMAND.format(option=opt75), "o_command": EMBRYO_O_COMMAND,
                      "expression": r75["expression"], "mix_lean": r75["mix_lean"],
                      "t_mix": None if t75 is None else float(t75),
                      "t_mix_exact": None if t75 is None else f"{t75.numerator}/{t75.denominator}",
                      "rule": EMBRYO_R75_RULES[f75]}
    if r78 is not None:  # row R78's file: the captain's order, the command and its exact t_mix
        f78, opt78, t78 = r78["file"], r78["option"], r78["t_mix"]
        out["r78"] = {"option": opt78, "file": f78, "captain_ruling": EMBRYO_R78_ROW,
                      "order": EMBRYO_R78_ORDER.format(option=opt78, file=f78, words=EMBRYO_R78_WORDS[f78],
                                                       what=EMBRYO_R78_WHAT[f78]),
                      "command": EMBRYO_R78_COMMAND.format(option=opt78), "o_command": EMBRYO_O_COMMAND,
                      "expression": r78["expression"], "mix_lean": r78["mix_lean"], "t_mix": float(t78),
                      "t_mix_exact": f"{t78.numerator}/{t78.denominator}", "rule": EMBRYO_R78_RULES[f78]}
    if rs_ruling is not None:  # row E57's file only: what it is, who ordered it and the plane; flat, no cell number
        out["neighbourhood"] = {"option": "--neighbourhood run_sorted", "file": EMBRYO_RS_FILE,
                                "captain_ruling": rs_ruling, "command": EMBRYO_RS_COMMAND.format(row=rs_ruling),
                                "o_command": EMBRYO_O_COMMAND, "rule": run_sort.RULE, "read": EMBRYO_RS_READ,
                                **{k: step.record[k] for k in ("normal_x", "normal_y", "normal_z",
                                                               "contrast_earlier")}}
    if bayes:  # row R71's file: what it is and where its order and gate stand
        name = _bf_file(a)
        out["shape_interp"]["ruling"] = EMBRYO_BF_ROW
        out["shape_interp"]["row"] = f"{EMBRYO_BF_ROW} (the captain's order); gate row {EMBRYO_BF_GATE_ROW}"
        out["bayes"] = {"option": "--shape-interp overlay_bayes", "file": name, "base": EMBRYO_BF_BASES[name],
                        "captain_ruling": EMBRYO_BF_ROW, "order": EMBRYO_BF_ORDER, "command": EMBRYO_BF_COMMANDS[name],
                        "gate": EMBRYO_BF_GATE, "rule": EMBRYO_BF_RULE,
                        "prior": "|".join(shape_interp.BAYES_PRIOR_STAGES)}
    return out


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="which", required=True)
    for name in ("t1", "heart", "embryo"):
        p = sub.add_parser(name)
        p.add_argument("--out", required=True)
        p.add_argument("--board", help="a validation board for a dry run (default: the test board)")
        p.add_argument("--seed", type=int, default=harness.BUILD_SEED_OFFSET + 1)
        if name == "t1":
            p.add_argument("--cm-share", type=float, help="cardiomyocyte share of heart cells (default: E10.5's own)")
            p.add_argument("--blood-share", type=float, help="blood share of chamber cells (default: E10.5's own)")
            p.add_argument("--herd", action="store_true",
                           help="herd each category's draw towards its pool's mean expression (T1 try T5 of "
                                "experiments/harness/TRIES.md; off: the entry as adopted so far)")
            p.add_argument("--hold-unbracketed", action=argparse.BooleanOptionalAction, default=True,
                           help="R1 (experiments/harness/TRIES.md; the captain's ruling of 2026-11-05): cells without "
                                "a Feng group keep their E10.5 expression (default); --no-hold-unbracketed gives "
                                "them B0's damped E9.5 -> E10.5 momentum, M1 as built before R1")
            p.add_argument("--feng-hygiene", action=argparse.BooleanOptionalAction, default=None,
                           help="R2 (experiments/harness/TRIES.md): also hold R2's gene modules flat in Feng's "
                                "direction (the union N12 wrote, pinned in configs/t1_hygiene.yaml; refused with "
                                "--no-hold-unbracketed or while Blood has a Feng direction). Not given: on with "
                                "--hold-unbracketed (the default), off with --no-hold-unbracketed")
            p.add_argument("--strong-null", action=argparse.BooleanOptionalAction, default=False,
                           help="T12a, T12b (experiments/harness/TRIES.md; the captain's rulings of 2026-11-10 (a) and "
                                "(d)): add the strong-null down-tail prior after M1's steps (off by default; refused "
                                "until configs/t1_strong_null.yaml says approved, and without R1 and R2)")
            p.add_argument("--base", choices=T1_BASES, default=T1_BASES[0],
                           help="the pool the cells are drawn from: pooled = every E10.5 cell (the default); lib3 = "
                                "E10.5's dissected-heart library alone at its own chamber mix (row R6 of "
                                "experiments/harness/TRIES.md; T1's first upload); mid = every library-3 cell and, "
                                "per chamber category, a seeded fraction of the other libraries' cells that puts "
                                "library 3's share half-way to 1 (row R7; T1's second upload, with --cm-share). The "
                                "count, Feng's direction, R1, R2 and g_mom's fit are the same in all three. With "
                                f"--feng-groups recomputed, pooled only with --captain-ruling {T1_R37_ROW} (row "
                                f"{T1_R37_ROW}'s one file: {T1_POOLED_COMMAND})")
            p.add_argument("--lib3-pair", action="store_true",
                           help="with --base lib3 and --cm-share: build the second file of the pair try T32 priced "
                                "(experiments/harness/TRIES.md), library 3's cells with the cardiomyocyte share "
                                "moved. Off by default and no entry: T1's second upload stays on the pooled base "
                                "unless the captain rules")
            p.add_argument("--g-mom", action="store_true",
                           help="add g_mom's within-type factors at test (T1 try T11 of experiments/harness/TRIES.md; "
                                "fitted on E9.5 -> E10.5, applied last, after Feng's direction, on cardiomyocytes and "
                                "endocardium only (R3b), R1's held labels left alone, each cell's expected depth "
                                "kept, no share, cell count or scope changed; refused with --no-hold-unbracketed; "
                                "off: the entry as built, "
                                "and only the captain's ruling turns it on for an upload)")
            p.add_argument("--feng-groups", choices=FENG_GROUPS, default=FENG_GROUPS_DEFAULT,
                           help="the cell groups Feng's direction is formed on: authors = Feng's own labels (the "
                                "default, every file as built); none = one file, the first upload's command "
                                f"({FENG_NONE_COMMAND}) with Feng's direction taken out and nothing put in its place: "
                                "every cell keeps its E10.5 expression before g_mom's edit and no Feng file, R2 list "
                                "or other outside data is read; g_mom's factors edit the cm and endocardium cells "
                                "whose label the readout config's regexes give a group (all but Early VIC, as in "
                                "the entry; the list of groups with a Feng direction is not read). Refused "
                                "without --g-mom and with --base pooled or mid, --cm-share, --blood-share, "
                                "--lib3-pair, --herd, --no-hold-unbracketed, --feng-hygiene or --strong-null. The "
                                "build record then says the direction step was removed. recomputed = the same one "
                                "command on the groups experiments/external/feng_relabel.py recomputes from the "
                                "permitted cells alone, refused unless its committed record routes to recomputed; "
                                f"with --base pooled only as row {T1_R37_ROW}'s one file ({T1_POOLED_COMMAND}). "
                                "Only the captain's ruling puts either in an upload")
            p.add_argument("--feng-none-ruling", default=None,
                           help="with --feng-groups none only, and only after a permitted-only relabelling whose "
                                "record routes to recomputed: a later docs/PLAN.md decision-log entry of the captain "
                                "(naming N54 and containing \"sets the N54 route aside\"; not the run's own "
                                "ruling; added after the record's commit) that sets the route aside")
            p.add_argument("--g-mom-damping", type=float, default=None,
                           help="candidate A of docs/PLAN.md decision log 2026-11-20 (i) (off by default, no ruling): "
                                f"g_mom's factors finished at 1.0 in place of {harness.B0_DAMPING:g} (practice arm "
                                "gmom_s1, rows T49 and T61). Takes 0.5 (the default, today's file) or 1.0; only with "
                                f"{T1_CANDIDATE_COMMAND} at the default seed on the test board. "
                                f"{G_MOM_DAMPINGS_R22[0]:g} and {G_MOM_DAMPINGS_R22[1]:g} are row {T1_R22_ROW}'s "
                                "second T1 search line (the captain's order of 2026-11-22) and need --captain-ruling "
                                f"{T1_R22_ROW}, with no other change; with --feng-fraction-mult "
                                f"{FENG_FRACTION_MULT_LINE5:g}, {G_MOM_DAMPINGS_LINE5[0]:g} and "
                                f"{G_MOM_DAMPINGS_LINE5[1]:g} are row {T1_R22_ROW}'s fifth search line and need "
                                f"--captain-ruling {T1_R22_ROW}; with --feng-fraction-mult "
                                f"{FENG_FRACTION_MULT_LINE6:g}, {G_MOM_DAMPINGS_LINE6[0]:g} and "
                                f"{G_MOM_DAMPINGS_LINE6[1]:g} are row {T1_R22_ROW}'s sixth search line (T1 line 6) "
                                f"and need --captain-ruling {T1_R22_ROW}; with --feng-fraction-mult "
                                f"{FENG_FRACTION_MULT_LINE7:g}, {G_MOM_DAMPING_LINE7:g} is row {T1_R22_ROW}'s seventh "
                                f"search line (T1 line 7) and needs --captain-ruling {T1_R22_ROW}")
            p.add_argument("--feng-fraction-mult", type=float, default=None,
                           help="candidate B of docs/PLAN.md decision log 2026-11-20 (i) (off by default, no ruling): "
                                "the outside atlas's direction at this multiple of the calendar fraction (2/5). Takes "
                                "1.0 (the default, today's file) or 1.5; only with the command "
                                f"{T1_CANDIDATE_COMMAND} at the default seed on the test board. "
                                f"{FENG_FRACTION_MULT_R20:g} is row {T1_R20_ROW}'s file (the captain's order of "
                                f"2026-11-21) and needs --captain-ruling {T1_R20_ROW}, with no other change; "
                                f"{FENG_FRACTION_MULTS_R22[0]:g} and {FENG_FRACTION_MULTS_R22[1]:g} are row "
                                f"{T1_R22_ROW}'s first search line (the captain's order of 2026-11-22) and need "
                                f"--captain-ruling {T1_R22_ROW}, with no other change; with --g-mom-damping 1.0, "
                                f"{FENG_FRACTION_MULTS_LINE3[0]:g} and {FENG_FRACTION_MULTS_LINE3[1]:g} are row "
                                f"{T1_R22_ROW}'s third search line and need --captain-ruling {T1_R22_ROW}; with "
                                f"--g-mom-damping {G_MOM_DAMPINGS_LINE5[0]:g} or {G_MOM_DAMPINGS_LINE5[1]:g}, "
                                f"{FENG_FRACTION_MULT_LINE5:g} is row {T1_R22_ROW}'s fifth search line and needs "
                                f"--captain-ruling {T1_R22_ROW}; with --g-mom-damping {G_MOM_DAMPINGS_LINE6[0]:g} or "
                                f"{G_MOM_DAMPINGS_LINE6[1]:g}, {FENG_FRACTION_MULT_LINE6:g} is row {T1_R22_ROW}'s "
                                f"sixth search line (T1 line 6); with --g-mom-damping {G_MOM_DAMPING_LINE7:g}, "
                                f"{FENG_FRACTION_MULT_LINE7:g} is row {T1_R22_ROW}'s seventh search line (T1 line 7); "
                                f"each needs --captain-ruling {T1_R22_ROW}")
            p.add_argument("--released-bet", choices=tuple(RELEASED_BETS), default=None,
                           help="try T61 of experiments/harness/TRIES.md (off by default): a variant of T1's last file "
                                "(--feng-groups none only) that changes which cells are drawn: share_trend (the "
                                "chamber mix moved along the E9.5 -> E10.5-portions trend) or advance (cells drawn "
                                "inside each category towards the released momentum axis). Composition, a lever "
                                "closed by the post-exposure block: refused unless --captain-ruling names the "
                                "captain's ruling row for the arm")
            p.add_argument("--direction-momentum", choices=DIRECTION_MOMENTUM_ARMS, default=None,
                           help=f"candidate MS (row {T1_R23_ROW} of experiments/harness/TRIES.md, the captain's order "
                                "of 2026-11-23; off by default): on the command "
                                f"{T1_CANDIDATE_COMMAND} at the default seed on the test board, the drawn cells whose "
                                "label's group has a recomputed direction also get B0's per-category damped mean "
                                "momentum (fitted from E9.5 to E10.5's libraries 1 and 2), before the atlas step and "
                                f"g_mom. Needs --captain-ruling {T1_R23_ROW}; refused with --g-mom-damping, "
                                "--feng-fraction-mult, --released-bet, --bet-strength or --no-hold-unbracketed. "
                                f"With --feng-fraction-mult {T1_MSB_MULT:g} and --captain-ruling {T1_R26_ROW} only: "
                                f"row {T1_R26_ROW}'s combined file, the same step on candidate B's command. With "
                                f"--feng-fraction-mult {T1_BANKED_MULT:g} and --captain-ruling <ROW> only, ROW a row "
                                "whose first cell is the captain's ruling naming "
                                f"{T1_MS25_OPTION} and `--captain-ruling <ROW>`: file {T1_MS25_FILE}, the same step on "
                                f"the banked file {T1_BANKED_FILE}'s command. With --captain-ruling {T1_R51_ROW}: row "
                                f"{T1_R51_ROW}'s seven files on {T1_R51_BASE_FILE}'s command (see --captain-ruling). "
                                f"With --captain-ruling {T1_R60_ROW} only: row {T1_R60_ROW}'s seven files (T1's line "
                                f"12), each one change to {T1_R60_BASE_FILE}'s command, exactly one of these sets: "
                                + "; ".join(f"{f} `{v[0]}`" for f, v in T1_R60_FILES.items())
                                + " (--direction-momentum-damping takes 1.0 and 2.0 there), while row "
                                f"{T1_R60_ROW}'s first cell is the captain's order naming the set's option string and "
                                f"`--captain-ruling {T1_R60_ROW}`. With --captain-ruling {T1_R46_ROW} only: row "
                                f"{T1_R46_ROW}'s pair {T1_R46_FILE} `{T1_R46_OPTION}` (Blood's step on "
                                f"{T1_R46_BASE_FILE}'s command), while row {T1_R46_ROW}'s first cell is the captain's "
                                f"and names {T1_R46_FILE}. With --captain-ruling {T1_R66_ROW} only: row "
                                f"{T1_R66_ROW}'s thirty files (T1's third-round step-on slots: for K1, K2 and K3, "
                                "two or three of line 12's changes together; for C1, g_mom's damping 2.5; for C3, "
                                f"Blood's step at 2.0; each on {T1_R66_BASE_FILE}'s command), exactly one of these "
                                "sets: " + "; ".join(f"{f} `{v[0]}`" for f, v in T1_R66_FILES.items())
                                + f", while row {T1_R66_ROW}'s first cell is the captain's order naming the set's "
                                f"option string and `--captain-ruling {T1_R66_ROW}`. With --captain-ruling "
                                f"{T1_R73_ROW} only: row {T1_R73_ROW}'s four files (T1's second refinement round, each "
                                f"one further change, g_mom's damping or Blood's step, to {T1_R73_BASE_FILE}'s "
                                "command), exactly one of these sets: "
                                + "; ".join(f"{f} `{v[0]}`" for f, v in T1_R73_FILES.items())
                                + f", while row {T1_R73_ROW}'s first cell is the captain's order naming the set's "
                                f"option string paired with its file and `--captain-ruling {T1_R73_ROW}`. Each is "
                                "refused under any other ruling")
            p.add_argument("--held-momentum", default=None, metavar="{0.5,1.0,2.0}",
                           help="Blood's released momentum on the banked file (off by default; files "
                                f"{', '.join(T1_HM_FILES.values())}): on the command {T1_CANDIDATE_COMMAND} "
                                f"--feng-fraction-mult {T1_BANKED_MULT:g} (the banked file {T1_BANKED_FILE}) at the "
                                "default seed on the test board, the drawn cells labelled Blood also get Blood's own "
                                "released E9.5 -> E10.5 change (B0's damped mean momentum fitted on Blood-labelled "
                                "cells of E9.5 and E10.5's libraries 1 and 2) at this damping, after the atlas step "
                                "and before g_mom; every other cell is the banked file's. Takes the literal values "
                                f"{', '.join(T1_HM_CHOICES)} only, each only with --captain-ruling <ROW>, ROW a row "
                                "whose first cell is the captain's ruling naming the value's option string (e.g. "
                                f"{T1_HM_OPTIONS[0.5]}) and `--captain-ruling <ROW>`; refused with "
                                "--direction-momentum, --g-mom-damping, --released-bet, --bet-strength or "
                                "--no-hold-unbracketed and with another multiplier. With --feng-fraction-mult "
                                f"{T1_MSB_MULT:g} --direction-momentum all and --captain-ruling {T1_R51_ROW}, 0.5 and "
                                f"1.0 are row {T1_R51_ROW}'s files T1_s2rcMSBH05 and T1_s2rcMSBH10 (see "
                                "--captain-ruling); with --feng-fraction-mult 2.0 --direction-momentum all and "
                                f"--captain-ruling {T1_R60_ROW}, 0.5 and 1.0 are row {T1_R60_ROW}'s files "
                                "T1_s2rcMS20H05 and T1_s2rcMS20H10, and with --feng-fraction-mult 2.5 "
                                f"--direction-momentum all and --captain-ruling {T1_R46_ROW}, 0.5 is row "
                                f"{T1_R46_ROW}'s file {T1_R46_FILE}; with --feng-fraction-mult 2.0 "
                                f"--direction-momentum all and --captain-ruling {T1_R66_ROW}, 0.5 and 1.0 are in "
                                f"row {T1_R66_ROW}'s two- and three-change files and 2.0 is its file "
                                "T1_s2rcMS20H20 (see --direction-momentum); with --feng-fraction-mult 1.5 "
                                "--direction-momentum all --direction-momentum-damping 1.0 and --captain-ruling "
                                f"{T1_R73_ROW}, 0.5 and 1.0 are row {T1_R73_ROW}'s files T1_s2rcMSBm10H05 and "
                                "T1_s2rcMSBm10H10 (see --direction-momentum)")
            p.add_argument("--direction-momentum-damping", default=None, metavar="{1.0}",
                           help=f"row {T1_R51_ROW}'s file {T1_R51_M10_FILE} (off by default): MS's own damping (the "
                                "per-category momentum of --direction-momentum, harness.B0_DAMPING in every other MS "
                                "file) at 1.0, on the command "
                                + T1_R51_COMMAND.format(option=T1_R51_FILES[T1_R51_M10_FILE][0])
                                + f" at the default seed on the test board only, while row {T1_R51_ROW}'s first cell "
                                "is the captain's ruling naming that option string; takes the literal value 1.0 only, "
                                "and is refused with every other command. Row "
                                f"{T1_R58_ROW}'s file {T1_R58_FILE}: the literal value {T1_R58_DMD}, only as "
                                f"{T1_R58_COMMAND} at the default seed on the test board, while row {T1_R58_ROW}'s "
                                f"first cell is the captain's order naming that option string. Row {T1_R67_ROW}'s "
                                "three files: "
                                + "; ".join(f"{f} the literal value {v[1]}" for f, v in T1_R67_FILES.items())
                                + ", each only as "
                                + T1_R67_COMMAND.format(option=T1_R67_PREFIX + "<value>")
                                + f" at the default seed on the test board, while row {T1_R67_ROW}'s first cell is "
                                "the captain's order naming that option string paired with the file. Row "
                                f"{T1_R73_ROW}'s four files: the literal value 1.0 with one further step each, g_mom's "
                                "damping or Blood's (see --direction-momentum)")
            p.add_argument("--late-mix", choices=LATE_MIX_CHOICES, default=None,
                           help=f"row {T1_R54_ROW}'s file {T1_R54_FILE} (off by default): T1_s2rcMSB's command with "
                                "library 3's own chamber mix moved to the shares row "
                                f"{LATE_MIX_ROW}'s step commits ({LATE_MIX_SHARES_PATH}: "
                                "cardiomyocytes of heart cells, then blood of chamber cells, the other categories in "
                                f"proportion), only as {T1_R54_COMMAND} at the default seed on the test board, while "
                                f"row {T1_R54_ROW}'s first cell is the captain's ruling naming "
                                f"{T1_R54_OPTION} and the shares file committed at HEAD has the sha256 row "
                                f"{LATE_MIX_ROW} names (`shares.json sha256 <hex>`); refused with every other command")
            p.add_argument("--bet-strength", type=float, default=None,
                           help="with --released-bet: lambda (share_trend) or tau (advance), one of row T61's values")
            p.add_argument("--captain-ruling", default=None,
                           help="with --released-bet: the ID of the row of experiments/harness/TRIES.md that is the "
                                "captain's ruling for the arm; with --feng-fraction-mult "
                                f"{FENG_FRACTION_MULT_R20:g}: {T1_R20_ROW}, the captain's order for that file; with "
                                f"{FENG_FRACTION_MULTS_R22[0]:g} or {FENG_FRACTION_MULTS_R22[1]:g}: {T1_R22_ROW}, the "
                                "captain's order for the first search line; with --g-mom-damping "
                                f"{G_MOM_DAMPINGS_R22[0]:g} or {G_MOM_DAMPINGS_R22[1]:g}: {T1_R22_ROW}, the captain's "
                                "order for the second T1 search line; with --g-mom-damping 1.0 and "
                                f"--feng-fraction-mult {FENG_FRACTION_MULTS_LINE3[0]:g} or "
                                f"{FENG_FRACTION_MULTS_LINE3[1]:g}: {T1_R22_ROW}, the captain's order for the third "
                                f"T1 search line; with --g-mom-damping {G_MOM_DAMPINGS_LINE5[0]:g} or "
                                f"{G_MOM_DAMPINGS_LINE5[1]:g} and --feng-fraction-mult {FENG_FRACTION_MULT_LINE5:g}: "
                                f"{T1_R22_ROW}, the captain's order for the fifth T1 search line; with "
                                f"--direction-momentum all: {T1_R23_ROW}, the captain's order "
                                f"for candidate MS; with --feng-fraction-mult {T1_MSB_MULT:g} --direction-momentum "
                                f"all: {T1_R26_ROW}, the captain's order for the combined MS file; with "
                                f"--g-mom-damping {G_MOM_DAMPINGS_LINE6[0]:g} or {G_MOM_DAMPINGS_LINE6[1]:g} and "
                                f"--feng-fraction-mult {FENG_FRACTION_MULT_LINE6:g}: {T1_R22_ROW}, the captain's "
                                "order for the sixth T1 search line; with --g-mom-damping "
                                f"{G_MOM_DAMPING_LINE7:g} and --feng-fraction-mult {FENG_FRACTION_MULT_LINE7:g}: "
                                f"{T1_R22_ROW}, the captain's order for the seventh T1 search line; with --base "
                                f"pooled --feng-groups recomputed: {T1_R37_ROW}, the captain's named exception for "
                                "R13's command on the pooled base; with --held-momentum, or with --feng-fraction-mult "
                                f"{T1_BANKED_MULT:g} --direction-momentum all: the ID of the row whose first cell is "
                                "the captain's ruling naming that file's option string and this flag with that ID; "
                                f"{T1_R51_ROW} with exactly one of row {T1_R51_ROW}'s seven option sets on "
                                f"{T1_CANDIDATE_COMMAND} (T1's line 10, each one change to {T1_R51_BASE_FILE}: "
                                + "; ".join(f"{f} `{v[0]}`" for f, v in T1_R51_FILES.items())
                                + f"), while that row's first cell is the captain's ruling naming the set's option "
                                f"string and `--captain-ruling {T1_R51_ROW}` as whole backticked strings. With "
                                f"{T1_R58_OPTION_STRING}: {T1_R58_ROW}, the captain's order for {T1_R58_FILE}. With "
                                f"--feng-fraction-mult {T1_MSB_MULT:g} --direction-momentum all --late-mix "
                                f"{LATE_MIX_CHOICES[0]}: {T1_R54_ROW}, the captain's named exception for "
                                f"{T1_R54_FILE}. With exactly one of row {T1_R66_ROW}'s thirty option sets (see "
                                f"--direction-momentum): {T1_R66_ROW}, the captain's order for T1's third-round "
                                f"step-on slots. With exactly one of row {T1_R67_ROW}'s three option sets on "
                                f"{T1_CANDIDATE_COMMAND} ("
                                + "; ".join(f"{f} `{v[0]}`" for f, v in T1_R67_FILES.items())
                                + f"): {T1_R67_ROW}, the captain's order for those three files, while that row's "
                                "first cell names the set's option string paired with its file and "
                                f"`--captain-ruling {T1_R67_ROW}`. With exactly one of row {T1_R73_ROW}'s four option "
                                f"sets on {T1_CANDIDATE_COMMAND} (see --direction-momentum): {T1_R73_ROW}, the "
                                "captain's order for those four files, while that row's first cell names the set's "
                                f"option string paired with its file and `--captain-ruling {T1_R73_ROW}`")
        if name == "heart":
            p.add_argument("--entry", choices=("s1", "s2", "s2_trim", "core_trim", "organ_scope"), required=True,
                           help="s1 = H1; s2 = H2 as first built (size only); s2_trim = H2 as a trim to half the "
                                "scope, the s2 of docs/S2_RULES.md; core_trim = the full trim (no ruling; only as "
                                f"{HEART_CORE_TRIM}, sized to the raw mode's own core forecast); organ_scope = the "
                                f"organ cut, file D (row R21; only as {HEART_ORGAN_SCOPE}, sized to the raw-mode core "
                                "forecast times rms(organ pool) / rms(core) on E10.5)")
            p.add_argument("--herd", action="store_true",
                           help="herd B0's draw towards the pool's mean expression (try H4; off: the entry so far)")
            p.add_argument("--depth", action="store_true",
                           help="thin B0's cells to the extrapolated E12.5 depth (heart try H3 of "
                                "experiments/harness/TRIES.md; off: the entry as adopted so far)")
            p.add_argument("--cfm", action="store_true",
                           help="filter each lineage's momentum by its own within-stage covariance (heart try H6 of "
                                "experiments/harness/TRIES.md; off: the entry as adopted so far)")
            p.add_argument("--trend", action="store_true",
                           help="b0_trend: each lineage's momentum the least-squares trend over E8.25_late, E8.5, "
                                "E8.75, E9.5 and E10.5, over B0's time budget (the captain's ruling of 2026-11-13, "
                                "docs/PLAN.md decision log (w); ledger row H11 of experiments/harness/TRIES.md); "
                                "with --entry s1 or s2_trim, refused with --entry s2, --depth, --herd and --cfm")
            p.add_argument("--expression", choices=("b0", "copy"), default="b0",
                           help="b0 = B0's momentum step (the default, the entry as built); copy = the copy baseline's "
                                "cells and expression, B0's own draw with no step (the captain's direction for the "
                                "heart's second upload, docs/PLAN.md decision log 2026-11-18 (k)); --entry s1, or "
                                "--entry s2_trim with --z-mode raw (the scope bet, the same draw from the half trim "
                                "sized for H2; no ruling, scope is a closed lever, only the captain's ruling puts it "
                                "in an upload), or --entry core_trim with --z-mode raw (the full trim, likewise no "
                                "ruling); refused with --entry s2, --trend, --depth, --herd and --cfm")
            p.add_argument("--z-mode", choices=("recorded", "raw"), default="recorded",
                           help="the size step's z mode: recorded = the mode the size JSON pre-registered (the "
                                "default); raw = the isotropic mode the JSON also logs, one factor on all three axes "
                                "to that raw mode's own target RMS for the hypothesis (z_modes.raw, not the recorded "
                                "mode's pre-registered target); refused when the size JSON is not valid (its "
                                "fallback is already isotropic: build with recorded then)")
            p.add_argument("--run-profile", default=None, metavar="SHA256",
                           help="file OXP125 (row R56 of experiments/harness/TRIES.md, the captain's order of "
                                "2026-11-24; off by default): OTS125's command with E10.5's run profile removed, "
                                f"{HEART_XPC_COMMAND.format(sha256='<sha256>')}, the sha256 of "
                                "experiments/heart_profile/results/xpc_profile.json (row H36's test profile). Refused "
                                "with any other command, seed, board or option, and unless row R56's first cell is "
                                "the captain's order naming --entry organ_scope, --trend-strength 1.25, --run-profile "
                                "and --captain-ruling R56")
            p.add_argument("--trend-strength", default=None, metavar="{0.5,1.0,0.75,1.25,1.5}",
                           help="the heart search line (row R28 of experiments/harness/TRIES.md, the captain's order "
                                "of 2026-11-23, \"start the heart search\"; off by default): the full trim C's command "
                                "with row H28's trend step at strength 0.5 (TS05) or 1.0 (TS10), nothing else changed: "
                                f"{HEART_TS_COMMANDS[0.5]}. Refused with any other value, entry or z mode, with "
                                "--expression copy, --trend, --depth, --herd and --cfm, on any board but the test "
                                "board, and unless row R28's first cell is the captain's order naming both values and "
                                "--captain-ruling R28. File D's line (row R29, off by default): "
                                f"{HEART_OTS_COMMANDS[0.5]} (OTS05) or {HEART_OTS_COMMANDS[1.0]} (OTS10), D's command "
                                "with the same step, refused alike and unless row R29's first cell is the captain's "
                                "order naming --entry organ_scope, both values and --captain-ruling R29. File D's "
                                "two further strengths (off by default): "
                                f"{HEART_FURTHER_COMMANDS[0.75].format(row=HEART_FURTHER_ROW)} (OTS075) or 1.25 "
                                "(OTS125), D's command with the same step, refused alike and unless TRIES.md holds "
                                f"exactly one row {HEART_FURTHER_ROW} whose first cell is the captain's order naming "
                                f"--entry organ_scope, the two values and --captain-ruling {HEART_FURTHER_ROW} (any "
                                "other row refused). File D's strength 1.5 (row R52, off by default): "
                                f"{HEART_OTS15_COMMAND} (OTS15), D's command with the same step, refused alike and "
                                "unless row R52's first cell is the captain's order naming --entry organ_scope, 1.5 "
                                "and --captain-ruling R52")
            p.add_argument("--captain-ruling", default=None,
                           help=f"with --trend-strength only: {HEART_TS_ROW} with --entry core_trim, {HEART_OTS_ROW} "
                                "with --entry organ_scope, the row that orders it; with --entry organ_scope "
                                f"--trend-strength 0.75 or 1.25: {HEART_FURTHER_ROW}, the ruling row that orders those "
                                f"two; with --entry organ_scope --trend-strength 1.5: {HEART_OTS15_ROW}")
            p.add_argument("--insitu-chords", default=None, metavar="SHA256",
                           help="file OIC125 (row R70 of experiments/harness/TRIES.md, the captain's order of "
                                "2026-11-25; off by default): OTK125's command with the in-situ chords added, "
                                f"{HEART_OIC_COMMAND.format(sha256='<sha256>')}, the sha256 of "
                                "experiments/heart_ext/results/chords_E10.5.json (row H38's G0). Refused with any "
                                "other command, seed, board or option, unless row R70's first cell is the captain's "
                                "order naming OIC125 and --captain-ruling R70, and unless row H38's outcome cell "
                                "holds G0's PASS statement with this sha256 and then the read's PASS")
            p.add_argument("--trend-form", default=None, metavar="{" + ",".join(HEART_TF_FORMS) + "}",
                           help=f"row {HEART_TF_ROW}'s trend-form line (off by default): OTS125's command with one "
                                "change, the direction of the trend step (vec.trend_forms): category (OTK125, the "
                                "Endothelial lineage split by chamber category), contrast (OTX125, the shared change "
                                "out, norms held), whole (OTW125, E10.5's point on the whole stage) or portion "
                                "(OTP125, the slope before E10.5, norms held); only as "
                                + HEART_TF_COMMANDS["category"].replace("category", "<F>")
                                + f" at the default seed on the test board, while row {HEART_TF_ROW}'s first cell is "
                                "the captain's order naming --entry organ_scope, --trend-strength 1.25, the form and "
                                f"--captain-ruling {HEART_TF_ROW} as whole backticked strings. Row {HEART_R74_ROW}'s "
                                "four files (off by default): "
                                + HEART_R74_COMMAND.format(option="<OPTION>")
                                + " with <OPTION> one of "
                                + ", ".join(f"{v[0]} ({f})" for f, v in HEART_R74_FILES.items())
                                + f" (its table), at the default seed on the test board, while row {HEART_R74_ROW}'s "
                                "first cell is the captain's order naming --entry organ_scope, --z-mode raw, the "
                                f"option string followed by its file and --captain-ruling {HEART_R74_ROW} as whole "
                                "backticked strings")
        if name == "embryo":
            p.add_argument("--method", choices=("pseudobulk_shift", *EMBRYO_T_POS), required=True,
                           help="fgw = s1; fgw_shape0 = the s2 of docs/S2_RULES.md (E7.5's geometry)")
            p.add_argument("--herd", action="store_true",
                           help="herd each bracket's draw towards its mean expression (try E9; off: the entry so far)")
            p.add_argument("--expression",
                           choices=("mixture", "depth", "blend", "latent", "series", "render_hi", "series_back"),
                           default="mixture",
                           help="fgw methods: the cells as released (the entry), or tries E7, E8, E14, E24 and E19c "
                                "of experiments/harness/TRIES.md: thinned to the predicted depth; blends; a count "
                                "model's latent interpolation (needs the gpu dependency group); the shallower "
                                "bracket's cells in the deeper bracket's run profile; every cell topped up to the "
                                "deeper bracket's depth. Neither of the last two is an entry until the captain rules. "
                                "series_back (row R48, the captain's order of 2026-11-24; off by default): the mirror "
                                "of series, the E8.0 rows in E7.5's run profile and depth and the E7.5 rows as "
                                f"released, only as {EMBRYO_OVSB_COMMAND} (O's command with it). With --shape-interp "
                                f"overlay, depth is row {EMBRYO_R75_ROW}'s file T2embryo_s2x_ovd.h5ad alone and series "
                                f"with --mix-lean later its file T2embryo_s2xs_ovt23.h5ad, each only with "
                                f"--captain-ruling {EMBRYO_R75_ROW} (see --mix-lean)")
            p.add_argument("--move-coupling", choices=t2_build.MOVE_COUPLINGS, default=t2_build.MOVE_COUPLINGS[0],
                           help="fgw methods: which coupling sets the moves: fgw (the entry) or aligned (try E10 of "
                                "experiments/harness/TRIES.md: an exact transport in the rotated common frame)")
            p.add_argument("--draw", choices=t2_build.DRAWS, default=t2_build.DRAWS[0],
                           help="fgw methods: which real cells of each bracket are drawn: uniform (the entry) or "
                                "midpoint (try E13 of experiments/harness/TRIES.md: tilted towards the local "
                                "displacement midpoint, every cell unedited at its moved position)")
            p.add_argument("--positions", choices=t2_build.POSITIONS, default=t2_build.POSITIONS[0],
                           help="--method fgw_shape0 only: where the cells sit: moved (the default, the entry) or "
                                "floor_exact (tries E47a and E47b of experiments/harness/TRIES.md: every cell on a "
                                "real E7.5 position, the E8.0 cells on unused ones by an exact expression assignment; "
                                "with --expression series the harness's fgw_shape0x_series). No entry until the "
                                "captain rules; refused with --herd and --draw midpoint")
            p.add_argument("--row-order", choices=t2_build.ROW_ORDERS, default=t2_build.ROW_ORDERS[0],
                           help="as_built (the default, every file as built) or copy_parity: the rows re-ordered so "
                                "the file's pre-fix SVD parity on the scorer's matched rows is the floor copy's "
                                "(vec.t2_build.copy_parity_order; branch campaign/embryo-occupancy). Same cells, "
                                "expression and positions. No try and no ruling: only the captain's ruling puts it "
                                "in an upload. It assumes the server scores the prediction's rows as the file holds "
                                "them (file order, full count, no subsample or sort); score_h5ad's _prediction does "
                                "this in the pinned copy, but the live server's handling is not shown")
            p.add_argument("--shape-interp",
                           choices=(SHAPE_INTERP_OFF, *BUILDABLE_SHAPES),
                           default=SHAPE_INTERP_OFF,
                           help="--method fgw_shape0 --positions floor_exact only: none (the default, every file as "
                                "built) or one shape form of vec.shape_interp between E7.5 and E8.0 at t = 1/2 on the "
                                "exact-position cloud, before the size rule (row E49 as drafted), or overlay_stretch "
                                "(row E55 as specified: the overlay's union stretched to the brackets' interpolated "
                                "profile). No try read and no ruling: only the captain's ruling puts it in an upload")
            p.add_argument("--shape-interp-test-ok", action="store_true",
                           help="--shape-interp on the test board: the captain's yes (off by default; without it the "
                                "option refuses on a test board before anything loads)")
            p.add_argument("--copy-rows", action="store_true",
                           help="--method fgw_shape0 --positions floor_exact only, with --expression mixture or "
                                "series (off by default; the captain's order of 2026-11-21, docs/PLAN.md decision log "
                                "2026-11-21 (h); row E50 as specified): the file takes the floor copy's own rows, "
                                "order, count and float32 coordinates (vec.evaluate.subsample of E7.5 as released, "
                                "0.1, seed 0), half its slots keep their own cell's expression and the rest take E8.0 "
                                "cells by an exact assignment (vec.t2_build.copy_rows_cloud). It departs from the "
                                "count rule and the size rule. Refused with --row-order copy_parity, --shape-interp, "
                                "--herd, --draw midpoint and --move-coupling aligned")
            p.add_argument("--copy-rows-test-ok", action="store_true",
                           help="--copy-rows on the test board: the captain's yes (off by default; without it the "
                                "option refuses on a test board before anything loads, because a refusal or a pass "
                                "there shows one bit about the board's limits)")
            p.add_argument("--copy-rows-size", action="store_true",
                           help="with --copy-rows and --expression mixture only (off by default; the captain's order "
                                "of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (v); row E52 as specified): the "
                                "copy-rows file with its coordinates scaled about the copy's own centroid to the "
                                "embryo size JSON's pre-registered full-cloud RMS (raw mode), the last step before "
                                "writing (vec.t2_build.copy_rows_sized)")
            p.add_argument("--copy-rows-size-test-ok", action="store_true",
                           help="--copy-rows-size on the test board: the captain's yes (off by default; a test build "
                                "needs it and --copy-rows-test-ok, before anything loads)")
            p.add_argument("--row-order-test-ok", action="store_true",
                           help="--row-order copy_parity on the test board: the captain's yes (off by default; "
                                "without it the option refuses on a test board before anything loads, because a "
                                "refusal or a success there shows one bit of the test count)")
            p.add_argument("--mix-lean", choices=(*MIX_LEANS, *EMBRYO_R75_NEW_LEANS, *EMBRYO_R78_NEW_LEANS),
                           default=None,
                           help=f"row {EMBRYO_OT23_ROW}'s file {EMBRYO_OT23_FILE} (the captain's order of 2026-11-24; "
                                "off by default): O's command with the later bracket's share of the drawn cells at "
                                "t_mix = 2/3 exactly (t = 1/2 still sets the moves), only as "
                                f"{EMBRYO_OT23_COMMAND}; it moves the overlay's form (the union leans to E8.0's "
                                f"specimen). Row {EMBRYO_R78_ROW}'s files (the captain's order of 2026-11-26): "
                                + "; ".join(f"`{o}` ({f})" for f, (o, *_) in EMBRYO_R78_FILES.items())
                                + f", each only as {EMBRYO_R78_COMMAND.format(option='<OPTION>')}, t_mix "
                                + ", ".join(f"{t.numerator}/{t.denominator}" for t in EMBRYO_R78_LEANS.values())
                                + f". Row {EMBRYO_R75_ROW}'s files (the captain's order of 2026-11-26): "
                                + "; ".join(f"`{o}` ({f})" for f, (o, *_) in EMBRYO_R75_FILES.items())
                                + f", each only as {EMBRYO_R75_COMMAND.format(option='<OPTION>')}, t_mix "
                                + ", ".join(f"{t.numerator}/{t.denominator}" for t in EMBRYO_R75_LEANS.values()))
            p.add_argument("--neighbourhood", choices=NEIGHBOURHOODS, default=NEIGHBOURHOODS[0],
                           help="as_built (the default, every file as built) or run_sorted (row E57's file "
                                f"{EMBRYO_RS_FILE}, off by default): O's rows and positions with each bracket's "
                                "expression rows moved to its own side of the earlier bracket's midline within fixed "
                                "lineage pools (vec.run_sort), the last step, on the final coordinates. Only as "
                                f"{EMBRYO_RS_COMMAND.format(row='<ROW>')}, with a ruling row of the captain naming "
                                f"{EMBRYO_RS_OPTION} and `--captain-ruling <ROW>`")
            p.add_argument("--captain-ruling", default=None, action=_OnceAction,
                           help=f"with --expression series_back: {EMBRYO_OVSB_ROW}, the row that orders it; with "
                                "--neighbourhood run_sorted: the captain's ruling row that orders it; with --mix-lean "
                                f"later: {EMBRYO_OT23_ROW}, the row that orders it; with --shape-interp overlay_bayes: "
                                f"{EMBRYO_BF_ROW}, the row that orders it; with one of row {EMBRYO_R78_ROW}'s two "
                                f"option strings on O's command: {EMBRYO_R78_ROW}; with one of row "
                                f"{EMBRYO_R75_ROW}'s five "
                                f"option strings on O's command: {EMBRYO_R75_ROW}; given once")
    return ap


def main(argv=None) -> None:
    a = build_parser().parse_args(argv)
    _check_seed(a.seed)
    if a.board not in (None, TEST_BOARDS[a.which]):  # a dry run
        if a.board not in DRY_BOARDS[a.which]:
            sys.exit(f"--board {a.board}: a {a.which} dry run takes only {sorted(DRY_BOARDS[a.which])}")
        if "dry" not in Path(a.out).name.lower():
            sys.exit(f"{a.out}: a dry run's file name must contain 'dry' (the file must never be uploaded)")
    elif "dry" in Path(a.out).name.lower():
        sys.exit(f"{a.out}: a test build's file name must not contain 'dry'")
    if getattr(a, "strong_null", False) and Path(a.out).exists():
        sys.exit(f"{a.out} exists: a build with --strong-null never overwrites a file (T12's row)")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    rec = {"t1": build_t1, "heart": build_heart, "embryo": build_embryo}[a.which](a)
    rec = _hide_counts(rec, load_boards()[rec["board"]])
    Path(a.out).with_suffix(".build.json").write_text(json.dumps(rec, indent=1, default=float))
    print(json.dumps({k: v for k, v in rec.items() if k not in ("shares", "external", "realised", "strong_null")},
                     default=float))


if __name__ == "__main__":
    main()
