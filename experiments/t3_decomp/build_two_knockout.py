#!/usr/bin/env python
"""The barred branch's T3 file with its Mab21l2 factor multiplied by the released Gata4 knockout's factor: the forms of
the two-knockout list (ledger row R17; the decision-log entry that fixes the list before any of its scores).

A thin builder, in the manner of build_barred_exponent.py. It imports the frozen builder (build_submission.py) and
script (t3_decompose.py) unchanged and builds a barred-branch file exactly as `build_submission.py build DECISIONS
--branch barred --submission s1|s2` does, with one change: the factor the frozen loader returns for the base's table is
replaced, for this process only, by a fixed combination of two committed factor tables, both loaded and sha256-checked
by the frozen loader (`build_submission.load_factor`) and checked again against the sha256 pinned below:

  phi_M  <run>/factors/transfer_mab21l2_all.csv   the released Mab21l2 knockout against WT E8.75 (the barred factor)
  phi_G  <run>/factors/gata4_pooled.csv            the released Gata4 knockout (both validation replicates) against
                                                   WT E8.75, strict non-Mesp1 labels (the allowed branch's factor)

Each table is used as loaded: t3_decompose.gene_factor already rescaled it (den-weighted total proportion unchanged)
and clipped it to [0.25, 4] when the run wrote it. After the combination nothing is rescaled or clipped again, as the
exponent builder (row R14) and the frozen build treat every exponent. So the combined factor is not rescaled: its
common shift is recorded in the sidecar (gene mean of ln phi, cos^2 of ln phi with the all-ones vector), not treated.

The forms (--form; nothing else is accepted, and no number is taken from the command line):

  identity-s1        phi_M                           = the frozen barred s1 byte for byte (3712c3f9988d...)
  identity-s1x1.5    phi_M ** 1.5                    = build_barred_exponent.py --exponent 1.5 (c060fcaa6c7f...)
  identity-s1x2.0    phi_M ** 2.0                    = build_barred_exponent.py --exponent 2.0 (7393810fbe3b...)
  identity-s2        phi_M                           = the frozen barred s2 (1b69521e3c28...)
  c1-s1, c1-s2       phi_M * phi_G                   (C1: base x = 1 only)
  c2-s1, c2-s1x1.5,  phi_M ** (x c) * phi_G ** c     (C2: held at the base's size)
  c2-s1x2.0, c2-s2   c = SD(x ln phi_M) / SD(x ln phi_M + ln phi_G), population SD over the panel's genes,
                     computed here from the two tables as loaded, never from a score

C1 on the exponent-1.5 and exponent-2.0 bases (c1-s1x1.5, c1-s1x2.0) is refused: the list fixes C1 on x = 1 only.

Where each quantity enters. The base fixes the slot (s1: the central composition at m = D7m, cells drawn on stream
"submission:barred:s1"; s2: the mild composition, stream "submission:barred:s2"), the "...:assemble" stream, the
seed and the board; all are the frozen builder's own, so the same wild-type cells are drawn in every form of a base.
The spec keeps the frozen factor {mab21l2, exponent 1.0}, so t3_decompose.build receives w = 1 and applies the
combined phi as given (`apply_gene_factor`, after the cells are drawn, before lognorm). The base's exponent x enters
as `phi_M ** x` here, the same numpy operation the frozen build applies for w = x, so each identity form gives its
base's bytes. No depth, stage or coverage step: the barred spec has none.

Refusals (exit status 2), before anything loads: another form (taken as typed, no whitespace stripped); on the command
line a seed other than 1000001 or an --n-cells other than 'max' (the list fixes both); t3_decompose.py, rules.yaml
or build_submission.py whose sha256 differs from the bytes pinned below (the same bytes build_barred_exponent.py
pins); a tree that is not clean by vec.provenance's rule; the frozen builder's decision-file check, with no allow
flag; then, when the factor loads, a table whose sha256 differs from the decision file's record (the frozen loader) or
from the value pinned here, a table whose genes differ from the reference's (order included; the frozen loader), or a
combination that is not finite and positive.

The sidecar is the frozen builder's plus one block `combined_build` (this script and its sha256, the form, the base,
x, the Gata4 weight, c, both tables with their paths and hashes, the SDs and correlation, the common-shift figures, the
pinned frozen hashes, the list's entry). It adds no cell number. `inputs.factor` stays the frozen loader's record of
phi_M; the factor applied is the one `combined_build` names, as `inputs.factor_note` says in every form.

  uv run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build DECISIONS.json \\
      --form c1-s1 --wt data/raw/T2/E8.75.h5ad --run-dir runs/t3_decomp/gata4_<STAMP> --board T3:bcatenin \\
      --seed 1000001 --out runs/<folder>/T3_c1_s1.h5ad

--gata4-exponent 2.0 --captain-ruling R22 (row R22's second T3 search line, the Gata4-table weight w on C1, as
amended by row R30; off by default: without both, every form, byte and behaviour above is unchanged and nothing below
is read or computed): C1 on the barred s2 with its Gata4 weight raised from 1 to 2, the one file T3_c1_s2_g2.h5ad,

  c1-s2 --gata4-exponent 2.0    phi_M ** 1 * phi_G ** 2.0     nothing rescaled or clipped after the product

applied where C1's factor is applied (bs.load_factor's replacement, at w = 1, once, after the cells are drawn).
Refused (exit 2, nothing built): any other value; any other form; --gata4-exponent without --captain-ruling R22;
--captain-ruling R22 without --gata4-exponent 2.0; any other ruling; and the pair while experiments/harness/TRIES.md
does not hold exactly one row R22 whose first cell passes build_outside_form.captain_order_names for the option
"`--gata4-exponent 2.0`". Before the build, from released WT E8.75 (den as build_outside_form.wt_den makes it) and
the two pinned tables only, the line's figures are computed and recorded: the den-weighted totals of phi_M, of
phi_M * phi_G (C1's) and of phi_M * phi_G ** 2, and the two relative shifts against phi_M's. By row R30 the line's old
5% limit gates nothing: no value of a shift stops the build (the old value is recorded as lifted). The factor the build
applies must be the very array those figures were computed from (its sha256, over WT E8.75's genes); if it is not,
exit 5 with nothing written. The sidecar's combined_build then also holds `gata4_line` (the ruling row, the line, the
weight, row R30's lifting, the figures and the sha256 of the factor applied).

--composition-from allowed-s1 --captain-ruling R36 (row R36, the captain's named exception of 2026-11-24; off by
default: without both, every form, byte and behaviour above is unchanged and nothing below is read or computed): the
file cD7, T3_c2_s2_cd7.h5ad, which is C2 on the barred s2 (`--form c2-s2`, recorded sha256 47111ac1d07e...) rebuilt by
its own command with one change: its composition block is replaced by submissions.allowed.s1.composition of the
decision file, read by code, so that the five scenario states go to 2^-1.3826521560990648 = 0.3835131214774091 of their
wild-type share (before the composition routine renormalises). Only that block is copied: the same entry's factor
(gata4) and depth step (depth: true, a closed lever) are never taken. C2's factor, row number (max_cells, read by
code), seed, streams and coordinates rule are unchanged, and C2 has no depth step. The replacement is made where the
frozen builder resolves the barred s2's spec (bs.resolve_spec, for this process only), so the frozen validation,
composition routine, draw and assembly run on it as on any spec; at the build the composition routine's result is
checked as G1 checks it, before any cell is drawn.
Refused (exit 2, nothing built): --captain-ruling R36 without the option; any value but allowed-s1; any form but
c2-s2; the option without --captain-ruling R36 or with another ruling; --gata4-exponent beside it; a ledger
(experiments/harness/TRIES.md) without exactly one row R36 whose first cell opens "R36 (the captain's named
exception", says "ruling, no run" and names `--composition-from allowed-s1`, `--form c2-s2` and `--captain-ruling
R36` (the check build_outside_form.captain_order_names makes of an order); a decision file whose sha256 is not the
pinned one (9440ae67...); a copied block whose states differ from rules.yaml's scenario states, or whose log2
magnitude differs from rule D7's m (the decision file's rules.D7.m, 1.3826521560990648), or that is not a central
composition at depletion factor 2 ** -m; a composition routine that does not give the five states at that factor of
their wild-type share. The sidecar's combined_build then also holds `composition_line` (the ruling row, the block and
the one it replaced, the five ratios).

--form c2-s2 --gata4-exponent 2.0 --captain-ruling R22 (row R22's third T3 search line, decision log 2026-11-24 (kk),
on the captain's own yes naming T3_c2_s2_g2.h5ad, ruling row R41; off by default: without the option on c2-s2, every
form, byte and behaviour above is unchanged, combine() included, and nothing below is read or computed): the file
T3_c2_s2_g2.h5ad, C2's command with the gata4_pooled weight g inside C2's size-held product raised from 1 to 2,

  c2-s2 --gata4-exponent 2.0   phi_M ** c * phi_G ** (g c), c = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), g = 2.0

population SD over the genes of the two pinned tables as loaded; nothing rescaled or clipped after the product; at
every g the SD of ln factor is phi_M's; g = 0 gives phi_M and g = 1 C2's factor (combine('c2', 1.0)) value for value.
On c2-s2 the option's value is the weight g, not a power: the powers applied are c on transfer_mab21l2_all and g c on
gata4_pooled. Cells, streams, composition, seed, count and board stay C2's (the barred s2's). Refused (exit 2, nothing
built): any other value on c2-s2 (3.0 included); the option without --captain-ruling R22 or with another ruling; a
ledger without exactly one row R22 whose first cell passes build_outside_form.captain_order_names for
"`--form c2-s2 --gata4-exponent 2.0`" (line 2's entry there does not open this line); and, as before, the option on
any form but c1-s2 (line 2) and c2-s2, or either flag alone. The guard, before the build, from released WT E8.75 (den
as build_outside_form.wt_den makes it) and the two pinned tables only: the den-weighted total of the file's factor
lies within 5% of C2's factor's; else exit 6, nothing written. C2's shift from phi_M's total and the file's shift from
phi_M's total are recorded and gate nothing. The factor the build applies must be the very array the guard describes
(its sha256, over WT E8.75's genes); else exit 5, nothing written. The sidecar's combined_build then also holds
`c2_gata4_line` (the ruling row, the line, g, c and the powers applied, the guard's figures and the factor's sha256).

--gata4-exponent inf and --size-scale 0.67 on --form c2-s2, each only with --captain-ruling R22 (the T3 design of the
design round of 2026-11-24, section 2, as its cold challenge fixes it; off by default: a command with none of
--size-scale, --den-ruling and --gata4-exponent inf on c2-s2 goes on exactly as above, combine() and c2_gata4_weighted()
included, and nothing below is read or computed): three files on C2's command, points of the family ln phi(g, s) =
s c_g (ln phi_M + g ln phi_G), c_g = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), population SD over the genes of the two
pinned tables as loaded, each with its option string,

  F1  T3_c2_s2_gonly.h5ad        `--form c2-s2 --gata4-exponent inf`                     g = inf, s = 1
  F2  T3_c2_s2_s067.h5ad         `--form c2-s2 --size-scale 0.67`                        g = 1 (C2's factor), s = 0.67
  F3  T3_c2_s2_gonly_s067.h5ad   `--form c2-s2 --gata4-exponent inf --size-scale 0.67`   g = inf, s = 0.67

g = inf is phi_G alone at C2's held size, ln phi = (SD(ln phi_M) / SD(ln phi_G)) ln phi_G (phi_M dropped: its power is
0); s multiplies ln phi after the size-held combination (phi ** s; s = 1 when the option is absent, and no power is then
taken); nothing rescaled or clipped after. Cells, streams, composition, seed, count and board stay C2's (the barred
s2's). Refused (exit 2, nothing built): --size-scale on another form or with any value but 0.67; --size-scale beside a
--gata4-exponent other than inf; either option without --captain-ruling R22, with another ruling or beside
--composition-from; --den-ruling on any other command; a ledger without exactly one row R22 whose first cell passes
build_outside_form.captain_order_names for the file's option string and holds it as one whole backtick-delimited
string. The whole-string match is needed because that check matches by containment and F1's string sits inside F3's;
no new string is or contains line 2's "`--gata4-exponent 2.0`" or line 3's "`--form c2-s2 --gata4-exponent 2.0`", and
those lines' own checks are unchanged. --gata4-exponent inf on any other form is refused as before (only 2.0). Two
checks, each a property of the file's own factor against C2's (the two pinned tables and released WT E8.75's den),
before any cell is drawn and before anything is written (exit 6, nothing written):
  range   the applied factor's max at most twice C2's factor's max and its min at least half C2's min
  5% den  line 3's check: the file's den-weighted total (den as build_outside_form.wt_den makes it) within 5% of C2's
          factor's. It gates unless --den-ruling ROW names a ruling row ("ruling, no run") whose first cell holds the
          file's option string as one whole backtick-delimited string and says "the 5% check gates nothing for it";
          then its figures are recorded and gate nothing (the lifting row R30 gave line 2's one file, keyed here to the
          row named on the command line and to the file by name). A --den-ruling row that does not is refused (exit 2).
The factor the build applies must be the very array the checks describe (its sha256, over WT E8.75's genes); else exit
5, nothing written. The sidecar's combined_build then also holds `c2_gonly_line` (the ruling row, the file, g, s, c and
the powers applied, both checks' figures, the lifting row if any and the factor's sha256).

--depth-from allowed-s1 on --form c2-s2, only with --captain-ruling R53 (the T3 ERROR design of the design round of
2026-11-25, section 3, proposal P1, dD, as its cold challenge fixes it; row R53, the captain's order of 2026-11-24: his
named exception on T3's depth for the three files, bound to row K30's read; off by default: a command with none of
--depth-from, --depth-exponent and --depth-form goes on exactly as above, every byte of C2 included, and nothing below
is read or computed): three files on C2's command (`--form c2-s2`), each C2 with one change, the frozen allowed
branch's per-label depth step (of the decision file's tables only submissions.allowed.s1 carries `depth: true`,
switched on by rule D4), each with its option string,

  dd   T3_c2_s2_dd.h5ad    `--depth-from allowed-s1`                        r_t, the frozen ratios
  ddh  T3_c2_s2_ddh.h5ad   `--depth-from allowed-s1 --depth-exponent 0.5`   r_t ** 0.5 (a new form)
  dg   T3_c2_s2_dg.h5ad    `--depth-from allowed-s1 --depth-form global`    one ratio for every label (a new form)

The step is added where the frozen builder resolves the barred s2's spec (bs.resolve_spec, for this process only, as
row R36's file replaces a composition block): C2's spec, its factor alone and a composition, gains `depth: true`, so
the frozen validation, bs.inputs_for (bs.run_parameters: the decision file's submissions.parameters) and
bs.prediction_inputs take it as the allowed s1 takes it. r_t is label t's ratio as the frozen bs.depth_ratios(ref.V,
params) reads it (depth_ratio_pooled; 1.0 for a label the frozen rule left unchanged); bs.depth_ratios, for this
process only, refuses ratios other than the pinned decision file's and returns the file's form; dg's one ratio is
exp(sum_t n_t ln r_t / sum_t n_t) over all of WT E8.75's labels (1.0 where none), n_t WT E8.75's cells in label t as
the build's reference holds them (duplicates dropped), computed here. The ratios enter as construct's `ratio`:
t3_decompose.build applies them after the factor (rescale_depth: binomial thinning below 1, Poisson top-up from the
kNN-smoothed profile above 1) on the same generator, so C2's draws (rows, factor) are untouched and only the last step
differs. Only the ratios are read: never the allowed entry's factor (gata4), composition, stage or coverage. Factor,
rows, composition, seed, count, streams, board and coordinates rule stay C2's. Refused (exit 2, nothing built):
--depth-exponent or --depth-form without --depth-from; any value of --depth-from but allowed-s1, of --depth-exponent
but 0.5 (1.0 included: dd takes no exponent), of --depth-form but global; --depth-exponent and --depth-form together;
any form but c2-s2; --composition-from, --gata4-exponent, --size-scale or --den-ruling beside it; no --captain-ruling,
or any but R53; a ledger without exactly one row R53 whose first cell passes build_outside_form.captain_order_names for
the file's option string (the captain's order, "ruling, no run", `--captain-ruling R53`) and holds that string and
`--form c2-s2` each as one whole backtick-delimited string (dd's sits inside ddh's and dg's); a decision file whose
sha256 is not the pinned one (9440ae67...), whose submissions.allowed.s1 lacks depth: true or whose ratios are not
finite and positive; a spec that is not C2's. Inside the build, before anything is written (exit 7, nothing written):
C2's cells are drawn again on a copy of the generator in the same state with no depth step, and the counts the depth
step receives must be C2's final counts row by row (the same rows, labels and sources: the same factor draws), and each
cell's ratio its label's. A cell whose thinning draw would be empty keeps its counts (t3_decompose.thin, the frozen
step's own rule): it stays unthinned, so no cell reaches lognorm empty through the step (lognorm's where=n > 0 guard
would write an empty row as zeros); those cells are counted on the very draw (a recording stand-in passes every call to
the real generator unchanged) and the number is recorded. The sidecar's combined_build then also holds `depth_line`
(the ruling row, the file, the exponent or form, the spec replaced, the ratios read and applied with their sha256, the
checks, the number of cells that stay unthinned) and its depth_step names it.

--captain-ruling R76 on --form c2-s2 (row R76, the captain's order of 2026-11-26: T3's refinement round on C2's
command, the T3 REF3 design as its cold challenge fixes it, the design's s11 struck; off by default: a command without
--captain-ruling R76 goes on exactly as above, every byte of C2 and of each earlier row's file included, and nothing
below is read or computed): five files on C2's command, each built from one of row R76's option strings, which one
literal table holds (R76_TABLE: the file, its option string as the row pairs them, and whether the row says "the 5%
check gates nothing for it" of it), each a point of the family ln phi(g, s) = s c_g (ln phi_M + g ln phi_G),
c_g = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), population SD over the genes of the two pinned tables as loaded,

  T3_c2_s2_s085.h5ad      `--form c2-s2 --size-scale 0.85`                        g = 1 (C2's), s = 0.85   lifted
  T3_c2_s2_s05.h5ad       `--form c2-s2 --size-scale 0.5`                         g = 1 (C2's), s = 0.5    lifted
  T3_c2_s2_g15.h5ad       `--form c2-s2 --gata4-exponent 1.5`                     g = 1.5, s = 1
  T3_c2_s2_g05.h5ad       `--form c2-s2 --gata4-exponent 0.5`                     g = 0.5, s = 1
  T3_c2_s2_g2_s067.h5ad   `--form c2-s2 --gata4-exponent 2.0 --size-scale 0.67`   g = 2 (line 3's), s = 0.67   lifted

The size-held part is line 3's c2_gata4_weighted at g, value for value (g = 1 gives C2's factor, combine('c2', 1.0), and
g = 2 T3_c2_s2_g2's); then ln phi is multiplied by s (** s; none at s = 1); nothing rescaled or clipped after. Cells,
streams, composition, seed, count and board stay C2's (the barred s2's). The R76 check runs first, before every other
line's, and only with --captain-ruling R76. Refused (exit 2, nothing built): a command whose --gata4-exponent and
--size-scale are not exactly one of the five strings (the options given and not given count: `--gata4-exponent 1.0
--size-scale 0.85` is not the first string, and line 3's `--form c2-s2 --gata4-exponent 2.0` and F2's `--form c2-s2
--size-scale 0.67` are not row R76's); any form but c2-s2; --composition-from or a --depth-* option beside it; for the
three lifted strings, a --den-ruling other than R76, or none; for the two weight-only strings, any --den-ruling (the row
lifts the 5% check for three strings only); a ledger without exactly one row R76 whose first cell passes
build_outside_form.captain_order_names for the file's string and holds the table: each of the five strings as one whole
backtick-delimited string followed by " (<file>)", as the row pairs them, and one part between semicolons that says
"the 5% check gates nothing for it" and names `--den-ruling R76` and exactly the three lifted strings, each whole. Each
of the five strings under --captain-ruling R22, or any ruling but R76, is refused as before (only 0.67 is a size scale
there, and only 2.0 a weight on c2-s2). The two checks of the options on C2's command, before any cell is drawn and
before anything is written (exit 6, nothing written): the range guard gates every file; the 5% den check gates the two
weight-only files and is recorded, gating nothing, for the three lifted ones. The factor the build applies must be the
very array the checks describe (exit 5, nothing written). The sidecar's combined_build then also holds `r76_line` (the
ruling record, the file, g, s, c and the powers applied, both checks' figures, the lifting if any and the factor's
sha256).

--captain-ruling R79 on --form c2-s2 (row R79, the captain's order of 2026-11-26: T3's refinement round 4 on C2's
command, the T3 REF4 design with its cold challenge's fixes, all five files sent at once; off by default: a command
without --captain-ruling R79 goes on exactly as above, every byte of C2 and of each earlier row's file included, and
nothing below is read or computed): five files on C2's command, each built from one of row R79's option strings, which
one literal table holds (R79_TABLE: the file, its option string as the row pairs them, and whether the row says "the 5%
check gates nothing for it" of it, which it says of all five: each string carries a size scale), each a point of row
R76's family ln phi(g, s) = s c_g (ln phi_M + g ln phi_G),

  T3_c2_s2_s035.h5ad        `--form c2-s2 --size-scale 0.35`                       g = 1 (C2's), s = 0.35   lifted
  T3_c2_s2_g05_s05.h5ad     `--form c2-s2 --gata4-exponent 0.5 --size-scale 0.5`   g = 0.5, s = 0.5         lifted
  T3_c2_s2_gonly_s05.h5ad   `--form c2-s2 --gata4-exponent inf --size-scale 0.5`   g = inf, s = 0.5         lifted
  T3_c2_s2_g2_s05.h5ad      `--form c2-s2 --gata4-exponent 2.0 --size-scale 0.5`   g = 2, s = 0.5           lifted
  T3_c2_s2_g15_s05.h5ad     `--form c2-s2 --gata4-exponent 1.5 --size-scale 0.5`   g = 1.5, s = 0.5         lifted

At a finite g the factor is row R76's r76_factor at (g, s), value for value (line 3's c2_gata4_weighted at g, then
** s), with R76's checks (r76_guards) and loader (r76_loader). At g = inf (the table's one string with "inf", read by
row R79's own parser, which takes exactly those three letters on --gata4-exponent only and tags them gonly, as row R47
names T3_c2_s2_gonly) it is the options' g = inf branch, c2_gonly_factor at (inf, s), value for value (phi_G alone at
C2's held size, then ** s), with gonly_guards and gonly_loader. Nothing rescaled or clipped after. Cells, streams,
composition, seed, count and board stay C2's (the barred s2's). The R79 check runs first, before R76's and every other
line's, and only with --captain-ruling R79. Refused (exit 2, nothing built): a command whose --gata4-exponent and
--size-scale are not exactly one of the five strings (the options given and not given count; the values are read as
every line reads them, as numbers, so `--form c2-s2 --gata4-exponent 1.0 --size-scale 0.35` is not the first string
and R76's `--form c2-s2 --size-scale 0.5` is not R79's); any form but c2-s2; --composition-from or a --depth-* option
beside it; a --den-ruling other than R79, or none; a ledger without exactly one row R79 whose first cell passes
build_outside_form.captain_order_names for the file's string and holds the table: each of the five strings as one
whole backtick-delimited string followed by " (<file>)", as the row pairs them, and exactly one part between
semicolons that says "the 5% check gates nothing for it" and names `--den-ruling R79` and exactly the five strings,
each whole (the row's quoted words say "the 5% check gates nothing." and are not such a part). Each of the five
strings under any other ruling is refused as before (under R22 only 0.67 is a size scale; under R76 only row R76's
strings), and each of row R76's strings under R79. The two checks before any cell is drawn and before anything is
written (exit 6, nothing written): the range guard gates every file; the 5% den check is recorded, gating nothing, for
each. The factor the build applies must be the very array the checks describe (exit 5, nothing written). The
sidecar's combined_build then also holds `r79_line` (the ruling record, the file, g, the string "inf" at g = inf, s, c
and the powers applied, both checks' figures, the lifting and the factor's sha256), and the command line's summary
prints that g as inf. The command line's option help is left as it was, byte for byte: it names the rows up to R76.

The two target-free gates of row R36's file (they build and write nothing):

  g1 DECISIONS.json --wt WT --composition-from allowed-s1
      before any build: the frozen composition routine (bs.composition), given C2's spec with the block replaced, on
      the reference the build draws from; prints the five states' shares over their wild-type shares before
      renormalising, one per line, and nothing else; exit 0 when each is 0.3835131214774091 within a relative 1e-12,
      1 otherwise, 2 refused, 3 not computed
  g2 DECISIONS.json --wt WT --run-dir RUN --board T3:bcatenin --seed 1000001 --composition-from allowed-s1
      before upload: C2's and cD7's gate bars 0.01 + 0.99 cos^2(dp, 1) on each of row N4's three gene sets, as N4's
      diagnostic computes them (experiments/research_2026_11/diagnostics.py t3-gate: its pseudobulk and cos2_ones; all
      genes, the genes feasible down or up, the top fifth by |dp|), each from the file's own build inputs (the cells
      the build draws, in memory; no built file is read); prints the six bars (C2's three, then cD7's) and nothing
      else; exit 0 when on each set cD7's bar is at most C2's plus 0.01, 1 otherwise, 2 refused, 3 not computed
"""

from __future__ import annotations

import copy
import hashlib
import io
import math
import re
import sys
from contextlib import contextmanager, nullcontext, redirect_stdout
from pathlib import Path

import numpy as np
import typer

import vec.provenance  # noqa: F401  (records this build's start, before the method code below is loaded)

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import build_submission as bs  # noqa: E402

from vec import provenance  # noqa: E402

SCRIPT = "experiments/t3_decomp/build_two_knockout.py"
SCRIPT_PATH = Path(__file__).resolve()
BRANCH = "barred"
FROZEN_FACTOR = {"source": "mab21l2", "exponent": 1.0}
TABLE_M, TABLE_G = "transfer_mab21l2_all", "gata4_pooled"
# the two committed factor tables of the T3 run gata4_20261028T122425Z, as the committed decision file records them
# (experiments/t3_decomp/results/gata4/gata4_decisions.json, submissions.parameters.factor_files)
TABLE_SHA256 = {
    TABLE_M: "64b1e0b8572c7193ed1e1f0b7339a26081e601f1e706e6b08193e1718dc5c7e8",
    TABLE_G: "36f2d65a36d2d4fd5307cddfadb56e9724c362c982795f30428eda752b4734ef",
}
# the frozen files' bytes on main at ca2ebc5 (unchanged since the third rehearsal, dc7b6ae), as build_barred_exponent
FROZEN_SHA256 = {
    "t3_decompose.py": "3e671c3bb5ff57bba5ceec80ea36fcc739b1f111399a1e2854400b5c9f5133d2",
    "rules.yaml": "cd2d94aa6e202529a2eead2a7ad60b312d9a2fd9ff0bd989d9320662de750c13",
    "build_submission.py": "1194a63993a27dde4dedfdb9180de035c5f8f2a57aaa7d97b38e7a15ea9e0c74",
}
# base -> (frozen slot, x, the reference sha256 of the base file built on the DGX)
BASES = {
    "s1": ("s1", 1.0, "3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2"),
    "s1x1.5": ("s1", 1.5, "c060fcaa6c7fdc6d6b1a58c8f850a1895a8295abbc50f297de71a32c97f9c420"),
    "s1x2.0": ("s1", 2.0, "7393810fbe3be2cde67cf23ed54f45d22da738d7d3c9e124780a8f6281075661"),
    "s2": ("s2", 1.0, "1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f"),
}
# form kind -> the Gata4 table's weight before c (c is computed for c2 only, 1 otherwise)
KINDS = {"identity": 0.0, "c1": 1.0, "c2": 1.0}
FORMS = tuple([f"identity-{b}" for b in BASES] + ["c1-s1", "c1-s2"] + [f"c2-{b}" for b in BASES])
REFUSED_FORMS = {"c1-s1x1.5": "C1 is fixed on the x = 1 bases only (barred s1, barred s2)",
                 "c1-s1x2.0": "C1 is fixed on the x = 1 bases only (barred s1, barred s2)"}
WHERE = ("the frozen loader's phi_M and phi_G (each already rescaled and clipped to [0.25, 4] by gene_factor) are "
         "combined gene by gene here; t3_decompose.build receives the combination at w = 1 and applies it through "
         "apply_gene_factor after the cells are drawn; nothing is rescaled or clipped after the combination")
LIST_ENTRY = "docs/PLAN.md decision log 2026-11-20 (q); experiments/harness/TRIES.md row R17"
LIST_SEED, LIST_N_CELLS = 1_000_001, "max"
LIST_INPUTS_REFUSAL = "the list fixes seed 1000001 and max_cells"
FACTOR_NOTE = "the frozen loader's record of phi_M only; the factor applied is combined_build.combination"
DEN_NOTE = ("not recorded: the factor files carry gene and factor only (t3_decompose.submission_parameters), so the "
            "den-weighted mean cannot be computed from the tables alone")

# Row R22's second T3 search line (the Gata4-table weight w on C1; decision log 2026-11-23 (ff), (hh)), as row R30
# (the captain's order of 2026-11-24, decision log 2026-11-24 (d)) amends its guard: the 5% limit stops nothing.
TRIES_FILE = HERE.parents[1] / "experiments" / "harness" / "TRIES.md"
G2_ROW = "R22"
G2_EXPONENTS = (2.0,)  # --gata4-exponent: the single value 2.0
G2_FORM = "c1-s2"
G2_OPTION = "`--gata4-exponent 2.0`"
G2_FILE = "T3_c1_s2_g2.h5ad"
G2_LINE = ("T3 search line 2 (R22): the Gata4-table weight w in C1's factor phi_M ** 1 * phi_G ** w on the barred s2; "
           "B (w = 0), C1 (w = 1), this file T3_c1_s2_g2.h5ad (w = 2.0); the line ends at 2.0")
G2_R30 = ("row R30 (the captain's order of 2026-11-24; decision log 2026-11-24 (d)): the limit of 5% on the "
          "den-weighted total gates nothing for this file; the figures are computed and recorded and the file is "
          "built whatever they are")
G2_OLD_LIMIT = 0.05  # the line's old guard limit; recorded only, lifted by row R30, gates nothing
G2_LIST_ENTRY = ("experiments/harness/TRIES.md rows R22 and R30; decision log 2026-11-23 (ff), (hh) and "
                 "2026-11-24 (d)")
G2_DEN_RULE = ("WT E8.75's mean within-cell proportions (build_outside_form.wt_den: load_stage, duplicates dropped); "
               "a den-weighted total is den . phi / sum den")
EXIT_REFUSED, EXIT_FACTOR_MISMATCH = 2, 5

# Row R36 (the captain's named exception of 2026-11-24, "Yes, cD7 on C2, my named exception, bound to the read";
# decision log 2026-11-24 (aa) and (bb)): the file cD7, C2 on the barred s2 with its composition block replaced by the
# decision file's submissions.allowed.s1.composition. Off by default; its own code below, apart from row R22's line.
R36_ROW = "R36"
R36_FORM, R36_BASE = "c2-s2", "s2"
R36_OPTION = "`--composition-from allowed-s1`"
R36_NAMES = (R36_OPTION, "`--form c2-s2`", "`--captain-ruling R36`")  # what row R36's first cell must name
R36_OPENER = "the captain's named exception"
R36_FILE = "T3_c2_s2_cd7.h5ad"
C2_S2_SHA256 = "47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7"  # C2's file (t3two_fb01487)
# --composition-from -> submissions.<branch>.<slot>; only that entry's composition block is copied, never its factor
# (gata4) or its depth step (depth: true; a closed lever): the factor stays C2's and C2 has no depth step
COMPOSITION_FROM = {"allowed-s1": ("allowed", "s1")}
# experiments/t3_decomp/results/gata4/gata4_decisions.json, the run of 2026-10-28 that the frozen rules decided on
R36_DECISION_SHA256 = "9440ae67625608112666d5235e6ceaa5dd9f47a4b0bbd30e16943f32f35c5a8c"
R36_LOG2_MAGNITUDE = 1.3826521560990648  # rule D7's m (the released validation knockout), the block's log2_magnitude
R36_RATIO = 2.0 ** -R36_LOG2_MAGNITUDE  # 0.3835131214774091: each scenario state's share over its wild-type share
R36_RATIO_RTOL = 1e-12  # G1: each of the five ratios within this relative distance of R36_RATIO
R36_N_STATES = 5
R36_G2_MARGIN = 0.01  # G2: on each of N4's three gene sets, cD7's bar at most C2's plus this
N4_SETS = ("all_genes", "feasible_down_or_up", "top20pct_abs_dp")  # row N4's gene sets, in its order
G2_FAILED_LINE = "G2 failed: no upload; C2 stays."  # the launcher's fixed line when G2 does not pass
R36_SOURCE = ("decision file submissions.allowed.s1.composition (--composition-from allowed-s1, row R36), replacing "
              "the barred s2's composition block")
R36_LINE = ("T3 file cD7 (row R36): C2 on the barred s2 (--form c2-s2) rebuilt by its own command with one change, "
            "its composition block replaced by the decision file's submissions.allowed.s1.composition (the five "
            "scenario states at 2^-1.3826521560990648 = 0.3835 of their wild-type share before renormalising); C2's "
            "factor, row number (max_cells, read by code), seed, streams and coordinates rule unchanged")
R36_LIST_ENTRY = "experiments/harness/TRIES.md row R36; decision log 2026-11-24 (aa) and (bb)"
N4_SCRIPT = HERE.parents[1] / "experiments" / "research_2026_11" / "diagnostics.py"  # row N4's t3-gate
EXIT_GATE_FAILED, EXIT_GATE_CRASHED = 1, 3  # g1 and g2: the gate does not pass; it could not be computed

# Row R22's third T3 search line (decision log 2026-11-24 (kk); the captain's own yes naming T3_c2_s2_g2.h5ad, ruling
# row R41): the gata4_pooled weight g inside C2's size-held product on the barred s2. Off by default; its own code
# below, beside combine(), which it leaves unchanged.
C2G_FORM = "c2-s2"
C2G_WEIGHTS = (2.0,)  # --gata4-exponent on c2-s2: the weight g, the single value 2.0
C2G_OPTION = "`--form c2-s2 --gata4-exponent 2.0`"  # what row R22's first cell must name for this line
C2G_FILE = "T3_c2_s2_g2.h5ad"
C2G_MAX_DEN_SHIFT = 0.05  # the guard: the file's den-weighted total within 5% of C2's factor's
C2G_LINE = ("T3 search line 3 (R22): the gata4_pooled weight g in C2's size-held product phi_M ** c * phi_G ** (g c), "
            "c = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), on the barred s2; B (g = 0), C2 (g = 1), this file "
            "T3_c2_s2_g2.h5ad (g = 2); the line ends at 2")
C2G_LIST_ENTRY = ("experiments/harness/TRIES.md rows R22 and R41; decision log 2026-11-24 (w), (cc) and (kk)")
C2G_GUARD_RULE = ("the den-weighted total of the file's factor lies within 5% of C2's factor's (den from released WT "
                  "E8.75 as build_outside_form.wt_den makes it, the two pinned tables only); C2's shift and the file's "
                  "shift from phi_M's total are recorded and gate nothing")
EXIT_GUARD_FAILED = 6

# Two options on C2's command (the T3 design of the design round of 2026-11-24, section 2, as its cold challenge fixes
# it: files F1, F2 and F3), each only with --captain-ruling R22 and only while row R22's first cell names the file's
# whole option string: on --form c2-s2, --gata4-exponent also takes inf (g = inf: the gata4_pooled table alone at C2's
# held size, phi_M dropped) and --size-scale takes 0.67 (ln phi multiplied by s after the size-held combination). The
# range guard and line 3's 5% den check gate each file; a ruling row named by --den-ruling lifts the 5% check for the
# files whose option strings it names. Off by default; its own code below, beside combine() and c2_gata4_weighted(),
# which it leaves unchanged.
GONLY_FORM = "c2-s2"
GONLY_SCALES = (0.67,)  # --size-scale: the single value 0.67
GONLY_POINTS = {"F1": (math.inf, 1.0), "F2": (1.0, 0.67), "F3": (math.inf, 0.67)}  # file -> (g, s)
# file -> its option string, which a ruling row's first cell must hold as one whole backtick-delimited string (F1's
# string sits inside F3's, so containment would let F3's entry open F1; the whole-string match does not)
GONLY_OPTIONS = {"F1": "--form c2-s2 --gata4-exponent inf", "F2": "--form c2-s2 --size-scale 0.67",
                 "F3": "--form c2-s2 --gata4-exponent inf --size-scale 0.67"}
GONLY_FILES = {"F1": "T3_c2_s2_gonly.h5ad", "F2": "T3_c2_s2_s067.h5ad", "F3": "T3_c2_s2_gonly_s067.h5ad"}
GONLY_RANGE_FACTOR = 2.0  # the range guard: max at most twice C2's max, min at least C2's min over this
GONLY_LIFT_WORDS = "the 5% check gates nothing for it"  # what a lifting row's first cell must say, with the string
GONLY_LIFT_ROW = re.compile(r"R[0-9]+")  # --den-ruling: a ruling row's id
GONLY_LINE = ("T3 options on C2's command (row R22): the family ln phi(g, s) = s c_g (ln phi_M + g ln phi_G), c_g = "
              "SD(ln phi_M) / SD(ln phi_M + g ln phi_G), on the barred s2; F1 T3_c2_s2_gonly.h5ad (g = inf: ln phi = "
              "(SD(ln phi_M) / SD(ln phi_G)) ln phi_G, phi_M dropped), F2 T3_c2_s2_s067.h5ad (g = 1, s = 0.67), F3 "
              "T3_c2_s2_gonly_s067.h5ad (g = inf, s = 0.67)")
GONLY_LIST_ENTRY = ("experiments/harness/TRIES.md row R22 (the captain's order naming `--form c2-s2 --gata4-exponent "
                    "inf`, `--form c2-s2 --size-scale 0.67` and `--form c2-s2 --gata4-exponent inf --size-scale "
                    "0.67`) and any row named by --den-ruling; the T3 design of the design round of 2026-11-24, "
                    "section 2, and its cold challenge")
GONLY_RANGE_RULE = ("the applied factor's max at most twice C2's factor's max and its min at least half C2's min "
                    "(C2's factor combine('c2', 1.0) from the same two pinned tables)")
GONLY_DEN_RULE = ("line 3's check: the den-weighted total of the file's factor lies within 5% of C2's factor's (den "
                  "from released WT E8.75 as build_outside_form.wt_den makes it, the two pinned tables only); it "
                  "gates the file unless the row named by --den-ruling names the file's option string and says \"the "
                  "5% check gates nothing for it\"; the shifts from phi_M's total are recorded and gate nothing")

Refused = bs.Refused


class FactorMismatch(Refused):
    """The factor the build applies is not the array the line's figures were computed from (exit status 5)."""


class GuardFailed(Refused):
    """Row R22's third line: the file's den-weighted total is not within 5% of C2's (exit status 6, nothing written)."""


class GonlyGuardFailed(GuardFailed):
    """The options on C2's command (F1, F2, F3): the range guard or the 5% den check (not lifted for the file) fails
    (exit status 6, nothing written)."""


class DepthCheckFailed(Refused):
    """The depth line (dd, ddh, dg): inside the build, the counts the depth step receives are not C2's (rows, labels,
    sources or counts), or a cell's ratio is not its label's (exit status 7, nothing written)."""


app = typer.Typer(no_args_is_help=True, add_completion=False)


def parse_form(form) -> tuple[str, str]:
    """--form: exactly one of FORMS; returns (kind, base), else Refused."""
    f = str(form) if form is not None else ""
    if f in REFUSED_FORMS:
        raise Refused(f"--form {form!r} is refused: {REFUSED_FORMS[f]}")
    if f not in FORMS:
        raise Refused(f"--form {form!r} is not one of {list(FORMS)}")
    kind, base = f.split("-", 1)
    return kind, base


def frozen_paths() -> dict[str, Path]:
    """The three frozen files as the frozen builder itself names them."""
    return {"t3_decompose.py": bs.T3D_PATH, "rules.yaml": bs.RULES_PATH, "build_submission.py": bs.BUILDER_PATH}


def check_frozen(paths: dict[str, Path] | None = None) -> dict[str, str]:
    """The frozen files' sha256 on disk against the pinned bytes; returns them."""
    paths = frozen_paths() if paths is None else paths
    got = {k: bs.sha256(p) for k, p in paths.items()}
    bad = [k for k in FROZEN_SHA256 if got.get(k) != FROZEN_SHA256[k]]
    if bad:
        raise Refused(f"frozen file(s) {bad} differ from the pinned bytes (main at ca2ebc5): not built")
    return got


def check_clean() -> dict:
    """The tree by vec.provenance's rule (the one the build stamp records); refuses unless clean."""
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    if commit == "unknown" or unclean is None:
        raise Refused("git cannot tell the tree's state: not built (the stamp would not say clean)")
    if unclean:
        raise Refused(f"the tree is not clean ({len(unclean)} path(s), e.g. {unclean[:3]}): commit first (the build "
                      "stamp must say clean)")
    return {"commit": commit, "clean": True}


def _sd(v: np.ndarray) -> float:
    return float(np.std(v, ddof=0))


def _shift(phi: np.ndarray) -> dict:
    """The common shift of a factor, from the factor alone: gene mean of ln phi and cos^2(ln phi, 1)."""
    lp = np.log(phi)
    ss = float(lp @ lp)
    return {"gene_mean_ln": float(lp.mean()), "cos2_ln_with_ones": float(lp.sum() ** 2 / (lp.size * ss)) if ss > 0
            else 0.0, "sd_ln": _sd(lp)}


def combine(kind: str, x: float, phi_m: np.ndarray, phi_g: np.ndarray) -> tuple[np.ndarray, dict]:
    """The form's factor from the two tables as loaded, gene by gene; returns (phi, record).

    identity: phi_M ** x (phi_M itself at x = 1); c1: phi_M ** x * phi_G; c2: phi_M ** (x c) * phi_G ** c with
    c = SD(x ln phi_M) / SD(x ln phi_M + ln phi_G). `**` is the operation the frozen build applies for w = x."""
    if kind not in KINDS:
        raise Refused(f"unknown form kind {kind!r}")
    phi_m = np.asarray(phi_m, dtype=np.float64)
    phi_g = np.asarray(phi_g, dtype=np.float64)
    if phi_m.shape != phi_g.shape or phi_m.ndim != 1:
        raise Refused("the two factor tables do not cover the same genes")
    for nm, p in ((TABLE_M, phi_m), (TABLE_G, phi_g)):
        if not np.all(np.isfinite(p) & (p > 0)):
            raise Refused(f"{nm}: factor values must be finite and positive")
    base = phi_m if x == 1.0 else phi_m ** x
    g = KINDS[kind]
    lm, lg = np.log(phi_m), np.log(phi_g)
    sd_base, sd_sum = _sd(x * lm), _sd(x * lm + lg)
    c = 1.0
    if kind == "identity":
        phi = base
    elif kind == "c1":
        phi = base * phi_g
    else:
        if not sd_sum > 0:
            raise Refused("SD(x ln phi_M + ln phi_G) is 0: c is undefined")
        c = sd_base / sd_sum
        phi = (phi_m ** (x * c)) * (phi_g ** c)
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused("the combined factor is not finite and positive")
    sd_m, sd_g = _sd(lm), _sd(lg)
    rec = {"kind": kind, "x": float(x), "gata4_weight": g * c, "mab21l2_exponent": float(x) * c, "c": float(c),
           "c_rule": "SD(x ln phi_M) / SD(x ln phi_M + ln phi_G), population SD over the genes" if kind == "c2"
           else "c = 1",
           "sd_ln_phi_m": sd_m, "sd_ln_phi_g": sd_g, "sd_ln_base": sd_base, "sd_ln_base_plus_g": sd_sum,
           "corr_ln_m_ln_g": float(np.corrcoef(lm, lg)[0, 1]) if sd_m > 0 and sd_g > 0 else None,
           "base_factor": _shift(base), "factor": _shift(phi), "den_weighted_mean": DEN_NOTE,
           "min": float(phi.min()), "max": float(phi.max()), "n_up": int((phi > 1).sum()),
           "n_down": int((phi < 1).sum())}
    return phi, rec


@contextmanager
def combined_loader(kind: str, x: float, out: dict):
    """bs.load_factor, for this process only: the frozen loader for phi_M (the barred spec's table) and for phi_G,
    each checked against TABLE_SHA256 too, then `combine`. `out` receives both tables' records and the combination."""
    original = bs.load_factor

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        phi_m, rec_m = original(run_dir, TABLE_M, genes, dec)
        phi_g, rec_g = original(run_dir, TABLE_G, genes, dec)
        for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
            if rec["sha256"] != TABLE_SHA256[nm]:
                raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                              f"({TABLE_SHA256[nm][:16]}): not built")
        phi, crec = combine(kind, x, phi_m, phi_g)
        out.update(tables={TABLE_M: rec_m, TABLE_G: rec_g}, combination=crec, n_loads=out.get("n_loads", 0) + 1)
        return phi, rec_m

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


# ------------------------------------------------------------------ row R22's second T3 line (off by default)


def _outside_form():
    """build_outside_form, imported only on the R22 line (its captain_order_names and wt_den are reused, not copied)."""
    import build_outside_form

    return build_outside_form


def check_gata4(form, gata4_exponent=None, captain_ruling=None, tries: Path | None = None) -> dict | None:
    """--gata4-exponent and --captain-ruling: None when neither is given (nothing is read); else the line's record,
    or Refused: any other ruling, the ruling without the option, any value but 2.0, any form but c1-s2, the option
    without the ruling, or a ledger whose row R22 does not pass build_outside_form.captain_order_names."""
    if gata4_exponent is None and captain_ruling is None:
        return None
    if captain_ruling is not None and captain_ruling != G2_ROW:
        raise Refused(f"--captain-ruling {captain_ruling!r} is refused: only row {G2_ROW} orders a file of this "
                      "builder")
    if gata4_exponent is None:
        raise Refused(f"--captain-ruling {G2_ROW} is taken only with --gata4-exponent 2.0 on --form {G2_FORM}")
    try:
        w = float(gata4_exponent)
    except (TypeError, ValueError):
        raise Refused(f"--gata4-exponent {gata4_exponent!r} is refused: only 2.0") from None
    if isinstance(gata4_exponent, bool) or w not in G2_EXPONENTS:
        raise Refused(f"--gata4-exponent {gata4_exponent!r} is refused: only 2.0")
    if form != G2_FORM:
        raise Refused(f"--gata4-exponent 2.0 is taken on --form {G2_FORM} only, not {form!r}")
    if captain_ruling is None:
        raise Refused(f"--gata4-exponent 2.0 is taken only with --captain-ruling {G2_ROW}")
    p = TRIES_FILE if tries is None else Path(tries)
    if not _outside_form().captain_order_names(G2_ROW, G2_OPTION, p):
        raise Refused(f"--gata4-exponent 2.0 --captain-ruling {G2_ROW} is refused: experiments/harness/TRIES.md has "
                      f"no single row {G2_ROW} whose first cell is the captain's order (\"ruling, no run\") naming "
                      f"{G2_OPTION} and `--captain-ruling {G2_ROW}`")
    return {"row": G2_ROW, "option": G2_OPTION, "gata4_exponent": w, "ledger": "experiments/harness/TRIES.md",
            "ledger_sha256": bs.sha256(p)}


def factor_sha256(phi) -> str:
    """The sha256 of a factor's float64 bytes, as build_outside_form's guard records it."""
    return hashlib.sha256(np.ascontiguousarray(phi, dtype=np.float64).tobytes()).hexdigest()


def c1_gata4_weighted(phi_m, phi_g, w: float) -> tuple[np.ndarray, dict]:
    """C1 on an x = 1 base with the Gata4 weight w: phi_M ** 1 * phi_G ** w (phi_M ** 1 is phi_M itself, as combine's
    base at x = 1), nothing rescaled or clipped after the product; returns (phi, record)."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    phi_g = np.asarray(phi_g, dtype=np.float64)
    if phi_m.shape != phi_g.shape or phi_m.ndim != 1:
        raise Refused("the two factor tables do not cover the same genes")
    for nm, p in ((TABLE_M, phi_m), (TABLE_G, phi_g)):
        if not np.all(np.isfinite(p) & (p > 0)):
            raise Refused(f"{nm}: factor values must be finite and positive")
    w = float(w)
    phi = phi_m * (phi_g ** w)
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused("the combined factor is not finite and positive")
    lm, lg = np.log(phi_m), np.log(phi_g)
    sd_m, sd_g = _sd(lm), _sd(lg)
    rec = {"kind": "c1", "x": 1.0, "gata4_weight": w, "mab21l2_exponent": 1.0, "c": 1.0,
           "c_rule": "c = 1; the Gata4 weight is the line's w", "sd_ln_phi_m": sd_m, "sd_ln_phi_g": sd_g,
           "sd_ln_base": _sd(lm), "sd_ln_base_plus_wg": _sd(lm + w * lg),
           "corr_ln_m_ln_g": float(np.corrcoef(lm, lg)[0, 1]) if sd_m > 0 and sd_g > 0 else None,
           "base_factor": _shift(phi_m), "factor": _shift(phi),
           "den_weighted_mean": "recorded in combined_build.gata4_line.figures",
           "min": float(phi.min()), "max": float(phi.max()), "n_up": int((phi > 1).sum()),
           "n_down": int((phi < 1).sum())}
    return phi, rec


def gata4_figures(phi_m, phi_g, w: float, den) -> dict:
    """The line's figures, a pure function of the two tables and den: the den-weighted totals of phi_M, of
    phi_M * phi_G (C1's) and of phi_M * phi_G ** w, and the relative shifts of the last two against phi_M's. They
    gate nothing (row R30); the old limit is recorded as lifted. Also the sha256 of the factor they describe."""
    phi, _ = c1_gata4_weighted(phi_m, phi_g, w)
    phi_m = np.asarray(phi_m, dtype=np.float64)
    phi_g = np.asarray(phi_g, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    if den.shape != phi_m.shape or not np.all(np.isfinite(den) & (den >= 0)) or not den.sum() > 0:
        raise Refused("den does not cover the tables' genes with finite, non-negative weights and a positive sum")
    s = float(den.sum())
    t_m, t_c1, t_w = (float(den @ f / s) for f in (phi_m, phi_m * phi_g, phi))
    return {"gata4_exponent": float(w), "den_total_phi_m": t_m, "den_total_c1": t_c1, "den_total_factor": t_w,
            "den_shift_c1": t_c1 / t_m - 1.0, "den_shift_factor": t_w / t_m - 1.0,
            "old_limit": {"max_abs_den_shift": G2_OLD_LIMIT, "status": "lifted by row R30: gates nothing"},
            "gates": None, "factor_sha256": factor_sha256(phi)}


def run_gata4_figures(dec: dict, wt: Path, run_dir: Path | None, w: float) -> dict:
    """The figures from released WT E8.75 and the two pinned tables only, before any build: den and genes by
    build_outside_form.wt_den, both tables through the frozen loader in that gene order (checked against
    TABLE_SHA256). Never refuses on a figure."""
    genes, den = _outside_form().wt_den(Path(wt))
    phi_m, rec_m = bs.load_factor(run_dir, TABLE_M, genes, dec)
    phi_g, rec_g = bs.load_factor(run_dir, TABLE_G, genes, dec)
    for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
        if rec["sha256"] != TABLE_SHA256[nm]:
            raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                          f"({TABLE_SHA256[nm][:16]}): not built")
    den = np.asarray(den, dtype=np.float64)
    g = gata4_figures(phi_m, phi_g, w, den)
    g |= {"wt": str(wt), "wt_sha256": bs.sha256(Path(wt)), "den_rule": G2_DEN_RULE, "den": den.tolist(),
          "den_sha256": factor_sha256(den), "n_genes": len(genes),
          "genes_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest(),
          "tables_sha256": {TABLE_M: rec_m["sha256"], TABLE_G: rec_g["sha256"]}, "genes": list(genes)}
    return g


@contextmanager
def gata4_loader(w: float, out: dict, figures: dict):
    """bs.load_factor, for this process only, on the R22 line: the frozen loader for phi_M and phi_G (checked against
    TABLE_SHA256 too), then c1_gata4_weighted at w, the hook C1's factor goes through. FactorMismatch (exit 5) unless
    the build's genes are the figures' and the factor's sha256 is the one the figures were computed from."""
    original = bs.load_factor

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        if list(genes) != figures["genes"]:
            raise FactorMismatch("the build's genes are not those the line's figures were computed over: not built")
        phi_m, rec_m = original(run_dir, TABLE_M, genes, dec)
        phi_g, rec_g = original(run_dir, TABLE_G, genes, dec)
        for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
            if rec["sha256"] != TABLE_SHA256[nm]:
                raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                              f"({TABLE_SHA256[nm][:16]}): not built")
        phi, crec = c1_gata4_weighted(phi_m, phi_g, w)
        sha = factor_sha256(phi)
        if sha != figures["factor_sha256"]:
            raise FactorMismatch(f"the factor the build would apply ({sha[:16]}) is not the one the line's figures "
                                 f"were computed from ({figures['factor_sha256'][:16]}): not built")
        out.update(tables={TABLE_M: rec_m, TABLE_G: rec_g}, combination=crec, factor_sha256=sha,
                   n_loads=out.get("n_loads", 0) + 1)
        return phi, rec_m

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


# ------------------------------------------------------------------ row R36: the file cD7 (off by default)


def ruling_row_names(row: str, opener: str, names, path: Path | None = None) -> bool:
    """Whether the ledger (TRIES_FILE) holds exactly one row `row` whose first cell starts with "<row> (<opener>",
    says "ruling, no run" and names every one of `names` (each given with its backticks): the check
    build_outside_form.captain_order_names makes of a captain's order, for a ruling row that opens otherwise (row R36
    opens "the captain's named exception"). A row that names them only outside its first cell does not pass."""
    p = TRIES_FILE if path is None else Path(path)
    if not p.is_file():
        return False
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {row} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    return first.startswith(f"{row} ({opener}") and "ruling, no run" in first and all(n in first for n in names)


def check_composition_from(form, composition_from=None, captain_ruling=None, gata4_exponent=None,
                           tries: Path | None = None) -> dict | None:
    """--composition-from and --captain-ruling R36: None when a command has neither the option nor the ruling R36
    (nothing is read; every such command goes on to check_gata4 as before); else the record, or Refused: the ruling
    R36 without the option, any value but allowed-s1, any form but c2-s2, the option without --captain-ruling R36,
    --gata4-exponent beside it, or a ledger whose row R36 does not pass ruling_row_names for R36_NAMES."""
    if composition_from is None and captain_ruling != R36_ROW:
        return None
    if composition_from is None:
        raise Refused(f"--captain-ruling {R36_ROW} is taken only with --composition-from allowed-s1 on --form "
                      f"{R36_FORM}")
    if not isinstance(composition_from, str) or composition_from not in COMPOSITION_FROM:
        raise Refused(f"--composition-from {composition_from!r} is refused: only allowed-s1")
    if form != R36_FORM:
        raise Refused(f"--composition-from allowed-s1 is taken on --form {R36_FORM} only, not {form!r}")
    if captain_ruling != R36_ROW:
        raise Refused(f"--composition-from allowed-s1 is taken only with --captain-ruling {R36_ROW}, not "
                      f"{captain_ruling!r}")
    if gata4_exponent is not None:
        raise Refused("--composition-from allowed-s1 takes no --gata4-exponent (cD7 is C2 with one change)")
    p = TRIES_FILE if tries is None else Path(tries)
    if not ruling_row_names(R36_ROW, R36_OPENER, R36_NAMES, p):
        raise Refused(f"--composition-from allowed-s1 --captain-ruling {R36_ROW} is refused: "
                      f"experiments/harness/TRIES.md has no single row {R36_ROW} whose first cell is the captain's "
                      f"named exception (\"ruling, no run\") naming {', '.join(R36_NAMES)}")
    return {"row": R36_ROW, "option": R36_OPTION, "composition_from": composition_from,
            "ledger": "experiments/harness/TRIES.md", "ledger_sha256": bs.sha256(p)}


def composition_block(dec: dict, rules: dict, decision: Path, composition_from: str = "allowed-s1") -> tuple[dict,
                                                                                                             dict]:
    """The block --composition-from names, read by code from the decision file (allowed-s1:
    submissions.allowed.s1.composition), and that block only: the entry's factor and depth step are never read.
    Refused unless the file is the pinned one (R36_DECISION_SHA256); the block's depleted and spared states are
    rules.yaml's scenario states, in order (five depleted); its log2 magnitude is rule D7's m, as the decision file
    records it (rules.D7.m) and as pinned here (R36_LOG2_MAGNITUDE); and it is a central composition whose depletion
    factor is 2 ** -magnitude. Returns (a copy of the block, its record)."""
    if not isinstance(composition_from, str) or composition_from not in COMPOSITION_FROM:
        raise Refused(f"--composition-from {composition_from!r} is refused: only allowed-s1")
    branch, slot = COMPOSITION_FROM[composition_from]
    where = f"submissions.{branch}.{slot}.composition"
    got = bs.sha256(Path(decision))
    if got != R36_DECISION_SHA256:
        raise Refused(f"the decision file's sha256 {got[:16]} is not the pinned one ({R36_DECISION_SHA256[:16]}): "
                      "not built")
    try:
        block = dec["submissions"][branch][slot]["composition"]
    except (KeyError, TypeError):
        raise Refused(f"the decision file has no {where}: not built") from None
    if not isinstance(block, dict):
        raise Refused(f"{where} is not a composition block: not built")
    scen = rules["presets"]["gata4"]["scenario"]
    if (block.get("depleted_states") != list(scen["depleted_states"]) or len(scen["depleted_states"]) != R36_N_STATES
            or block.get("spared_states") != list(scen["spared_states"])):
        raise Refused(f"the states of {where} are not rules.yaml's scenario states: not built")
    d7 = (dec.get("rules") or {}).get("D7")
    m = d7.get("m") if isinstance(d7, dict) else None
    l2 = block.get("log2_magnitude")
    if not (type(m) is float and m == R36_LOG2_MAGNITUDE and type(l2) is float and l2 == m):
        raise Refused(f"the log2 magnitude of {where} is not rule D7's m ({R36_LOG2_MAGNITUDE!r}): not built")
    if block.get("kind") != "central" or block.get("depletion_factor") != R36_RATIO:
        raise Refused(f"{where} is not the central composition at depletion factor {R36_RATIO!r}: not built")
    return copy.deepcopy(block), {"source": f"decision file {where}", "decision": str(decision),
                                  "decision_sha256": got, "pinned_sha256": R36_DECISION_SHA256,
                                  "rule_d7_m": m, "taken": "the composition block only (not the entry's factor or "
                                  "depth step)"}


def with_composition(spec: dict, block: dict) -> dict:
    """`spec` with its composition block replaced by a copy of `block`; every other key as it was."""
    return {**spec, "composition": copy.deepcopy(block)}


def composition_ratios(V, f, info: dict, states) -> list[float]:
    """Each of `states`' share in the composition f over its wild-type share V.p, before the routine renormalised:
    f / p on the state divided by f / p on the states the routine left alone (whose one common value is the
    routine's renormalising factor). nan for a state the routine did not find or whose wild-type share is not
    positive, and for every state when the states left alone do not share one value (the routine changed others)."""
    p = np.asarray(V.p, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    nan = [float("nan")] * len(states)
    pos = {t: i for i, t in enumerate(V.labels)}
    found = set(info.get("states_found") or []) - set(info.get("states_missing") or [])
    named = {pos[s] for s in states if s in pos}
    rest = np.array([i for i in range(p.size) if i not in named and p[i] > 0], dtype=int)
    if f.shape != p.shape or not rest.size:
        return nan
    z = f[rest] / p[rest]
    if not (np.all(np.isfinite(z)) and z.min() > 0 and z.max() - z.min() <= 1e-12 * z.max()):
        return nan
    zc = float(z.mean())
    return [float(f[pos[s]] / p[pos[s]] / zc) if s in found and s in pos and p[pos[s]] > 0 else float("nan")
            for s in states]


def ratios_pass(ratios) -> bool:
    """G1's verdict: five ratios, each finite and within R36_RATIO_RTOL of R36_RATIO (0.3835)."""
    r = np.asarray(ratios, dtype=np.float64)
    return (r.shape == (R36_N_STATES,) and bool(np.all(np.isfinite(r)))
            and bool(np.all(np.abs(r - R36_RATIO) <= R36_RATIO_RTOL * R36_RATIO)))


@contextmanager
def r36_composition(block: dict, out: dict):
    """bs.resolve_spec and bs.composition, for this process only, on row R36's line: the frozen resolve_spec, then
    the barred s2's composition block replaced by `block` (any other slot is refused); the frozen composition routine,
    then G1's check of its result (Refused unless the five states are at R36_RATIO of their wild-type share before
    renormalising), before any cell is drawn. `out` receives the block replaced, the ratios and the call counts."""
    orig_resolve, orig_comp = bs.resolve_spec, bs.composition
    states = list(block["depleted_states"])

    def resolve(dec, rules, branch, submission, hashes):
        spec, info = orig_resolve(dec, rules, branch, submission, hashes)
        if (branch, submission) != (BRANCH, R36_BASE):
            raise Refused(f"--composition-from replaces the barred s2's composition only, not {branch} {submission}")
        if set(spec) != {"factor", "composition"} or spec["factor"] != FROZEN_FACTOR:
            raise Refused("the barred s2's spec is not C2's (its factor alone, no depth, stage or coverage step): "
                          "not built")
        out.update(replaced=copy.deepcopy(spec.get("composition")), n_resolves=out.get("n_resolves", 0) + 1)
        return with_composition(spec, block), {**info, "composition_source": R36_SOURCE}

    def comp(V, spec, params):
        f, info = orig_comp(V, spec, params)
        r = composition_ratios(V, f, info, states)
        if not ratios_pass(r):
            raise Refused(f"the composition routine does not give the {R36_N_STATES} scenario states at "
                          f"{R36_RATIO:.4f} of their wild-type share (G1): not built")
        out.update(ratios=dict(zip(states, r, strict=True)), n_compositions=out.get("n_compositions", 0) + 1)
        return f, info

    bs.resolve_spec, bs.composition = resolve, comp
    try:
        yield
    finally:
        bs.resolve_spec, bs.composition = orig_resolve, orig_comp


def run_g1(decision: Path, wt: Path, composition_from: str = "allowed-s1", threads: int = 16) -> list[float]:
    """G1, before any build: the frozen builder's composition routine (bs.composition), given C2's spec (the barred
    s2's, as bs.barred_table gives it) with its composition block replaced as the build replaces it, on the reference
    the build draws from (bs.load_reference on `wt`, the frozen builder's preset); returns the five states' ratios to
    their wild-type shares before renormalising (composition_ratios), in the block's order. Nothing is built or
    written."""
    check_frozen()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)
    block, _ = composition_block(dec, rules, Path(decision), composition_from)
    spec = with_composition(bs.barred_table(dec, rules)[R36_BASE], block)
    bs.validate_spec(spec, rules)
    ref = bs.load_reference(Path(wt), dec["preset"], rules, None, threads)
    f, info = bs.composition(ref.V, spec, {})
    return composition_ratios(ref.V, f, info, block["depleted_states"])


def _n4():
    """Row N4's diagnostic (experiments/research_2026_11/diagnostics.py), imported by G2 only: its pseudobulk and
    cos2_ones are used as they are, not copied."""
    import importlib.util

    name = "_n4_diagnostics"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, N4_SCRIPT)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        try:
            spec.loader.exec_module(mod)
        except BaseException:
            sys.modules.pop(name, None)
            raise
    return sys.modules[name]


def n4_sets(dp, pb_wt) -> dict[str, np.ndarray]:
    """Row N4's three gene sets for one file's dp, as its t3-gate builds them for each file (N4_SETS): every gene;
    the genes feasible down (pb_wt >= 0.25) or with dp > 0; the fifth of the genes with the largest |dp|."""
    dp = np.asarray(dp, dtype=np.float64)
    feas = np.asarray(pb_wt, dtype=np.float64) >= 0.25
    top = np.argsort(-np.abs(dp))[: max(1, dp.size // 5)]
    sets = {"all_genes": np.arange(dp.size), "feasible_down_or_up": np.flatnonzero(feas | (dp > 0)),
            "top20pct_abs_dp": top}
    return {s: sets[s] for s in N4_SETS}


def n4_bars(dp, pb_wt) -> dict[str, float]:
    """Row N4's gate bar 0.01 + 0.99 cos^2(dp, 1) on each of its three gene sets (n4_sets), cos2_ones its own."""
    dp = np.asarray(dp, dtype=np.float64)
    return {s: 0.01 + 0.99 * _n4().cos2_ones(dp[idx]) for s, idx in n4_sets(dp, pb_wt).items()}


def g2_passes(bars_c2: dict, bars_cd7: dict) -> bool:
    """G2's verdict: on each of N4's three gene sets both bars are finite and cD7's is at most C2's plus
    R36_G2_MARGIN."""
    return all(s in bars_c2 and s in bars_cd7 and bool(np.isfinite(bars_c2[s]) and np.isfinite(bars_cd7[s]))
               and bars_cd7[s] <= bars_c2[s] + R36_G2_MARGIN for s in N4_SETS)


def run_g2(decision: Path, wt: Path, run_dir: Path | None, board, n_cells: int, seed: int = 1_000_001,
           composition_from: str = "allowed-s1", threads: int = 16) -> dict:
    """G2, before upload: C2's and cD7's gate bars on N4's three gene sets (n4_bars), as row N4's diagnostic computes
    them: dp is the pseudobulk of the drawn cells over the board's panel minus that of `wt` as released
    (vec.prep.load_stage, as the scorer loads the reference), both by N4's pseudobulk. Each file's cells come from its
    own build inputs, through the build's own code up to the draw: bs.preflight on the barred s2 (cD7's under
    r36_composition), bs.load_reference, the run records check, bs.inputs_for with C2's factor (combined_loader, c2 at
    x = 1) and bs.construct at the seed and the stream "submission:barred:s2". Nothing is assembled, written or read
    from a built file. Returns both files' three bars and the verdict (g2_passes)."""
    from vec.prep import load_stage

    n4 = _n4()
    check_frozen()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)
    block, _ = composition_block(dec, rules, Path(decision), composition_from)
    kind, base = parse_form(R36_FORM)
    slot, x, _ = BASES[base]
    ref, dps, hook = None, {}, {}
    for name in ("c2", "cd7"):
        got: dict = {}
        hooked = r36_composition(block, hook) if name == "cd7" else nullcontext()
        with hooked, combined_loader(kind, x, got):
            pf = bs.preflight(Path(decision), BRANCH, slot, board, n_cells)
            if ref is None:
                ref = bs.load_reference(Path(wt), pf["preset"], rules, None, threads)
                bs.check_inputs_against_run(bs.run_records(run_dir, pf["preset"]), pf["preset"],
                                            ref.record["sha256"], {})
                panel = list(board.genes)
                gpos = {g: i for i, g in enumerate(ref.genes)}
                missing = [g for g in panel if g not in gpos]
                if missing:
                    raise Refused(f"{len(missing)} board genes are not in the WT reference")
                cols = [gpos[g] for g in panel]
                wt_rel = load_stage(Path(wt))
                wpos = {g: i for i, g in enumerate(wt_rel.var_names)}
                pb_wt = n4.pseudobulk(wt_rel.X.tocsr()[:, [wpos[g] for g in panel]])
            kw, _ = bs.inputs_for(pf["spec"], pf["dec"], ref, run_dir)
            built = bs.construct(ref, n_cells, bs.t3d.stable_rng(seed, f"submission:{BRANCH}:{slot}"), **kw)
        if got.get("n_loads") != 1 or kw["w"] != 1.0:
            raise RuntimeError("the draw did not apply C2's combined factor once at w = 1")
        dps[name] = n4.pseudobulk(built.X[:, cols]) - pb_wt
        del built
    if hook.get("n_resolves") != 1 or hook.get("n_compositions") != 1:
        raise RuntimeError("cD7's composition was not replaced and checked exactly once")
    bars = {name: n4_bars(dp, pb_wt) for name, dp in dps.items()}
    return {"c2": bars["c2"], "cd7": bars["cd7"], "margin": R36_G2_MARGIN, "sets": list(N4_SETS),
            "passes": g2_passes(bars["c2"], bars["cd7"])}


# ------------------------------------------------------------- row R22's third T3 line: T3_c2_s2_g2 (off by default)


def check_c2_gata4(form, gata4_exponent=None, captain_ruling=None, tries: Path | None = None) -> dict | None:
    """--gata4-exponent on --form c2-s2 (row R22's third T3 line): None unless the form is c2-s2 and the option is
    given (nothing is read; every other command, the ruling alone on c2-s2 included, goes on to check_gata4 as
    before); else the line's record, or Refused: a ruling other than R22, a value other than 2.0, the option without
    the ruling, or a ledger whose row R22 does not pass build_outside_form.captain_order_names for C2G_OPTION (line
    2's entry there, `--gata4-exponent 2.0` on `--form c1-s2`, does not open this line)."""
    if form != C2G_FORM or gata4_exponent is None:
        return None
    if captain_ruling is not None and captain_ruling != G2_ROW:
        raise Refused(f"--captain-ruling {captain_ruling!r} is refused: only row {G2_ROW} orders a file of this "
                      "builder")
    try:
        g = float(gata4_exponent)
    except (TypeError, ValueError):
        raise Refused(f"--gata4-exponent {gata4_exponent!r} is refused: only 2.0") from None
    if isinstance(gata4_exponent, bool) or g not in C2G_WEIGHTS:
        raise Refused(f"--gata4-exponent {gata4_exponent!r} is refused: only 2.0")
    if captain_ruling is None:
        raise Refused(f"--gata4-exponent 2.0 on --form {C2G_FORM} is taken only with --captain-ruling {G2_ROW}")
    p = TRIES_FILE if tries is None else Path(tries)
    if not _outside_form().captain_order_names(G2_ROW, C2G_OPTION, p):
        raise Refused(f"--form {C2G_FORM} --gata4-exponent 2.0 --captain-ruling {G2_ROW} is refused: "
                      f"experiments/harness/TRIES.md has no single row {G2_ROW} whose first cell is the captain's "
                      f"order (\"ruling, no run\") naming {C2G_OPTION} and `--captain-ruling {G2_ROW}`")
    return {"row": G2_ROW, "option": C2G_OPTION, "gata4_weight": g, "ledger": "experiments/harness/TRIES.md",
            "ledger_sha256": bs.sha256(p)}


def c2_gata4_weighted(phi_m, phi_g, g: float) -> tuple[np.ndarray, dict]:
    """C2's size-held product with the gata4_pooled weight g, from the two tables as loaded, gene by gene:
    phi_M ** c * phi_G ** (g c), c = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), population SD over the genes; nothing
    is rescaled or clipped after the product. g = 0 gives phi_M and g = 1 combine('c2', 1.0, ...)'s array, value for
    value (the same operations in the same order); at every g the SD of ln factor is phi_M's. Returns (phi, record);
    the record names g (the option's value on c2-s2: a weight, not a power) and the powers applied, c and g c."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    phi_g = np.asarray(phi_g, dtype=np.float64)
    if phi_m.shape != phi_g.shape or phi_m.ndim != 1:
        raise Refused("the two factor tables do not cover the same genes")
    for nm, p in ((TABLE_M, phi_m), (TABLE_G, phi_g)):
        if not np.all(np.isfinite(p) & (p > 0)):
            raise Refused(f"{nm}: factor values must be finite and positive")
    g = float(g)
    if not np.isfinite(g):
        raise Refused(f"the gata4_pooled weight {g!r} is not finite")
    lm, lg = np.log(phi_m), np.log(phi_g)
    sd_base, sd_sum = _sd(lm), _sd(lm + g * lg)
    if not sd_sum > 0:
        raise Refused("SD(ln phi_M + g ln phi_G) is 0: c is undefined")
    c = sd_base / sd_sum
    phi = (phi_m ** c) * (phi_g ** (g * c))
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused("the combined factor is not finite and positive")
    sd_m, sd_g = _sd(lm), _sd(lg)
    rec = {"kind": "c2", "x": 1.0, "g": g, "gata4_weight": g * c, "mab21l2_exponent": float(c), "c": float(c),
           "c_rule": "SD(ln phi_M) / SD(ln phi_M + g ln phi_G), population SD over the genes; on c2-s2 the option's "
                     "value is the weight g, not a power: the powers applied are c (transfer_mab21l2_all) and g c "
                     "(gata4_pooled)",
           "sd_ln_phi_m": sd_m, "sd_ln_phi_g": sd_g, "sd_ln_base": sd_base, "sd_ln_base_plus_wg": sd_sum,
           "corr_ln_m_ln_g": float(np.corrcoef(lm, lg)[0, 1]) if sd_m > 0 and sd_g > 0 else None,
           "base_factor": _shift(phi_m), "factor": _shift(phi),
           "den_weighted_mean": "recorded in combined_build.c2_gata4_line.guard",
           "min": float(phi.min()), "max": float(phi.max()), "n_up": int((phi > 1).sum()),
           "n_down": int((phi < 1).sum())}
    return phi, rec


def c2g_guard(phi_m, phi_g, g: float, den) -> dict:
    """The line's guard, a pure function of the two tables and den: the den-weighted totals (den . phi / sum den) of
    phi_M (B's factor), of C2's factor (combine('c2', 1.0), as C2's build applies it) and of the file's
    (c2_gata4_weighted at g); passes when the file's total lies within C2G_MAX_DEN_SHIFT of C2's. C2's shift and the
    file's shift from phi_M's total are recorded and gate nothing. Also the sha256 of the factor it describes."""
    phi, rec = c2_gata4_weighted(phi_m, phi_g, g)
    c2, _ = combine("c2", 1.0, phi_m, phi_g)
    phi_m = np.asarray(phi_m, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    if den.shape != phi_m.shape or not np.all(np.isfinite(den) & (den >= 0)) or not den.sum() > 0:
        raise Refused("den does not cover the tables' genes with finite, non-negative weights and a positive sum")
    s = float(den.sum())
    t_m, t_c2, t_f = (float(den @ f / s) for f in (phi_m, c2, phi))
    shift = t_f / t_c2 - 1.0
    return {"gata4_weight_g": float(g), "c": rec["c"],
            "powers_applied": {TABLE_M: rec["mab21l2_exponent"], TABLE_G: rec["gata4_weight"]},
            "den_total_phi_m": t_m, "den_total_c2": t_c2, "den_total_factor": t_f, "den_shift_vs_c2": shift,
            "max_den_shift": C2G_MAX_DEN_SHIFT, "passed": bool(abs(shift) <= C2G_MAX_DEN_SHIFT),
            "den_shift_c2_vs_phi_m": t_c2 / t_m - 1.0, "den_shift_factor_vs_phi_m": t_f / t_m - 1.0,
            "vs_phi_m": "recorded; gates nothing", "rule": C2G_GUARD_RULE, "factor_sha256": factor_sha256(phi)}


def run_c2g_guard(dec: dict, wt: Path, run_dir: Path | None, g: float) -> dict:
    """The guard from released WT E8.75 and the two pinned tables only, before any build: den and genes by
    build_outside_form.wt_den, both tables through the frozen loader in that gene order (checked against
    TABLE_SHA256). Returns the guard's record whatever it says; the caller refuses on it."""
    genes, den = _outside_form().wt_den(Path(wt))
    phi_m, rec_m = bs.load_factor(run_dir, TABLE_M, genes, dec)
    phi_g, rec_g = bs.load_factor(run_dir, TABLE_G, genes, dec)
    for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
        if rec["sha256"] != TABLE_SHA256[nm]:
            raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                          f"({TABLE_SHA256[nm][:16]}): not built")
    den = np.asarray(den, dtype=np.float64)
    out = c2g_guard(phi_m, phi_g, g, den)
    out |= {"wt": str(wt), "wt_sha256": bs.sha256(Path(wt)), "den_rule": G2_DEN_RULE, "den": den.tolist(),
            "den_sha256": factor_sha256(den), "n_genes": len(genes),
            "genes_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest(),
            "tables_sha256": {TABLE_M: rec_m["sha256"], TABLE_G: rec_g["sha256"]}, "genes": list(genes)}
    return out


@contextmanager
def c2g_loader(g: float, out: dict, guard: dict):
    """bs.load_factor, for this process only, on row R22's third line: the frozen loader for phi_M and phi_G (checked
    against TABLE_SHA256 too), then c2_gata4_weighted at g, the hook C2's factor goes through. FactorMismatch (exit 5)
    unless the build's genes are the guard's and the factor's sha256 is the one the guard describes."""
    original = bs.load_factor

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        if list(genes) != guard["genes"]:
            raise FactorMismatch("the build's genes are not those the guard was computed over: not built")
        phi_m, rec_m = original(run_dir, TABLE_M, genes, dec)
        phi_g, rec_g = original(run_dir, TABLE_G, genes, dec)
        for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
            if rec["sha256"] != TABLE_SHA256[nm]:
                raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                              f"({TABLE_SHA256[nm][:16]}): not built")
        phi, crec = c2_gata4_weighted(phi_m, phi_g, g)
        sha = factor_sha256(phi)
        if sha != guard["factor_sha256"]:
            raise FactorMismatch(f"the factor the build would apply ({sha[:16]}) is not the one the guard describes "
                                 f"({guard['factor_sha256'][:16]}): not built")
        out.update(tables={TABLE_M: rec_m, TABLE_G: rec_g}, combination=crec, factor_sha256=sha,
                   n_loads=out.get("n_loads", 0) + 1)
        return phi, rec_m

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


# ------------------------------------------------- the options on C2's command: F1, F2 and F3 (off by default)


def _is_pos_inf(v) -> bool:
    """Whether a --gata4-exponent value is +inf, read as the option reads every value (float(); never a bool)."""
    if v is None or isinstance(v, bool):
        return False
    try:
        return float(v) == math.inf
    except (TypeError, ValueError, OverflowError):
        return False


def is_gonly_command(form, gata4_exponent=None, size_scale=None) -> bool:
    """Whether a command is this line's: any --size-scale, or --gata4-exponent inf on --form c2-s2. Every other command,
    --gata4-exponent inf on another form included, never reaches the line's code and goes on as before."""
    return size_scale is not None or (form == GONLY_FORM and _is_pos_inf(gata4_exponent))


def ledger_first_cell(row: str, path: Path) -> str | None:
    """The first cell of the ledger's single row `row`, read as build_outside_form.captain_order_names reads it (the
    lines that start "| <row> "); None when the ledger holds no such row or more than one."""
    p = Path(path)
    if not p.is_file():
        return None
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {row} ")]
    return lines[0][2:].split("|", 1)[0] if len(lines) == 1 else None


def names_whole(first: str | None, option: str) -> bool:
    """Whether a first cell holds `option` (given without its backticks) as one whole backtick-delimited string. A
    string that only contains it (F3's contains F1's) or that it contains does not count."""
    return first is not None and option in re.findall(r"`([^`]*)`", first)


def check_den_lift(name: str, den_ruling, tries: Path | None = None) -> dict:
    """--den-ruling ROW for the file `name` (F1, F2 or F3): the lifting's record, or Refused unless ROW is a ruling
    row's id (R and digits) and the ledger holds exactly one row ROW whose first cell says "ruling, no run" and "the 5%
    check gates nothing for it" and holds the file's option string as one whole backtick-delimited string. Row R30
    lifted line 2's check for its one file by name; here the row is named on the command line and must name the file."""
    if not isinstance(den_ruling, str) or not GONLY_LIFT_ROW.fullmatch(den_ruling):
        raise Refused(f"--den-ruling {den_ruling!r} is refused: it takes a ruling row's id (R and digits)")
    p = TRIES_FILE if tries is None else Path(tries)
    first = ledger_first_cell(den_ruling, p)
    opt = GONLY_OPTIONS[name]
    if first is None or "ruling, no run" not in first or GONLY_LIFT_WORDS not in first or not names_whole(first, opt):
        raise Refused(f"--den-ruling {den_ruling} is refused: experiments/harness/TRIES.md has no single row "
                      f"{den_ruling} whose first cell is a ruling (\"ruling, no run\") naming `{opt}` as one whole "
                      f"string and saying \"{GONLY_LIFT_WORDS}\"")
    return {"row": den_ruling, "option": f"`{opt}`", "words": GONLY_LIFT_WORDS,
            "ledger": "experiments/harness/TRIES.md", "ledger_sha256": bs.sha256(p)}


def check_c2_gonly(form, gata4_exponent=None, size_scale=None, captain_ruling=None, tries: Path | None = None,
                   composition_from=None, den_ruling=None) -> dict | None:
    """--gata4-exponent inf and --size-scale on --form c2-s2 (and --den-ruling): None unless is_gonly_command or
    --den-ruling is given (nothing is read); else the line's record, or Refused: --den-ruling on any other command; a
    ruling other than R22; --size-scale on another form or with a value other than 0.67; --size-scale beside a
    --gata4-exponent other than inf; either option without --captain-ruling R22 or beside --composition-from; a ledger
    whose row R22 does not pass build_outside_form.captain_order_names for the file's option string or does not hold
    that string as one whole backtick-delimited string; a --den-ruling that check_den_lift refuses."""
    line = is_gonly_command(form, gata4_exponent, size_scale)
    if not line and den_ruling is None:
        return None
    if not line:
        raise Refused("--den-ruling is taken only with the options of F1, F2 or F3 (--gata4-exponent inf and/or "
                      f"--size-scale 0.67 on --form {GONLY_FORM})")
    if captain_ruling is not None and captain_ruling != G2_ROW:
        raise Refused(f"--captain-ruling {captain_ruling!r} is refused: only row {G2_ROW} orders a file of this "
                      "builder")
    if form != GONLY_FORM:
        raise Refused(f"--size-scale is taken on --form {GONLY_FORM} only, not {form!r}")
    s = 1.0
    if size_scale is not None:
        try:
            s = float(size_scale)
        except (TypeError, ValueError):
            raise Refused(f"--size-scale {size_scale!r} is refused: only 0.67") from None
        if isinstance(size_scale, bool) or s not in GONLY_SCALES:
            raise Refused(f"--size-scale {size_scale!r} is refused: only 0.67")
    if gata4_exponent is not None and not _is_pos_inf(gata4_exponent):
        raise Refused(f"--size-scale 0.67 is taken alone (C2's g = 1) or with --gata4-exponent inf, not with "
                      f"--gata4-exponent {gata4_exponent!r}")
    g = math.inf if gata4_exponent is not None else 1.0
    name = next(k for k, v in GONLY_POINTS.items() if v == (g, s))
    flags = {"F1": "--gata4-exponent inf", "F2": "--size-scale 0.67", "F3": "--gata4-exponent inf --size-scale 0.67"}
    if composition_from is not None:
        raise Refused(f"{flags[name]} takes no --composition-from (the file is C2's command with its own options)")
    if captain_ruling is None:
        raise Refused(f"{flags[name]} on --form {GONLY_FORM} is taken only with --captain-ruling {G2_ROW}")
    p = TRIES_FILE if tries is None else Path(tries)
    opt = GONLY_OPTIONS[name]
    if not (_outside_form().captain_order_names(G2_ROW, f"`{opt}`", p) and names_whole(ledger_first_cell(G2_ROW, p),
                                                                                         opt)):
        raise Refused(f"{opt} --captain-ruling {G2_ROW} is refused: experiments/harness/TRIES.md has no single row "
                      f"{G2_ROW} whose first cell is the captain's order (\"ruling, no run\") naming `{opt}` as one "
                      f"whole string and `--captain-ruling {G2_ROW}`")
    lift = None if den_ruling is None else check_den_lift(name, den_ruling, p)
    return {"row": G2_ROW, "file": GONLY_FILES[name], "name": name, "option": f"`{opt}`",
            "g": "inf" if g == math.inf else g, "size_scale": s, "den_ruling": lift,
            "ledger": "experiments/harness/TRIES.md", "ledger_sha256": bs.sha256(p)}


def gonly_point(rec: dict) -> tuple[float, float]:
    """(g, s) of a record of check_c2_gonly (g = math.inf for "inf")."""
    return GONLY_POINTS[rec["name"]]


def c2_gonly_factor(phi_m, phi_g, g: float, s: float = 1.0) -> tuple[np.ndarray, dict]:
    """The family ln phi(g, s) = s c_g (ln phi_M + g ln phi_G) on C2's base at this line's points, from the two tables
    as loaded, gene by gene: g = inf is phi_G alone at C2's held size, phi_G ** c_inf with c_inf = SD(ln phi_M) /
    SD(ln phi_G) (phi_M dropped), and g = 1 is C2's factor (combine('c2', 1.0), value for value); then, unless s is
    1, ln phi multiplied by s after that size-held combination (** s). Population SD over the genes; nothing rescaled
    or clipped after. Returns (phi, record)."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    phi_g = np.asarray(phi_g, dtype=np.float64)
    if phi_m.shape != phi_g.shape or phi_m.ndim != 1:
        raise Refused("the two factor tables do not cover the same genes")
    for nm, p in ((TABLE_M, phi_m), (TABLE_G, phi_g)):
        if not np.all(np.isfinite(p) & (p > 0)):
            raise Refused(f"{nm}: factor values must be finite and positive")
    g, s = float(g), float(s)
    if not (np.isfinite(s) and s > 0):
        raise Refused(f"the size scale {s!r} is not finite and positive")
    lm, lg = np.log(phi_m), np.log(phi_g)
    sd_m, sd_g = _sd(lm), _sd(lg)
    if g == math.inf:
        if not sd_g > 0:
            raise Refused("SD(ln phi_G) is 0: c_inf is undefined")
        c = sd_m / sd_g
        held = phi_g ** c
        powers = (0.0, c)
        rule = "g = inf: c_inf = SD(ln phi_M) / SD(ln phi_G); ln phi = c_inf ln phi_G (phi_M dropped)"
    elif g == 1.0:
        held, crec = combine("c2", 1.0, phi_m, phi_g)
        c = crec["c"]
        powers = (c, c)
        rule = "g = 1: C2's factor, c = SD(ln phi_M) / SD(ln phi_M + ln phi_G)"
    else:
        raise Refused(f"the gata4_pooled weight {g!r} is not this line's (inf, or 1 under a size scale)")
    phi = held if s == 1.0 else held ** s
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused("the combined factor is not finite and positive")
    rec = {"kind": "c2", "x": 1.0, "g": "inf" if g == math.inf else g, "size_scale": s,
           "gata4_weight": s * powers[1], "mab21l2_exponent": s * powers[0], "c": float(c),
           "c_rule": f"{rule}; population SD over the genes; then ln phi multiplied by the size scale s (** s, none "
                     "at s = 1)",
           "sd_ln_phi_m": sd_m, "sd_ln_phi_g": sd_g, "sd_ln_size_held": _sd(np.log(held)),
           "base_factor": _shift(phi_m), "size_held_factor": _shift(held), "factor": _shift(phi),
           "den_weighted_mean": "recorded in combined_build.c2_gonly_line.guards",
           "min": float(phi.min()), "max": float(phi.max()), "n_up": int((phi > 1).sum()),
           "n_down": int((phi < 1).sum())}
    return phi, rec


def gonly_guards(phi_m, phi_g, g: float, s: float, den, lift: dict | None = None) -> dict:
    """The line's two checks, a pure function of the two tables and den: each is a property of the file's own factor
    (c2_gonly_factor) against C2's (combine('c2', 1.0), as C2's build applies it); no cell is drawn. Range: the file's
    max at most GONLY_RANGE_FACTOR times C2's max and its min at least C2's min over GONLY_RANGE_FACTOR. Den, line 3's
    check (c2g_guard's): the file's den-weighted total (den . phi / sum den) within C2G_MAX_DEN_SHIFT of C2's; it gates
    unless `lift` (check_den_lift's record of the row that names the file) is given, when it is recorded and gates
    nothing. The shifts from phi_M's total are recorded and gate nothing. Also the sha256 of the factor described."""
    phi, rec = c2_gonly_factor(phi_m, phi_g, g, s)
    c2, _ = combine("c2", 1.0, phi_m, phi_g)
    phi_m = np.asarray(phi_m, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    if den.shape != phi_m.shape or not np.all(np.isfinite(den) & (den >= 0)) or not den.sum() > 0:
        raise Refused("den does not cover the tables' genes with finite, non-negative weights and a positive sum")
    lo, hi = float(c2.min()) / GONLY_RANGE_FACTOR, float(c2.max()) * GONLY_RANGE_FACTOR
    fmin, fmax = float(phi.min()), float(phi.max())
    range_ok = bool(fmin >= lo and fmax <= hi)
    tot = float(den.sum())
    t_m, t_c2, t_f = (float(den @ f / tot) for f in (phi_m, c2, phi))
    shift = t_f / t_c2 - 1.0
    within = bool(abs(shift) <= C2G_MAX_DEN_SHIFT)
    den_ok = bool(within or lift is not None)
    status = "gates" if lift is None else f"lifted by row {lift['row']}: {GONLY_LIFT_WORDS}"
    return {"g": rec["g"], "size_scale": rec["size_scale"], "c": rec["c"],
            "powers_applied": {TABLE_M: rec["mab21l2_exponent"], TABLE_G: rec["gata4_weight"]},
            "range": {"factor_min": fmin, "factor_max": fmax, "c2_min": float(c2.min()), "c2_max": float(c2.max()),
                      "bound_min": lo, "bound_max": hi, "factor": GONLY_RANGE_FACTOR, "passed": range_ok,
                      "rule": GONLY_RANGE_RULE},
            "den_check": {"den_total_phi_m": t_m, "den_total_c2": t_c2, "den_total_factor": t_f,
                          "den_shift_vs_c2": shift, "max_den_shift": C2G_MAX_DEN_SHIFT, "within": within,
                          "limit": {"max_abs_den_shift": C2G_MAX_DEN_SHIFT, "status": status}, "lifted_by": lift,
                          "passed": den_ok, "den_shift_c2_vs_phi_m": t_c2 / t_m - 1.0,
                          "den_shift_factor_vs_phi_m": t_f / t_m - 1.0, "vs_phi_m": "recorded; gates nothing",
                          "rule": GONLY_DEN_RULE},
            "passed": bool(range_ok and den_ok), "factor_sha256": factor_sha256(phi)}


def gonly_failure(guards: dict) -> str:
    """What a failed check says (its own figures; no cell is drawn)."""
    says = []
    r = guards["range"]
    if not r["passed"]:
        says.append(f"the range guard: the factor spans [{r['factor_min']:.4g}, {r['factor_max']:.4g}], not within "
                    f"[{r['bound_min']:.4g}, {r['bound_max']:.4g}] (half C2's min, twice C2's max)")
    d = guards["den_check"]
    if not d["passed"]:
        says.append(f"the 5% den check: the file's den-weighted total is {d['den_shift_vs_c2']:+.2%} against C2's "
                    f"factor's, beyond {C2G_MAX_DEN_SHIFT:.0%}, and no --den-ruling row lifts it for this file")
    return "; ".join(says) + ": not built"


def run_gonly_guards(dec: dict, wt: Path, run_dir: Path | None, g: float, s: float, lift: dict | None = None) -> dict:
    """The checks from released WT E8.75's den and the two pinned tables only, before any build: den and genes by
    build_outside_form.wt_den, both tables through the frozen loader in that gene order (checked against TABLE_SHA256).
    Returns the record whatever it says; the caller refuses on it."""
    genes, den = _outside_form().wt_den(Path(wt))
    phi_m, rec_m = bs.load_factor(run_dir, TABLE_M, genes, dec)
    phi_g, rec_g = bs.load_factor(run_dir, TABLE_G, genes, dec)
    for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
        if rec["sha256"] != TABLE_SHA256[nm]:
            raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                          f"({TABLE_SHA256[nm][:16]}): not built")
    den = np.asarray(den, dtype=np.float64)
    out = gonly_guards(phi_m, phi_g, g, s, den, lift)
    out |= {"wt": str(wt), "wt_sha256": bs.sha256(Path(wt)), "den_rule": G2_DEN_RULE, "den": den.tolist(),
            "den_sha256": factor_sha256(den), "n_genes": len(genes),
            "genes_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest(),
            "tables_sha256": {TABLE_M: rec_m["sha256"], TABLE_G: rec_g["sha256"]}, "genes": list(genes)}
    return out


@contextmanager
def gonly_loader(g: float, s: float, out: dict, guards: dict):
    """bs.load_factor, for this process only, on this line: the frozen loader for phi_M and phi_G (checked against
    TABLE_SHA256 too), then c2_gonly_factor at (g, s), the hook C2's factor goes through. FactorMismatch (exit 5) unless
    the build's genes are the guards' and the factor's sha256 is the one the guards describe."""
    original = bs.load_factor

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        if list(genes) != guards["genes"]:
            raise FactorMismatch("the build's genes are not those the guards were computed over: not built")
        phi_m, rec_m = original(run_dir, TABLE_M, genes, dec)
        phi_g, rec_g = original(run_dir, TABLE_G, genes, dec)
        for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
            if rec["sha256"] != TABLE_SHA256[nm]:
                raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                              f"({TABLE_SHA256[nm][:16]}): not built")
        phi, crec = c2_gonly_factor(phi_m, phi_g, g, s)
        sha = factor_sha256(phi)
        if sha != guards["factor_sha256"]:
            raise FactorMismatch(f"the factor the build would apply ({sha[:16]}) is not the one the guards describe "
                                 f"({guards['factor_sha256'][:16]}): not built")
        out.update(tables={TABLE_M: rec_m, TABLE_G: rec_g}, combination=crec, factor_sha256=sha,
                   n_loads=out.get("n_loads", 0) + 1)
        return phi, rec_m

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


# --------------------------------------------- the depth step on C2's command: dd, ddh and dg (off by default)

# The T3 ERROR design of the design round of 2026-11-25 (proposal P1, dD), as its cold challenge fixes it: C2's command
# with one change, the frozen allowed branch's per-label depth step (only submissions.allowed.s1 carries `depth: true`,
# switched on by rule D4), added to the barred s2's spec where the frozen builder resolves it (as r36_composition
# replaces a composition block), its ratios read from the pinned decision file, dd's as they are and ddh's and dg's in
# two new forms, each file only with --captain-ruling R53 while row R53's first cell is the captain's order naming the
# file's whole option string (row R53: the captain's order of 2026-11-24, his named exception on T3's depth for the
# three files, bound to row K30's read). Off by default; its own code below, beside the other lines, which it leaves
# unchanged.
DEPTH_FORM, DEPTH_BASE = "c2-s2", "s2"
DEPTH_ROW = "R53"  # the one ruling row that orders the line's files
DEPTH_SOURCES = {"allowed-s1": ("allowed", "s1")}  # --depth-from -> the decision file's entry that carries the step
DEPTH_EXPONENTS = (0.5,)  # --depth-exponent: the single value 0.5 (ddh)
DEPTH_FORMS = ("global",)  # --depth-form: the single value global (dg)
# file -> its option string, which row R53's first cell holds as one whole backtick-delimited string (dd's string sits
# inside ddh's and dg's, so containment would let either entry open dd; the whole-string match does not); the row names
# the form, `--form c2-s2`, as a string of its own
DEPTH_OPTIONS = {"dd": "--depth-from allowed-s1",
                 "ddh": "--depth-from allowed-s1 --depth-exponent 0.5",
                 "dg": "--depth-from allowed-s1 --depth-form global"}
DEPTH_FORM_OPTION = "--form c2-s2"
DEPTH_FILES = {"dd": "T3_c2_s2_dd.h5ad", "ddh": "T3_c2_s2_ddh.h5ad", "dg": "T3_c2_s2_dg.h5ad"}
# experiments/t3_decomp/results/gata4/gata4_decisions.json, the frozen run of 2026-10-28 (the file row R36 pins too)
DEPTH_DECISION_SHA256 = "9440ae67625608112666d5235e6ceaa5dd9f47a4b0bbd30e16943f32f35c5a8c"
DEPTH_PARAMS_SOURCE = "decision file submissions.parameters"  # bs.run_parameters' name for the decision file's own
EXIT_DEPTH_CHECK = 7
DEPTH_LINE = ("T3 depth line (the T3 ERROR design of 2026-11-25, proposal P1, dD, as its cold challenge fixes it): C2 "
              "on the barred s2 (--form c2-s2) with one change, the frozen allowed branch's per-label depth step "
              "(decision file submissions.allowed.s1, depth: true; its ratios "
              "submissions.parameters.depth_ratio_pooled) added to the barred s2's spec and applied as construct's "
              "ratio after C2's factor on the same generator; dd T3_c2_s2_dd.h5ad (the frozen ratios as they are), "
              "ddh T3_c2_s2_ddh.h5ad (the ratios ** 0.5, a new form) and dg T3_c2_s2_dg.h5ad (one ratio for every "
              "label, a new form); row R53")
DEPTH_SOURCE = ("the depth step of decision file submissions.allowed.s1 (depth: true; ratios "
                "submissions.parameters.depth_ratio_pooled), added to the barred s2's spec by --depth-from allowed-s1")
DEPTH_GLOBAL_RULE = ("one ratio for every label, exp(sum_t n_t ln r_t / sum_t n_t): r_t label t's ratio as "
                     "build_submission.depth_ratios gives it (1.0 for a label the frozen rule left unchanged), n_t WT "
                     "E8.75's cells in label t as the build's reference holds them (Vocab.count, duplicates dropped)")
DEPTH_EMPTY_RULE = ("as the frozen step: t3_decompose.rescale_depth, unchanged, thins through t3_decompose.thin, in "
                    "which a cell whose binomial draw would be empty keeps its counts (it stays unthinned), so no cell "
                    "reaches lognorm empty through the step (a cell empty before it, none expected, stays empty and "
                    "lognorm's where=n > 0 guard writes it as zeros); the cells that stay unthinned are counted on the "
                    "very draw (a recording stand-in passes every call to the real generator unchanged)")
DEPTH_WHERE = ("bs.resolve_spec (for this process only) adds depth: true to the barred s2's spec, so the frozen "
               "inputs_for, run_parameters and prediction_inputs read the decision file's parameters and call "
               "bs.depth_ratios (for this process only: the frozen reading, then the file's form); t3_decompose.build "
               "applies the ratios after the factor, on the same generator (rescale_depth: binomial thinning below 1, "
               "Poisson top-up from the kNN-smoothed profile above 1), so C2's draws (rows, factor) are untouched; "
               "inside the build, before anything is written, C2's cells are drawn again on a copy of the generator in "
               "the same state with no depth step, and the counts the depth step receives must be C2's final counts "
               "row by row (same rows, labels and sources) and each cell's ratio its label's")
DEPTH_LIST_ENTRY = ("experiments/harness/TRIES.md row R53 (the captain's order of 2026-11-24: his named exception on "
                    "T3's depth for T3_c2_s2_dd, T3_c2_s2_ddh and T3_c2_s2_dg, bound to row K30's read; decision log "
                    "2026-11-24 (rrr)); the T3 ERROR design of the design round of 2026-11-25, section 3 (proposal P1, "
                    "dD), and its cold challenge")


def is_depth_command(depth_from=None, depth_exponent=None, depth_form=None) -> bool:
    """Whether a command is the depth line's: any of --depth-from, --depth-exponent and --depth-form. Every other
    command never reaches the line's code and goes on as before."""
    return depth_from is not None or depth_exponent is not None or depth_form is not None


def check_depth_line(form, depth_from=None, depth_exponent=None, depth_form=None, captain_ruling=None,
                     tries: Path | None = None, composition_from=None, gata4_exponent=None, size_scale=None,
                     den_ruling=None) -> dict | None:
    """--depth-from allowed-s1 on --form c2-s2, with --depth-exponent 0.5 or --depth-form global as the variants:
    None unless is_depth_command (nothing is read); else the line's record (the file, its option string, the ruling
    row), or Refused: --depth-exponent or --depth-form without --depth-from; any value of --depth-from but allowed-s1,
    of --depth-exponent but 0.5, of --depth-form but global; any form but c2-s2; --depth-exponent and --depth-form
    together; --composition-from, --gata4-exponent, --size-scale or --den-ruling beside it; no --captain-ruling, or one
    other than R53; a ledger whose row R53 does not pass build_outside_form.captain_order_names for the file's option
    string (the captain's order, "ruling, no run", `--captain-ruling R53`) or does not hold that string and
    `--form c2-s2` each as one whole backtick-delimited string."""
    if not is_depth_command(depth_from, depth_exponent, depth_form):
        return None
    if depth_from is None:
        raise Refused("--depth-exponent and --depth-form are taken only with --depth-from allowed-s1")
    if not isinstance(depth_from, str) or depth_from not in DEPTH_SOURCES:
        raise Refused(f"--depth-from {depth_from!r} is refused: only allowed-s1")
    if form != DEPTH_FORM:
        raise Refused(f"--depth-from allowed-s1 is taken on --form {DEPTH_FORM} only, not {form!r}")
    e = None
    if depth_exponent is not None:
        try:
            e = float(depth_exponent)
        except (TypeError, ValueError, OverflowError):
            raise Refused(f"--depth-exponent {depth_exponent!r} is refused: only 0.5") from None
        if isinstance(depth_exponent, bool) or e not in DEPTH_EXPONENTS:
            raise Refused(f"--depth-exponent {depth_exponent!r} is refused: only 0.5")
    if depth_form is not None and (not isinstance(depth_form, str) or depth_form not in DEPTH_FORMS):
        raise Refused(f"--depth-form {depth_form!r} is refused: only global")
    if e is not None and depth_form is not None:
        raise Refused("--depth-exponent and --depth-form are not taken together (one variant per file)")
    for flag, value in (("--composition-from", composition_from), ("--gata4-exponent", gata4_exponent),
                        ("--size-scale", size_scale), ("--den-ruling", den_ruling)):
        if value is not None:
            raise Refused(f"--depth-from allowed-s1 takes no {flag} (the file is C2's command with one change)")
    if captain_ruling is None:
        raise Refused(f"--depth-from allowed-s1 is taken only with --captain-ruling {DEPTH_ROW}")
    if captain_ruling != DEPTH_ROW:
        raise Refused(f"--captain-ruling {captain_ruling!r} is refused on the depth line: only row {DEPTH_ROW} orders "
                      "its files")
    arm = "dg" if depth_form is not None else "ddh" if e is not None else "dd"
    opt = DEPTH_OPTIONS[arm]
    p = TRIES_FILE if tries is None else Path(tries)
    first = ledger_first_cell(DEPTH_ROW, p)
    if not (_outside_form().captain_order_names(DEPTH_ROW, f"`{opt}`", p) and names_whole(first, opt)
            and names_whole(first, DEPTH_FORM_OPTION)):
        raise Refused(f"--form {DEPTH_FORM} {opt} --captain-ruling {DEPTH_ROW} is refused: experiments/harness/"
                      f"TRIES.md has no single row {DEPTH_ROW} whose first cell is the captain's order (\"ruling, no "
                      f"run\") naming `{opt}` and `{DEPTH_FORM_OPTION}` each as one whole string and `--captain-ruling "
                      f"{DEPTH_ROW}`")
    return {"row": DEPTH_ROW, "arm": arm, "file": DEPTH_FILES[arm], "option": f"`{opt}`",
            "form_option": f"`{DEPTH_FORM_OPTION}`", "depth_from": depth_from, "exponent": e,
            "form": "global" if depth_form is not None else "per label", "ledger": "experiments/harness/TRIES.md",
            "ledger_sha256": bs.sha256(p)}


def depth_parameters(dec: dict, decision: Path, depth_from: str = "allowed-s1") -> tuple[dict, dict]:
    """The frozen depth ratios, and nothing else of the decision file's parameters: submissions.parameters.
    depth_ratio_pooled as bs.run_parameters reads it from the decision file itself (no run folder is read). Refused
    unless the decision file is the pinned one (DEPTH_DECISION_SHA256), the entry --depth-from names
    (submissions.allowed.s1) carries the depth step (depth: true) and every ratio is a finite positive number. The
    entry's factor (gata4), composition, stage and coverage are never read. Returns ({"depth_ratio_pooled": ratios},
    its record)."""
    if not isinstance(depth_from, str) or depth_from not in DEPTH_SOURCES:
        raise Refused(f"--depth-from {depth_from!r} is refused: only allowed-s1")
    branch, slot = DEPTH_SOURCES[depth_from]
    got = bs.sha256(Path(decision))
    if got != DEPTH_DECISION_SHA256:
        raise Refused(f"the decision file's sha256 {got[:16]} is not the pinned one ({DEPTH_DECISION_SHA256[:16]}): "
                      "not built")
    entry = ((dec.get("submissions") or {}).get(branch) or {}).get(slot) if isinstance(dec.get("submissions"),
                                                                                         dict) else None
    if not isinstance(entry, dict) or entry.get("depth") is not True:
        raise Refused(f"submissions.{branch}.{slot} of the decision file does not carry the depth step (depth: true): "
                      "not built")
    params, src = bs.run_parameters(dec, None)
    pooled = params.get("depth_ratio_pooled") if src == DEPTH_PARAMS_SOURCE else None
    if not isinstance(pooled, dict) or not pooled:
        raise Refused("the decision file has no submissions.parameters.depth_ratio_pooled: not built")
    for lab, v in pooled.items():
        if isinstance(v, bool) or not isinstance(v, int | float) or not math.isfinite(v) or not v > 0:
            raise Refused(f"the depth ratio {v!r} of {lab!r} is not a finite positive number: not built")
    pooled = {str(k): float(v) for k, v in pooled.items()}
    text = "\n".join(f"{k}\t{v!r}" for k, v in pooled.items())
    return {"depth_ratio_pooled": pooled}, {
        "source": f"decision file submissions.parameters.depth_ratio_pooled (the step submissions.{branch}.{slot} "
                  "carries, depth: true)", "entry": f"submissions.{branch}.{slot}", "decision": str(decision),
        "decision_sha256": got, "pinned_sha256": DEPTH_DECISION_SHA256, "n_labels": len(pooled), "pooled": pooled,
        "pooled_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "pooled_sha256_rule": "label, tab, repr(ratio), one line per label in the decision file's order"}


def depth_ratios_applied(V, params: dict, exponent=None, form=None, reader=None) -> tuple[np.ndarray, dict]:
    """The ratios the line passes as construct's `ratio`, one per label of the reference vocabulary V, from the frozen
    ratios as the frozen builder reads them (`reader`, default bs.depth_ratios(V, params): 1.0 for a label without a
    pooled ratio): as they are (dd), to the power 0.5 (ddh), or one ratio for every label (dg), DEPTH_GLOBAL_RULE
    (all of WT E8.75's labels, 1.0 where none; weights V.count, WT E8.75's cells per label in the reference the build
    draws from). Returns (ratios, record); nothing in the record is a count of the file's cells."""
    r = np.asarray((bs.depth_ratios if reader is None else reader)(V, params), dtype=np.float64)
    if r.shape != (V.L,) or not np.all(np.isfinite(r) & (r > 0)):
        raise Refused("the frozen depth ratios are not finite and positive on the reference's labels: not built")
    with_ratio = np.zeros(V.L, dtype=bool)
    with_ratio[V.index(list(params["depth_ratio_pooled"]))] = True
    rec = {"labels": list(V.labels), "labels_sha256": hashlib.sha256("\n".join(V.labels).encode()).hexdigest(),
           "n_labels": int(V.L), "n_labels_with_a_pooled_ratio": int(with_ratio.sum()),
           "wt_share_of_labels_left_at_1": float(np.asarray(V.p, dtype=np.float64)[~with_ratio].sum()),
           "frozen": [float(x) for x in r], "frozen_sha256": factor_sha256(r),
           "sha256_rule": "sha256 of the float64 bytes, one value per label in the reference vocabulary's order"}
    if form == "global":
        n = np.asarray(V.count, dtype=np.float64)
        if n.shape != r.shape or not np.all(np.isfinite(n) & (n >= 0)) or not n.sum() > 0:
            raise Refused("WT E8.75's cells per label do not give weights for the global ratio: not built")
        g = float(np.exp(n @ np.log(r) / n.sum()))
        applied = np.full(V.L, g)
        rec |= {"form": "global", "exponent": None, "global_ratio": g, "rule": DEPTH_GLOBAL_RULE}
    elif exponent is not None:
        e = float(exponent)
        if e not in DEPTH_EXPONENTS:
            raise Refused(f"the depth exponent {exponent!r} is not this line's (0.5)")
        applied = r ** e
        rec |= {"form": "per label", "exponent": e, "global_ratio": None, "rule": "r_t ** 0.5, label by label"}
    elif form is None:
        applied = r
        rec |= {"form": "per label", "exponent": 1.0, "global_ratio": None,
                "rule": "r_t, the frozen ratios as they are"}
    else:
        raise Refused(f"the depth form {form!r} is not this line's (global)")
    if not np.all(np.isfinite(applied) & (applied > 0)):
        raise Refused("the depth ratios applied are not finite and positive: not built")
    rec |= {"applied": [float(x) for x in applied], "applied_sha256": factor_sha256(applied)}
    return applied, rec


class _DrawRecorder:
    """The generator as t3_decompose.rescale_depth sees it on the depth line: every call goes to the real generator
    unchanged (the same method, arguments and order, so the stream is the frozen step's own); the row sums of each
    binomial draw are kept, from which the cells the frozen thinning would empty are counted."""

    def __init__(self, rng):
        self._rng = rng
        self.binomial_row_sums: list = []

    def binomial(self, *args, **kwargs):
        out = self._rng.binomial(*args, **kwargs)
        self.binomial_row_sums.append(np.asarray(out).sum(axis=-1))
        return out

    def __getattr__(self, name):
        return getattr(self._rng, name)


def _twin_generator(rng):
    """A copy of the generator in the same state: C2's draws are made again on it, the real one is left to the file."""
    return copy.deepcopy(rng)


@contextmanager
def _observe(name: str, seen: list, wrap):
    """t3_decompose.<name>, for this process only, replaced by wrap(original, seen)."""
    original = getattr(bs.t3d, name)
    setattr(bs.t3d, name, wrap(original, seen))
    try:
        yield
    finally:
        setattr(bs.t3d, name, original)


def _lognorm_spy(original, seen: list):
    def spy(C):
        seen.append(C)
        return original(C)

    return spy


def _rescale_spy(original, seen: list):
    def spy(C, ratio, P, rng):
        rec = _DrawRecorder(rng)
        res = original(C, ratio, P, rec)
        seen.append({"pre": C, "ratio": np.asarray(ratio, dtype=np.float64), "post": res,
                     "draws": rec.binomial_row_sums})
        return res

    return spy


def depth_checks(c2, c2_counts, built, call: dict, ratios) -> dict:
    """The line's checks of one build, a pure function of C2's cells (drawn again with no depth step: Built and the
    counts its lognorm received) and the file's (Built and what rescale_depth received and returned): whether the
    counts the depth step receives are C2's final counts row by row, whether the rows, labels and sources are C2's,
    and whether each cell's ratio is its label's; and the counts of cells the thinning would empty, which keep their
    counts (t3_decompose.thin: they stay unthinned), of cells empty before the step and after it."""
    pre, post = np.asarray(call["pre"]), np.asarray(call["post"])
    c2_counts = np.asarray(c2_counts)
    draws = call["draws"]
    out = {"pre_depth_counts_are_c2s": bool(pre.shape == c2_counts.shape and np.array_equal(pre, c2_counts)),
           "rows_are_c2s": bool(np.array_equal(built.row, c2.row)),
           "labels_are_c2s": bool(np.array_equal(built.lab, c2.lab)),
           "sources_are_c2s": bool(np.array_equal(built.source, c2.source)),
           "ratio_is_the_label_s": bool(call["ratio"].shape == built.lab.shape
                                        and np.array_equal(call["ratio"], np.asarray(ratios)[built.lab])),
           "one_thinning_draw": len(draws) == 1}
    kept = ((draws[0] == 0) & (pre.sum(1) > 0)) if len(draws) == 1 else np.zeros(pre.shape[0], dtype=bool)
    out["counts"] = {"stay_unthinned": int(kept.sum()), "empty_before": int((pre.sum(1) == 0).sum()),
                     "empty_after": int((post.sum(1) == 0).sum())}
    out["passed"] = all(out[k] for k in ("pre_depth_counts_are_c2s", "rows_are_c2s", "labels_are_c2s",
                                         "sources_are_c2s", "ratio_is_the_label_s", "one_thinning_draw"))
    return out


@contextmanager
def depth_step(params: dict, exponent, form, out: dict):
    """bs.resolve_spec, bs.depth_ratios and bs.construct, for this process only, on the depth line. resolve_spec: the
    frozen one, then the barred s2's spec with depth: true added (any other slot refused, and a spec that is not C2's:
    its factor alone and a composition, no depth, stage or coverage step), as r36_composition replaces a composition
    block; the frozen validation, bs.inputs_for (bs.run_parameters: the decision file's submissions.parameters) and
    bs.prediction_inputs then take the step as the allowed s1 takes it. depth_ratios: refused unless the ratios the
    build reads are the pinned decision file's (`params`); then the frozen reading in the file's form
    (depth_ratios_applied). construct: C2's cells first, on a copy of the generator in the same state and with no depth
    step (the counts lognorm receives captured), then the file's on the real generator with
    t3_decompose.rescale_depth observed (_DrawRecorder); DepthCheckFailed (exit 7) before anything is written unless
    depth_checks passes. `out` receives the spec replaced, the ratios' record, the checks and the call counts."""
    orig_resolve, orig_ratios, orig_construct = bs.resolve_spec, bs.depth_ratios, bs.construct

    def resolve(dec, rules, branch, submission, hashes):
        if (branch, submission) != (BRANCH, DEPTH_BASE):
            raise Refused(f"--depth-from adds the depth step to the barred s2 only, not {branch} {submission}")
        spec, info = orig_resolve(dec, rules, branch, submission, hashes)
        if set(spec) != {"factor", "composition"} or spec["factor"] != FROZEN_FACTOR:
            raise Refused("the barred s2's spec is not C2's (its factor alone and a composition, no depth, stage or "
                          "coverage step): not built")
        out.update(replaced=copy.deepcopy(spec), n_resolves=out.get("n_resolves", 0) + 1)
        return {**spec, "depth": True}, {**info, "depth_source": DEPTH_SOURCE}

    def ratios(V, p):
        read = p.get("depth_ratio_pooled") if isinstance(p, dict) else None
        if not isinstance(read, dict) or read != params["depth_ratio_pooled"]:
            raise Refused("the depth ratios the build reads are not the pinned decision file's: not built")
        ratio, rec = depth_ratios_applied(V, p, exponent, form, reader=orig_ratios)
        out.update(ratio=ratio, record=rec, n_ratios=out.get("n_ratios", 0) + 1)
        return ratio

    def construct(ref, n, rng, f, phi=None, w=1.0, ratio=None, alphas=None):
        if ratio is None or ratio is not out.get("ratio") or alphas is not None:
            raise DepthCheckFailed("construct did not receive the line's ratios alone: nothing written")
        lognormed: list = []
        with _observe("lognorm", lognormed, _lognorm_spy):
            c2 = orig_construct(ref, n, _twin_generator(rng), f, phi=phi, w=w, ratio=None, alphas=None)
        calls: list = []
        with _observe("rescale_depth", calls, _rescale_spy):
            built = orig_construct(ref, n, rng, f, phi=phi, w=w, ratio=ratio, alphas=None)
        if len(lognormed) != 1 or len(calls) != 1:
            raise DepthCheckFailed("C2's cells or the depth step were not observed exactly once: nothing written")
        checks = depth_checks(c2, lognormed[0], built, calls[0], ratio)
        out.update(checks=checks, n_constructs=out.get("n_constructs", 0) + 1)
        if not checks["passed"]:
            bad = [k for k, v in checks.items() if v is False]
            raise DepthCheckFailed(f"the depth step's checks failed ({', '.join(bad)}): nothing written")
        return built

    bs.resolve_spec, bs.depth_ratios, bs.construct = resolve, ratios, construct
    try:
        yield
    finally:
        bs.resolve_spec, bs.depth_ratios, bs.construct = orig_resolve, orig_ratios, orig_construct


# ------------------------------------------ row R76: T3's refinement round on C2's command, five files (off by default)

# Row R76 (the captain's order of 2026-11-26, ruling, no run: his named exception on T3's depth for each of the five
# files; for T3_c2_s2_g15, T3_c2_s2_g05 and T3_c2_s2_g2_s067 his named exception for the Gata4 weight, as R41's for
# T3_c2_s2_g2; for T3_c2_s2_s085, T3_c2_s2_s05 and T3_c2_s2_g2_s067 "the 5% check gates nothing for it"): T3's
# refinement round on C2's command (the T3 REF3 design, as its cold challenge fixes it; the design's s11 struck). One
# literal table: each file, its option string as row R76 pairs them, and whether the row lifts the 5% den check for it.
# The builder takes exactly these five strings, each only with --captain-ruling R76, with --den-ruling R76 for the
# three the row lifts and never for the two weight-only strings. The R76 check runs first and only with
# --captain-ruling R76: every other command goes on exactly as before. Off by default; its own code below, beside the
# other lines, which it leaves unchanged.
R76_ROW = "R76"
R76_FORM = "c2-s2"
R76_TABLE = (
    ("T3_c2_s2_s085", "--form c2-s2 --size-scale 0.85", True),
    ("T3_c2_s2_s05", "--form c2-s2 --size-scale 0.5", True),
    ("T3_c2_s2_g15", "--form c2-s2 --gata4-exponent 1.5", False),
    ("T3_c2_s2_g05", "--form c2-s2 --gata4-exponent 0.5", False),
    ("T3_c2_s2_g2_s067", "--form c2-s2 --gata4-exponent 2.0 --size-scale 0.67", True),
)
R76_FLAGS = (("--gata4-exponent", "g"), ("--size-scale", "s"))  # the options a string may give, in this order
R76_DEN_OPTION = "--den-ruling R76"  # what the row's lifting part names, beside the lifted strings


def _r76_tag(value: str) -> str:
    """A value's part of a file name, as the record names T3's files: its point dropped, and a trailing .0 with it
    (2.0 -> 2, 1.5 -> 15, 0.85 -> 085)."""
    return value[:-2] if value.endswith(".0") else value.replace(".", "")


def _r76_values(name: str, option: str) -> tuple:
    """A row R76 option string's --gata4-exponent and --size-scale values, None for an option it does not give, read
    off the string: `--form c2-s2`, then one or both of R76_FLAGS in that order, each with a decimal value; and the
    file's name must be T3_c2_s2_ and the string's tags (g and s, each with its value's tag). ValueError otherwise."""
    words = option.split()
    flags = [flag for flag, _ in R76_FLAGS]
    tail = words[2:]
    given = tail[::2]
    if (words[:2] != ["--form", R76_FORM] or not tail or len(tail) % 2 or any(f not in flags for f in given)
            or [flags.index(f) for f in given] != sorted({flags.index(f) for f in given})
            or not all(re.fullmatch(r"[0-9]+\.[0-9]+", v) for v in tail[1::2])):
        raise ValueError(f"not a row {R76_ROW} option string: {option}")
    values = dict(zip(given, tail[1::2], strict=True))
    tags = [f"{tag}{_r76_tag(values[flag])}" for flag, tag in R76_FLAGS if flag in values]
    if name != "T3_c2_s2_" + "_".join(tags):
        raise ValueError(f"row {R76_ROW}'s file {name} does not match its option string {option}")
    return tuple(float(values[flag]) if flag in values else None for flag, _ in R76_FLAGS)


# file (without .h5ad, as row R76 names it) -> (its option string, the --gata4-exponent value or None, the
# --size-scale value or None, whether the row lifts the 5% den check for it), read off the literal table
R76_FILES = {name: (option, *_r76_values(name, option), lifted) for name, option, lifted in R76_TABLE}
R76_LINE = ("T3's refinement round on C2's command (row R76): the family ln phi(g, s) = s c_g (ln phi_M + g ln phi_G), "
            "c_g = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), on the barred s2, at row R76's five points: "
            "T3_c2_s2_s085 (g = 1, s = 0.85), T3_c2_s2_s05 (g = 1, s = 0.5), T3_c2_s2_g15 (g = 1.5), T3_c2_s2_g05 "
            "(g = 0.5) and T3_c2_s2_g2_s067 (g = 2, s = 0.67); the size-held part line 3's c2_gata4_weighted at g")
R76_LIST_ENTRY = ("experiments/harness/TRIES.md row R76 (the captain's order of 2026-11-26: his named exception on "
                  "T3's depth for each of the five files; for T3_c2_s2_g15, T3_c2_s2_g05 and T3_c2_s2_g2_s067 his "
                  "named exception for the Gata4 weight, as R41's; for T3_c2_s2_s085, T3_c2_s2_s05 and "
                  "T3_c2_s2_g2_s067 the 5% check gates nothing; decision log 2026-11-26 (a)); the T3 REF3 design of "
                  "2026-11-25 and its cold challenge (section 5, the builder plan)")
R76_CLOSED_LEVER = ("T3's depth (closed after an exposure): the factor moves each label's expected library and the "
                    "cells' totals against C2's; built only under the captain's named exception on T3's depth for "
                    "this file (row R76)")
R76_DEN_RULE = ("line 3's check: the den-weighted total of the file's factor lies within 5% of C2's factor's (den from "
                "released WT E8.75 as build_outside_form.wt_den makes it, the two pinned tables only); it gates the "
                "two weight-only files and is recorded, gating nothing, for the three files row R76 lifts it for "
                "(--den-ruling R76); the shifts from phi_M's total are recorded and gate nothing")


class R76GuardFailed(GonlyGuardFailed):
    """Row R76's files: the range guard or the 5% den check (not lifted for the file) fails (exit status 6, nothing
    written)."""


def r76_name(gata4_exponent=None, size_scale=None) -> str:
    """The file of row R76 whose option string a command gives exactly: its --gata4-exponent and --size-scale values,
    None for an option not given, read as the options read every value (float(); never a bool), equal to one row of
    R76_FILES, the options given and not given included; Refused otherwise (`--gata4-exponent 1.0 --size-scale 0.85`
    is not `--size-scale 0.85`, and 0.850001 is not 0.85)."""
    got, said = [], []
    for flag, v in (("--gata4-exponent", gata4_exponent), ("--size-scale", size_scale)):
        if v is None:
            got.append(None)
            continue
        said.append(f"{flag} {v!r}")
        try:
            if isinstance(v, bool):
                raise TypeError
            got.append(float(v))
        except (TypeError, ValueError, OverflowError):
            raise Refused(f"{flag} {v!r} is refused under --captain-ruling {R76_ROW}: not one of row {R76_ROW}'s "
                          "option strings") from None
    names = [name for name, (_, g, s, _) in R76_FILES.items() if (g, s) == tuple(got)]
    if len(names) != 1:
        raise Refused(f"--form {R76_FORM} {' '.join(said) or '(no option)'} --captain-ruling {R76_ROW} is refused: not "
                      f"one of row {R76_ROW}'s five option strings exactly ("
                      + "; ".join(f"`{opt}`" for opt, *_ in R76_FILES.values()) + ")")
    return names[0]


def r76_lift_part(first: str | None) -> list[str] | None:
    """The whole backtick-delimited strings of the one part of a first cell (its text between semicolons) that says
    GONLY_LIFT_WORDS ("the 5% check gates nothing for it"); None when no part or more than one says it."""
    if first is None:
        return None
    parts = [part for part in first.split(";") if GONLY_LIFT_WORDS in part]
    return re.findall(r"`([^`]*)`", parts[0]) if len(parts) == 1 else None


def r76_row_holds_the_table(first: str | None) -> bool:
    """Whether row R76's first cell holds the literal table: each file's option string as one whole backtick-delimited
    string, followed there by " (<file>)" as the row pairs them, and one lifting part (r76_lift_part) naming
    `--den-ruling R76` and exactly the strings the table lifts, as whole strings (a part naming another string, or
    missing one, fails: check_den_lift's reading, the words anywhere in the cell, would lift every string)."""
    if first is None:
        return False
    for name, (opt, *_) in R76_FILES.items():
        if not (names_whole(first, opt) and f"`{opt}` ({name})" in first):
            return False
    part = r76_lift_part(first)
    lifted = {opt for opt, _, _, lift in R76_FILES.values() if lift}
    return part is not None and R76_DEN_OPTION in part and set(part) - {R76_DEN_OPTION} == lifted


def check_r76(form, gata4_exponent=None, size_scale=None, captain_ruling=None, tries: Path | None = None,
              composition_from=None, den_ruling=None, depth_from=None, depth_exponent=None, depth_form=None
              ) -> dict | None:
    """--captain-ruling R76 (row R76's five files): None unless the ruling is R76 (nothing is read; every other command
    goes on to the other lines' checks as before); else the file's record, or Refused: --composition-from or a
    --depth-* option beside it; any form but c2-s2; options that are not exactly one of the five strings (r76_name);
    for a string the table lifts, a --den-ruling other than R76, or none; for a weight-only string, any --den-ruling; a
    ledger whose row R76 does not pass build_outside_form.captain_order_names for the file's string or does not hold
    the table (r76_row_holds_the_table)."""
    if captain_ruling != R76_ROW:
        return None
    for flag, value in (("--composition-from", composition_from), ("--depth-from", depth_from),
                        ("--depth-exponent", depth_exponent), ("--depth-form", depth_form)):
        if value is not None:
            raise Refused(f"--captain-ruling {R76_ROW} takes no {flag} (each of its files is C2's command with its own "
                          "option string)")
    if form != R76_FORM:
        raise Refused(f"--captain-ruling {R76_ROW} is taken on --form {R76_FORM} only, not {form!r}")
    name = r76_name(gata4_exponent, size_scale)
    opt, g, s, lifted = R76_FILES[name]
    if lifted and den_ruling != R76_ROW:
        raise Refused(f"`{opt}` --captain-ruling {R76_ROW} is taken only with --den-ruling {R76_ROW} (row {R76_ROW}: "
                      f"the 5% check gates nothing for it), not with {den_ruling!r}")
    if not lifted and den_ruling is not None:
        raise Refused(f"--den-ruling {den_ruling!r} is refused for `{opt}`: row {R76_ROW} lifts the 5% check for three "
                      "strings only, not this one, so the check gates this file")
    p = TRIES_FILE if tries is None else Path(tries)
    first = ledger_first_cell(R76_ROW, p)
    if not (_outside_form().captain_order_names(R76_ROW, f"`{opt}`", p) and r76_row_holds_the_table(first)):
        raise Refused(f"`{opt}` --captain-ruling {R76_ROW} is refused: experiments/harness/TRIES.md has no single row "
                      f"{R76_ROW} whose first cell is the captain's order (\"ruling, no run\") naming "
                      f"`--captain-ruling {R76_ROW}` and each of the table's five option strings as one whole string "
                      f"followed by its file's name, with one part saying \"{GONLY_LIFT_WORDS}\" that names "
                      f"`{R76_DEN_OPTION}` and exactly the three strings the table lifts")
    sha = bs.sha256(p)
    lift = None if not lifted else {"row": R76_ROW, "option": f"`{opt}`", "words": GONLY_LIFT_WORDS,
                                    "ledger": "experiments/harness/TRIES.md", "ledger_sha256": sha}
    return {"row": R76_ROW, "name": name, "file": f"{name}.h5ad", "option": f"`{opt}`",
            "g": 1.0 if g is None else g, "s": 1.0 if s is None else s, "den_ruling": lift,
            "ledger": "experiments/harness/TRIES.md", "ledger_sha256": sha}


def r76_factor(phi_m, phi_g, g: float, s: float = 1.0) -> tuple[np.ndarray, dict]:
    """Row R76's factor at (g, s), from the two tables as loaded, gene by gene: line 3's size-held product
    c2_gata4_weighted(phi_M, phi_G, g), value for value (g = 1 gives combine('c2', 1.0)'s array and g = 2
    T3_c2_s2_g2's), then, unless s is 1, ln phi multiplied by s (** s); nothing rescaled or clipped after. Returns
    (phi, record); g must be finite (line 3's function refuses any other) and s finite and positive."""
    s = float(s)
    if not (np.isfinite(s) and s > 0):
        raise Refused(f"the size scale {s!r} is not finite and positive")
    held, hrec = c2_gata4_weighted(phi_m, phi_g, g)
    phi = held if s == 1.0 else held ** s
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused("the combined factor is not finite and positive")
    rec = {"kind": "c2", "x": 1.0, "g": hrec["g"], "size_scale": s, "gata4_weight": s * hrec["gata4_weight"],
           "mab21l2_exponent": s * hrec["mab21l2_exponent"], "c": hrec["c"],
           "c_rule": "line 3's c = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), population SD over the genes "
                     "(c2_gata4_weighted at g, value for value); then ln phi multiplied by the size scale s (** s, "
                     "none at s = 1); the powers applied are s c (transfer_mab21l2_all) and s g c (gata4_pooled)",
           "sd_ln_phi_m": hrec["sd_ln_phi_m"], "sd_ln_phi_g": hrec["sd_ln_phi_g"], "sd_ln_size_held": _sd(np.log(held)),
           "base_factor": hrec["base_factor"], "size_held_factor": hrec["factor"], "factor": _shift(phi),
           "den_weighted_mean": "recorded in combined_build.r76_line.guards",
           "min": float(phi.min()), "max": float(phi.max()), "n_up": int((phi > 1).sum()),
           "n_down": int((phi < 1).sum())}
    return phi, rec


def r76_guards(phi_m, phi_g, g: float, s: float, den, lift: dict | None = None) -> dict:
    """The two checks of the options on C2's command (gonly_guards' checks and figures) for row R76's factor
    (r76_factor), a pure function of the two tables and den; no cell is drawn. Range: the file's max at most
    GONLY_RANGE_FACTOR times C2's max and its min at least C2's min over GONLY_RANGE_FACTOR. Den, line 3's check: the
    file's den-weighted total (den . phi / sum den) within C2G_MAX_DEN_SHIFT of C2's; it gates unless `lift` (row R76's
    lifting of the file, check_r76's record) is given, when it is recorded and gates nothing. The shifts from phi_M's
    total are recorded and gate nothing. Also the sha256 of the factor described."""
    phi, rec = r76_factor(phi_m, phi_g, g, s)
    c2, _ = combine("c2", 1.0, phi_m, phi_g)
    phi_m = np.asarray(phi_m, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    if den.shape != phi_m.shape or not np.all(np.isfinite(den) & (den >= 0)) or not den.sum() > 0:
        raise Refused("den does not cover the tables' genes with finite, non-negative weights and a positive sum")
    lo, hi = float(c2.min()) / GONLY_RANGE_FACTOR, float(c2.max()) * GONLY_RANGE_FACTOR
    fmin, fmax = float(phi.min()), float(phi.max())
    range_ok = bool(fmin >= lo and fmax <= hi)
    tot = float(den.sum())
    t_m, t_c2, t_f = (float(den @ f / tot) for f in (phi_m, c2, phi))
    shift = t_f / t_c2 - 1.0
    within = bool(abs(shift) <= C2G_MAX_DEN_SHIFT)
    den_ok = bool(within or lift is not None)
    status = "gates" if lift is None else f"lifted by row {lift['row']}: {GONLY_LIFT_WORDS}"
    return {"g": rec["g"], "size_scale": rec["size_scale"], "c": rec["c"],
            "powers_applied": {TABLE_M: rec["mab21l2_exponent"], TABLE_G: rec["gata4_weight"]},
            "range": {"factor_min": fmin, "factor_max": fmax, "c2_min": float(c2.min()), "c2_max": float(c2.max()),
                      "bound_min": lo, "bound_max": hi, "factor": GONLY_RANGE_FACTOR, "passed": range_ok,
                      "rule": GONLY_RANGE_RULE},
            "den_check": {"den_total_phi_m": t_m, "den_total_c2": t_c2, "den_total_factor": t_f,
                          "den_shift_vs_c2": shift, "max_den_shift": C2G_MAX_DEN_SHIFT, "within": within,
                          "limit": {"max_abs_den_shift": C2G_MAX_DEN_SHIFT, "status": status}, "lifted_by": lift,
                          "passed": den_ok, "den_shift_c2_vs_phi_m": t_c2 / t_m - 1.0,
                          "den_shift_factor_vs_phi_m": t_f / t_m - 1.0, "vs_phi_m": "recorded; gates nothing",
                          "rule": R76_DEN_RULE},
            "passed": bool(range_ok and den_ok), "factor_sha256": factor_sha256(phi)}


def run_r76_guards(dec: dict, wt: Path, run_dir: Path | None, g: float, s: float, lift: dict | None = None) -> dict:
    """The checks from released WT E8.75's den and the two pinned tables only, before any build: den and genes by
    build_outside_form.wt_den, both tables through the frozen loader in that gene order (checked against TABLE_SHA256).
    Returns the record whatever it says; the caller refuses on it."""
    genes, den = _outside_form().wt_den(Path(wt))
    phi_m, rec_m = bs.load_factor(run_dir, TABLE_M, genes, dec)
    phi_g, rec_g = bs.load_factor(run_dir, TABLE_G, genes, dec)
    for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
        if rec["sha256"] != TABLE_SHA256[nm]:
            raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                          f"({TABLE_SHA256[nm][:16]}): not built")
    den = np.asarray(den, dtype=np.float64)
    out = r76_guards(phi_m, phi_g, g, s, den, lift)
    out |= {"wt": str(wt), "wt_sha256": bs.sha256(Path(wt)), "den_rule": G2_DEN_RULE, "den": den.tolist(),
            "den_sha256": factor_sha256(den), "n_genes": len(genes),
            "genes_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest(),
            "tables_sha256": {TABLE_M: rec_m["sha256"], TABLE_G: rec_g["sha256"]}, "genes": list(genes)}
    return out


@contextmanager
def r76_loader(g: float, s: float, out: dict, guards: dict):
    """bs.load_factor, for this process only, on row R76's line: the frozen loader for phi_M and phi_G (checked against
    TABLE_SHA256 too), then r76_factor at (g, s), the hook C2's factor goes through. FactorMismatch (exit 5) unless the
    build's genes are the guards' and the factor's sha256 is the one the guards describe."""
    original = bs.load_factor

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        if list(genes) != guards["genes"]:
            raise FactorMismatch("the build's genes are not those the guards were computed over: not built")
        phi_m, rec_m = original(run_dir, TABLE_M, genes, dec)
        phi_g, rec_g = original(run_dir, TABLE_G, genes, dec)
        for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
            if rec["sha256"] != TABLE_SHA256[nm]:
                raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                              f"({TABLE_SHA256[nm][:16]}): not built")
        phi, crec = r76_factor(phi_m, phi_g, g, s)
        sha = factor_sha256(phi)
        if sha != guards["factor_sha256"]:
            raise FactorMismatch(f"the factor the build would apply ({sha[:16]}) is not the one the guards describe "
                                 f"({guards['factor_sha256'][:16]}): not built")
        out.update(tables={TABLE_M: rec_m, TABLE_G: rec_g}, combination=crec, factor_sha256=sha,
                   n_loads=out.get("n_loads", 0) + 1)
        return phi, rec_m

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


def build_r76(decision: Path, base: str, wt: Path, board, n_cells: int, out: Path, r76: dict,
              run_dir: Path | None = None, seed: int = 1_000_001, overwrite: bool = False, threads: int = 16) -> dict:
    """Row R76's file (check_r76's record) on C2's command through bs.build_submission: the frozen and clean checks
    and the frozen builder's own, then the two checks before any cell is drawn (R76GuardFailed, exit 6, nothing
    written), then the build with r76_loader; returns the sidecar (written next to the file, the frozen builder's
    plus `combined_build` with its `r76_line`)."""
    kind = "c2"
    slot, x, base_sha = BASES[base]
    frozen = check_frozen()
    tree = check_clean()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)  # the frozen builder's own check
    spec = bs.barred_table(dec, rules)[slot]
    if spec["factor"] != FROZEN_FACTOR:
        raise Refused(f"the frozen barred {slot} factor is {spec['factor']}, not {FROZEN_FACTOR}: not built")
    out_p = Path(out)
    if not overwrite and (out_p.exists() or out_p.with_suffix(".json").exists()):
        raise Refused(f"{out_p} or its sidecar exists; pass --overwrite")
    g, s = r76["g"], r76["s"]
    guards = run_r76_guards(dec, Path(wt), run_dir, g, s, r76["den_ruling"])
    if not guards["passed"]:
        raise R76GuardFailed(gonly_failure(guards))
    got: dict = {}
    with r76_loader(g, s, got, guards):
        side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                   run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
    if got.get("factor_sha256") != guards["factor_sha256"]:
        raise RuntimeError("the factor applied is not the one the guards describe")
    if got.get("n_loads") != 1 or side["construction"]["factor_exponent"] != 1.0 or \
            side["spec"]["factor"] != FROZEN_FACTOR:
        raise RuntimeError("the build did not apply the combined factor once at w = 1")
    if (set(side["spec"]) != {"factor", "composition"} or side["construction"]["depth_ratios"] is not None
            or side["construction"]["stage_alphas"] is not None):
        raise RuntimeError("the build took more than C2's command with its factor")
    side["inputs"]["factor_note"] = FACTOR_NOTE
    side["combined_build"] = {
        "script": SCRIPT, "script_sha256": bs.sha256(SCRIPT_PATH), "form": f"{kind}-{base}",
        "base": {"name": base, "slot": slot, "x": x, "reference_sha256": base_sha,
                 "rng_streams": side["inputs"]["rng_streams"]},
        "tables": {k: {"file": v["file"], "sha256": v["sha256"], "recorded_in": v["recorded_in"],
                       "pinned_sha256": TABLE_SHA256[k]} for k, v in got["tables"].items()},
        "combination": got["combination"], "where": WHERE, "list_entry": LIST_ENTRY,
        "frozen_files_sha256": frozen, "tree_at_start": tree, "identity": None, "depth_step": None,
    }
    crec = got["combination"]
    side["combined_build"]["r76_line"] = {
        "captain_ruling": r76, "line": R76_LINE, "file": r76["file"], "option": r76["option"], "g": crec["g"],
        "size_scale": crec["size_scale"], "c": crec["c"],
        "powers_applied": {TABLE_M: crec["mab21l2_exponent"], TABLE_G: crec["gata4_weight"]},
        "option_value": "on --form c2-s2 the value of --gata4-exponent is the gata4_pooled weight g, not a power "
                        "(g = 1, C2's, when absent); --size-scale s multiplies ln phi after the size-held combination "
                        "(s = 1 when absent)",
        "factor_rule": "ln phi = s c_g (ln phi_M + g ln phi_G), c_g = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), "
                       "population SD over the genes (line 3's c2_gata4_weighted at g, then ** s); nothing rescaled or "
                       "clipped after",
        "closed_lever": R76_CLOSED_LEVER, "guards": {k: v for k, v in guards.items() if k != "genes"},
        "den_ruling": r76["den_ruling"], "factor_sha256_applied": got["factor_sha256"],
        "c2_reference_sha256": C2_S2_SHA256, "list_entry": R76_LIST_ENTRY,
    }
    bs.dump_json(side, Path(out).with_suffix(".json"))
    return side


# ---------------------------------- row R79: T3's refinement round 4 on C2's command, five files (off by default)

# Row R79 (the captain's order of 2026-11-26, ruling, no run: the captain's named exception on T3's depth for each of
# the five files; for T3_c2_s2_g05_s05, T3_c2_s2_g2_s05 and T3_c2_s2_g15_s05 the named exception for the Gata4 weight,
# as R41's for T3_c2_s2_g2; for T3_c2_s2_gonly_s05 the named exception for the Gata4 weight, as R47's for
# T3_c2_s2_gonly; for each of the five strings "the 5% check gates nothing for it"): T3's refinement round 4 on C2's
# command (the T3 REF4 design with its cold challenge's fixes; all five sent at once). One literal table: each file,
# its option string as row R79 pairs them, and whether the row lifts the 5% den check for it (all five: each string
# carries a size scale).
# The builder takes exactly these five strings, each only with --captain-ruling R79 and --den-ruling R79. The R79 check
# runs first, before R76's, and only with --captain-ruling R79: every other command goes on exactly as before. At a
# finite g the file goes through row R76's functions (r76_factor, run_r76_guards, r76_loader), at g = inf through the
# options' g = inf branch (c2_gonly_factor, run_gonly_guards, gonly_loader); row R79 computes no factor and no check of
# its own. Off by default; its own code below, beside the other lines, which it leaves unchanged.
R79_ROW = "R79"
R79_FORM = "c2-s2"
R79_TABLE = (
    ("T3_c2_s2_s035", "--form c2-s2 --size-scale 0.35", True),
    ("T3_c2_s2_g05_s05", "--form c2-s2 --gata4-exponent 0.5 --size-scale 0.5", True),
    ("T3_c2_s2_gonly_s05", "--form c2-s2 --gata4-exponent inf --size-scale 0.5", True),
    ("T3_c2_s2_g2_s05", "--form c2-s2 --gata4-exponent 2.0 --size-scale 0.5", True),
    ("T3_c2_s2_g15_s05", "--form c2-s2 --gata4-exponent 1.5 --size-scale 0.5", True),
)
R79_DEN_OPTION = "--den-ruling R79"  # what the row's lifting part names, beside the lifted strings
R79_INF = "inf"  # the one --gata4-exponent value of a row R79 string that is not a decimal: exactly these three letters
R79_INF_TAG = "only"  # its part of the file's name after the g, as row R47 names T3_c2_s2_gonly


def _r79_values(name: str, option: str) -> tuple:
    """A row R79 option string's --gata4-exponent and --size-scale values, None for an option it does not give, read
    off the string as _r76_values reads row R76's (`--form c2-s2`, single spaces, then one or both of R76_FLAGS in that
    order, each with a decimal value; the file's name T3_c2_s2_ and the string's tags, _r76_tag's), with one value more:
    --gata4-exponent R79_INF, exactly the letters "inf" (not "Inf", "INF", "+inf" or "infinity", and never on
    --size-scale), read as math.inf and tagged "gonly" (as row R47 names T3_c2_s2_gonly). ValueError otherwise."""
    words = option.split()
    flags = [flag for flag, _ in R76_FLAGS]
    tail = words[2:]
    given = tail[::2]
    if (option != " ".join(words) or words[:2] != ["--form", R79_FORM] or not tail or len(tail) % 2
            or any(f not in flags for f in given)
            or [flags.index(f) for f in given] != sorted({flags.index(f) for f in given})
            or not all(re.fullmatch(r"[0-9]+\.[0-9]+", v) or (f == "--gata4-exponent" and v == R79_INF)
                       for f, v in zip(given, tail[1::2], strict=True))):
        raise ValueError(f"not a row {R79_ROW} option string: {option}")
    values = dict(zip(given, tail[1::2], strict=True))
    tags = [f"{tag}{R79_INF_TAG if values[flag] == R79_INF else _r76_tag(values[flag])}" for flag, tag in R76_FLAGS
            if flag in values]
    if name != "T3_c2_s2_" + "_".join(tags):
        raise ValueError(f"row {R79_ROW}'s file {name} does not match its option string {option}")
    return tuple((math.inf if values[flag] == R79_INF else float(values[flag])) if flag in values else None
                 for flag, _ in R76_FLAGS)


# file (without .h5ad, as row R79 names it) -> (its option string, the --gata4-exponent value or None (math.inf for
# "inf"), the --size-scale value or None, whether the row lifts the 5% den check for it), read off the literal table
R79_FILES = {name: (option, *_r79_values(name, option), lifted) for name, option, lifted in R79_TABLE}
R79_LINE = ("T3's refinement round 4 on C2's command (row R79): the family ln phi(g, s) = s c_g (ln phi_M + g ln "
            "phi_G), c_g = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), on the barred s2, at row R79's five points: "
            "T3_c2_s2_s035 (g = 1, s = 0.35), T3_c2_s2_g05_s05 (g = 0.5, s = 0.5), T3_c2_s2_gonly_s05 (g = inf: "
            "ln phi = s (SD(ln phi_M) / SD(ln phi_G)) ln phi_G, phi_M dropped; s = 0.5), T3_c2_s2_g2_s05 (g = 2, "
            "s = 0.5) and T3_c2_s2_g15_s05 (g = 1.5, s = 0.5); at a finite g row R76's r76_factor, at g = inf the "
            "options' c2_gonly_factor, each value for value")
R79_LIST_ENTRY = ("experiments/harness/TRIES.md row R79 (the captain's order of 2026-11-26: the named exception on "
                  "T3's depth for each of the five files; for T3_c2_s2_g05_s05, T3_c2_s2_g2_s05 and T3_c2_s2_g15_s05 "
                  "the named exception for the Gata4 weight, as R41's for T3_c2_s2_g2; for T3_c2_s2_gonly_s05 the "
                  "named exception for the Gata4 weight, as R47's for T3_c2_s2_gonly; for each of the five strings the "
                  "5% check gates nothing; decision log 2026-11-26 (h)); the T3 REF4 design of 2026-11-26 and its cold "
                  "challenge (sections 5 and 7)")
R79_CLOSED_LEVER = ("T3's depth (closed after an exposure): the factor moves each label's expected library and the "
                    "cells' totals against C2's; built only under the captain's named exception on T3's depth for "
                    "this file (row R79); for the four files that move g off C2's 1, also the named exception for the "
                    "Gata4 weight (as R41's at g = 0.5, 2 and 1.5; as R47's at g = inf)")
R79_DEN_RULE = ("line 3's check: the den-weighted total of the file's factor lies within 5% of C2's factor's (den from "
                "released WT E8.75 as build_outside_form.wt_den makes it, the two pinned tables only); row R79 lifts "
                "it for each of its five files (--den-ruling R79), so its figures are recorded and gate nothing; the "
                "range guard gates; the figures are r76_guards' at a finite g and gonly_guards' at g = inf, as "
                "computed (the den_check's own rule text is that function's)")


class R79GuardFailed(GonlyGuardFailed):
    """Row R79's files: the range guard fails (the 5% den check is lifted for each of the five: recorded, gating
    nothing) (exit status 6, nothing written)."""


def r79_name(gata4_exponent=None, size_scale=None) -> str:
    """The file of row R79 whose option string a command gives exactly: its --gata4-exponent and --size-scale values,
    None for an option not given, read as the options read every value (float(), so the command line's inf is
    math.inf; never a bool), equal to one row of R79_FILES, the options given and not given included; Refused
    otherwise (`--gata4-exponent 1.0 --size-scale 0.35` is not `--size-scale 0.35`, 0.350001 is not 0.35, and -inf
    or nan is no weight of the table)."""
    got, said = [], []
    for flag, v in (("--gata4-exponent", gata4_exponent), ("--size-scale", size_scale)):
        if v is None:
            got.append(None)
            continue
        said.append(f"{flag} {v!r}")
        try:
            if isinstance(v, bool):
                raise TypeError
            got.append(float(v))
        except (TypeError, ValueError, OverflowError):
            raise Refused(f"{flag} {v!r} is refused under --captain-ruling {R79_ROW}: not one of row {R79_ROW}'s "
                          "option strings") from None
    names = [name for name, (_, g, s, _) in R79_FILES.items() if (g, s) == tuple(got)]
    if len(names) != 1:
        raise Refused(f"--form {R79_FORM} {' '.join(said) or '(no option)'} --captain-ruling {R79_ROW} is refused: not "
                      f"one of row {R79_ROW}'s five option strings exactly ("
                      + "; ".join(f"`{opt}`" for opt, *_ in R79_FILES.values()) + ")")
    return names[0]


def r79_row_holds_the_table(first: str | None) -> bool:
    """Whether row R79's first cell holds the literal table, as r76_row_holds_the_table reads row R76's: each file's
    option string as one whole backtick-delimited string, followed there by " (<file>)" as the row pairs them, and one
    lifting part (r76_lift_part: exactly one part between semicolons says "the 5% check gates nothing for it"; the
    row's quoted words, "the 5% check gates nothing.", are no such part) naming `--den-ruling R79` and exactly the
    strings the table lifts, all five, as whole strings."""
    if first is None:
        return False
    for name, (opt, *_) in R79_FILES.items():
        if not (names_whole(first, opt) and f"`{opt}` ({name})" in first):
            return False
    part = r76_lift_part(first)
    lifted = {opt for opt, _, _, lift in R79_FILES.values() if lift}
    return part is not None and R79_DEN_OPTION in part and set(part) - {R79_DEN_OPTION} == lifted


def check_r79(form, gata4_exponent=None, size_scale=None, captain_ruling=None, tries: Path | None = None,
              composition_from=None, den_ruling=None, depth_from=None, depth_exponent=None, depth_form=None
              ) -> dict | None:
    """--captain-ruling R79 (row R79's five files): None unless the ruling is R79 (nothing is read; every other command
    goes on to the other lines' checks as before); else the file's record, or Refused: --composition-from or a
    --depth-* option beside it; any form but c2-s2; options that are not exactly one of the five strings (r79_name);
    for a string the table lifts (all five), a --den-ruling other than R79, or none; a ledger whose row R79 does not
    pass build_outside_form.captain_order_names for the file's string or does not hold the table
    (r79_row_holds_the_table). The record's g is the string "inf" at g = inf (as the options' records name it)."""
    if captain_ruling != R79_ROW:
        return None
    for flag, value in (("--composition-from", composition_from), ("--depth-from", depth_from),
                        ("--depth-exponent", depth_exponent), ("--depth-form", depth_form)):
        if value is not None:
            raise Refused(f"--captain-ruling {R79_ROW} takes no {flag} (each of its files is C2's command with its own "
                          "option string)")
    if form != R79_FORM:
        raise Refused(f"--captain-ruling {R79_ROW} is taken on --form {R79_FORM} only, not {form!r}")
    name = r79_name(gata4_exponent, size_scale)
    opt, g, s, lifted = R79_FILES[name]
    if lifted and den_ruling != R79_ROW:
        raise Refused(f"`{opt}` --captain-ruling {R79_ROW} is taken only with --den-ruling {R79_ROW} (row {R79_ROW}: "
                      f"the 5% check gates nothing for it), not with {den_ruling!r}")
    if not lifted and den_ruling is not None:
        raise Refused(f"--den-ruling {den_ruling!r} is refused for `{opt}`: row {R79_ROW} does not lift the 5% check "
                      "for this string, so the check gates this file")
    p = TRIES_FILE if tries is None else Path(tries)
    first = ledger_first_cell(R79_ROW, p)
    if not (_outside_form().captain_order_names(R79_ROW, f"`{opt}`", p) and r79_row_holds_the_table(first)):
        raise Refused(f"`{opt}` --captain-ruling {R79_ROW} is refused: experiments/harness/TRIES.md has no single row "
                      f"{R79_ROW} whose first cell is the captain's order (\"ruling, no run\") naming "
                      f"`--captain-ruling {R79_ROW}` and each of the table's five option strings as one whole string "
                      f"followed by its file's name, with exactly one part saying \"{GONLY_LIFT_WORDS}\", which names "
                      f"`{R79_DEN_OPTION}` and exactly the five strings the table lifts")
    sha = bs.sha256(p)
    lift = None if not lifted else {"row": R79_ROW, "option": f"`{opt}`", "words": GONLY_LIFT_WORDS,
                                    "ledger": "experiments/harness/TRIES.md", "ledger_sha256": sha}
    return {"row": R79_ROW, "name": name, "file": f"{name}.h5ad", "option": f"`{opt}`",
            "g": R79_INF if g == math.inf else (1.0 if g is None else g), "s": 1.0 if s is None else s,
            "den_ruling": lift, "ledger": "experiments/harness/TRIES.md", "ledger_sha256": sha}


def r79_point(rec: dict) -> tuple[float, float]:
    """(g, s) of a record of check_r79, as the literal table holds them: math.inf for "inf", 1.0 for an option the
    string does not give."""
    _, g, s, _ = R79_FILES[rec["name"]]
    return (1.0 if g is None else g, 1.0 if s is None else s)


def r79_functions(g: float) -> dict:
    """The functions row R79's file at weight g goes through, none of them row R79's own: at a finite g row R76's
    (r76_factor, r76_guards, run_r76_guards, r76_loader), at g = inf the options' g = inf branch (c2_gonly_factor,
    gonly_guards, run_gonly_guards, gonly_loader). Each is looked up when called."""
    if float(g) == math.inf:
        return {"branch": "g = inf: c2_gonly_factor, gonly_guards, run_gonly_guards, gonly_loader",
                "factor": c2_gonly_factor, "guards": gonly_guards, "run_guards": run_gonly_guards,
                "loader": gonly_loader}
    return {"branch": "finite g: r76_factor, r76_guards, run_r76_guards, r76_loader", "factor": r76_factor,
            "guards": r76_guards, "run_guards": run_r76_guards, "loader": r76_loader}


def r79_factor(phi_m, phi_g, g: float, s: float) -> tuple[np.ndarray, dict]:
    """Row R79's factor at (g, s), from the two tables as loaded: r79_functions(g)'s factor, the very function and
    its record (r76_factor(phi_M, phi_G, g, s) at a finite g; c2_gonly_factor(phi_M, phi_G, inf, s) at g = inf)."""
    return r79_functions(g)["factor"](phi_m, phi_g, g, s)


def build_r79(decision: Path, base: str, wt: Path, board, n_cells: int, out: Path, r79: dict,
              run_dir: Path | None = None, seed: int = 1_000_001, overwrite: bool = False, threads: int = 16) -> dict:
    """Row R79's file (check_r79's record) on C2's command through bs.build_submission, as build_r76 builds row
    R76's: the frozen and clean checks and the frozen builder's own, then the file's branch's two checks before any
    cell is drawn (R79GuardFailed, exit 6, nothing written), then the build with that branch's loader
    (r79_functions); returns the sidecar (written next to the file, the frozen builder's plus `combined_build` with
    its `r79_line`)."""
    kind = "c2"
    slot, x, base_sha = BASES[base]
    frozen = check_frozen()
    tree = check_clean()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)  # the frozen builder's own check
    spec = bs.barred_table(dec, rules)[slot]
    if spec["factor"] != FROZEN_FACTOR:
        raise Refused(f"the frozen barred {slot} factor is {spec['factor']}, not {FROZEN_FACTOR}: not built")
    out_p = Path(out)
    if not overwrite and (out_p.exists() or out_p.with_suffix(".json").exists()):
        raise Refused(f"{out_p} or its sidecar exists; pass --overwrite")
    g, s = r79_point(r79)
    fns = r79_functions(g)
    guards = fns["run_guards"](dec, Path(wt), run_dir, g, s, r79["den_ruling"])
    if not guards["passed"]:
        raise R79GuardFailed(gonly_failure(guards))
    got: dict = {}
    with fns["loader"](g, s, got, guards):
        side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                   run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
    if got.get("factor_sha256") != guards["factor_sha256"]:
        raise RuntimeError("the factor applied is not the one the guards describe")
    if got.get("n_loads") != 1 or side["construction"]["factor_exponent"] != 1.0 or \
            side["spec"]["factor"] != FROZEN_FACTOR:
        raise RuntimeError("the build did not apply the combined factor once at w = 1")
    if (set(side["spec"]) != {"factor", "composition"} or side["construction"]["depth_ratios"] is not None
            or side["construction"]["stage_alphas"] is not None):
        raise RuntimeError("the build took more than C2's command with its factor")
    side["inputs"]["factor_note"] = FACTOR_NOTE
    side["combined_build"] = {
        "script": SCRIPT, "script_sha256": bs.sha256(SCRIPT_PATH), "form": f"{kind}-{base}",
        "base": {"name": base, "slot": slot, "x": x, "reference_sha256": base_sha,
                 "rng_streams": side["inputs"]["rng_streams"]},
        "tables": {k: {"file": v["file"], "sha256": v["sha256"], "recorded_in": v["recorded_in"],
                       "pinned_sha256": TABLE_SHA256[k]} for k, v in got["tables"].items()},
        "combination": got["combination"], "where": WHERE, "list_entry": LIST_ENTRY,
        "frozen_files_sha256": frozen, "tree_at_start": tree, "identity": None, "depth_step": None,
    }
    crec = got["combination"]
    side["combined_build"]["r79_line"] = {
        "captain_ruling": r79, "line": R79_LINE, "file": r79["file"], "option": r79["option"], "g": crec["g"],
        "size_scale": crec["size_scale"], "c": crec["c"],
        "powers_applied": {TABLE_M: crec["mab21l2_exponent"], TABLE_G: crec["gata4_weight"]},
        "option_value": "on --form c2-s2 the value of --gata4-exponent is the gata4_pooled weight g, not a power "
                        "(g = 1, C2's, when absent; inf keeps phi_G alone at C2's held size); --size-scale s "
                        "multiplies ln phi after the size-held combination",
        "factor_rule": "ln phi = s c_g (ln phi_M + g ln phi_G), c_g = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), "
                       "population SD over the genes (row R76's r76_factor at a finite g); at g = inf ln phi = "
                       "s (SD(ln phi_M) / SD(ln phi_G)) ln phi_G (the options' c2_gonly_factor); nothing rescaled or "
                       "clipped after",
        "functions": fns["branch"], "closed_lever": R79_CLOSED_LEVER, "den_rule": R79_DEN_RULE,
        "guards": {k: v for k, v in guards.items() if k != "genes"}, "den_ruling": r79["den_ruling"],
        "factor_sha256_applied": got["factor_sha256"], "c2_reference_sha256": C2_S2_SHA256,
        "list_entry": R79_LIST_ENTRY,
    }
    bs.dump_json(side, Path(out).with_suffix(".json"))
    return side


def r79_g_text(g) -> str:
    """g as the command line's summary prints it: the string "inf" as it is, a number in the general format (0.5, 2,
    1.5, 1)."""
    return g if isinstance(g, str) else f"{g:g}"


def build_form(decision: Path, form, wt: Path, board, n_cells: int, out: Path, run_dir: Path | None = None,
               seed: int = 1_000_001, overwrite: bool = False, threads: int = 16, gata4_exponent=None,
               captain_ruling=None, tries: Path | None = None, composition_from=None, size_scale=None,
               den_ruling=None, depth_from=None, depth_exponent=None, depth_form=None) -> dict:
    """The named form through bs.build_submission; returns the sidecar (written next to the file, the frozen
    builder's plus `combined_build`). gata4_exponent 2.0 with captain_ruling 'R22' (default None, today's paths
    unchanged) builds row R22's second T3 line's file, c1-s2 at Gata4 weight 2.0, under the ruling row of TRIES.md
    (tries, default TRIES_FILE). composition_from 'allowed-s1' with captain_ruling 'R36' (default None) builds row
    R36's file cD7, c2-s2 with the decision file's allowed s1 composition block, under the ruling row of TRIES.md.
    gata4_exponent 2.0 with captain_ruling 'R22' on c2-s2 builds row R22's third T3 line's file, C2 at gata4_pooled
    weight g = 2.0, after its guard, under the ruling row of TRIES.md. On c2-s2 with captain_ruling 'R22',
    gata4_exponent inf and/or size_scale 0.67 (default None) build F1, F2 or F3 after the range guard and the 5% den
    check, under the ruling row of TRIES.md; den_ruling (default None) names the row that lifts the 5% check for the
    file. On c2-s2, depth_from 'allowed-s1' (default None) with captain_ruling ROW builds the depth line's dd, with
    depth_exponent 0.5 ddh and with depth_form 'global' dg (each default None), under the ruling row ROW of TRIES.md.
    captain_ruling 'R76' (checked first, before every other line) builds row R76's file whose option string
    gata4_exponent and size_scale give exactly, den_ruling 'R76' for the three strings the row lifts the 5% check for.
    captain_ruling 'R79' (checked before R76's and every other line) builds row R79's file whose option string
    gata4_exponent (inf for T3_c2_s2_gonly_s05) and size_scale give exactly, den_ruling 'R79' for each of its five.
    A command with none of these goes on exactly as before."""
    kind, base = parse_form(form)
    r79 = check_r79(f"{kind}-{base}", gata4_exponent, size_scale, captain_ruling, tries, composition_from, den_ruling,
                    depth_from, depth_exponent, depth_form)
    if r79 is not None:  # row R79's file, on its own path; without --captain-ruling R79 nothing of it is read
        return build_r79(decision, base, wt, board, n_cells, out, r79, run_dir=run_dir, seed=seed,
                         overwrite=overwrite, threads=threads)
    r76 = check_r76(f"{kind}-{base}", gata4_exponent, size_scale, captain_ruling, tries, composition_from, den_ruling,
                    depth_from, depth_exponent, depth_form)
    if r76 is not None:  # row R76's file, on its own path; without --captain-ruling R76 nothing of it is read
        return build_r76(decision, base, wt, board, n_cells, out, r76, run_dir=run_dir, seed=seed,
                         overwrite=overwrite, threads=threads)
    depth = check_depth_line(f"{kind}-{base}", depth_from, depth_exponent, depth_form, captain_ruling, tries,
                             composition_from, gata4_exponent, size_scale, den_ruling)
    gonly = None if depth is not None else check_c2_gonly(f"{kind}-{base}", gata4_exponent, size_scale,
                                                          captain_ruling, tries, composition_from, den_ruling)
    cd7 = None if gonly is not None or depth is not None else check_composition_from(
        f"{kind}-{base}", composition_from, captain_ruling, gata4_exponent, tries)
    c2g = None if cd7 is not None or gonly is not None or depth is not None else check_c2_gata4(
        f"{kind}-{base}", gata4_exponent, captain_ruling, tries)
    ruled = None if cd7 is not None or c2g is not None or gonly is not None or depth is not None else check_gata4(
        f"{kind}-{base}", gata4_exponent, captain_ruling, tries)
    slot, x, base_sha = BASES[base]
    frozen = check_frozen()
    tree = check_clean()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)  # the frozen builder's own check
    spec = bs.barred_table(dec, rules)[slot]
    if spec["factor"] != FROZEN_FACTOR:
        raise Refused(f"the frozen barred {slot} factor is {spec['factor']}, not {FROZEN_FACTOR}: not built")
    got: dict = {}
    if depth is not None:
        out_p = Path(out)
        if not overwrite and (out_p.exists() or out_p.with_suffix(".json").exists()):
            raise Refused(f"{out_p} or its sidecar exists; pass --overwrite")
        params, prec = depth_parameters(dec, Path(decision), depth["depth_from"])
        dhook: dict = {}
        with combined_loader(kind, x, got), depth_step(params, depth["exponent"],
                                                        "global" if depth["form"] == "global" else None, dhook):
            side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                       run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
        drec = dhook.get("record") or {}
        want = {t: v for t, v in zip(drec.get("labels", []), drec.get("applied", []), strict=True) if v != 1.0}
        if (dhook.get("n_resolves") != 1 or dhook.get("n_ratios") != 1 or dhook.get("n_constructs") != 1
                or not dhook["checks"]["passed"] or side["construction"]["depth_ratios"] != want):
            raise RuntimeError("the build did not take the depth line's step, checked, exactly once")
        if (side["spec"] != {**dhook["replaced"], "depth": True} or side["construction"]["stage_alphas"] is not None
                or side["inputs"]["parameters_source"] != DEPTH_PARAMS_SOURCE):
            raise RuntimeError("the build took more than the depth step")
    elif cd7 is not None:
        block, brec = composition_block(dec, rules, Path(decision), composition_from)
        hook: dict = {}
        with combined_loader(kind, x, got), r36_composition(block, hook):
            side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                       run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
        if hook.get("n_resolves") != 1 or hook.get("n_compositions") != 1 or side["spec"]["composition"] != block:
            raise RuntimeError("the build did not replace and check the composition block exactly once")
        if set(side["spec"]) != {"factor", "composition"} or side["construction"]["depth_ratios"] is not None:
            raise RuntimeError("the build took more than the composition block")
    elif c2g is not None:
        out_p = Path(out)
        if not overwrite and (out_p.exists() or out_p.with_suffix(".json").exists()):
            raise Refused(f"{out_p} or its sidecar exists; pass --overwrite")
        guard = run_c2g_guard(dec, Path(wt), run_dir, c2g["gata4_weight"])
        if not guard["passed"]:
            raise GuardFailed(f"the file's den-weighted total is {guard['den_shift_vs_c2']:+.2%} against C2's "
                              f"factor's, beyond {C2G_MAX_DEN_SHIFT:.0%}: not built")
        with c2g_loader(c2g["gata4_weight"], got, guard):
            side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                       run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
        if got.get("factor_sha256") != guard["factor_sha256"]:
            raise RuntimeError("the factor applied is not the one the guard describes")
    elif gonly is not None:
        out_p = Path(out)
        if not overwrite and (out_p.exists() or out_p.with_suffix(".json").exists()):
            raise Refused(f"{out_p} or its sidecar exists; pass --overwrite")
        g_s = gonly_point(gonly)
        guards = run_gonly_guards(dec, Path(wt), run_dir, *g_s, gonly["den_ruling"])
        if not guards["passed"]:
            raise GonlyGuardFailed(gonly_failure(guards))
        with gonly_loader(*g_s, got, guards):
            side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                       run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
        if got.get("factor_sha256") != guards["factor_sha256"]:
            raise RuntimeError("the factor applied is not the one the guards describe")
    elif ruled is None:
        with combined_loader(kind, x, got):
            side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                       run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
    else:
        out_p = Path(out)
        if not overwrite and (out_p.exists() or out_p.with_suffix(".json").exists()):
            raise Refused(f"{out_p} or its sidecar exists; pass --overwrite")
        figs = run_gata4_figures(dec, Path(wt), run_dir, ruled["gata4_exponent"])
        with gata4_loader(ruled["gata4_exponent"], got, figs):
            side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                       run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
        if got.get("factor_sha256") != figs["factor_sha256"]:
            raise RuntimeError("the factor applied is not the one the line's figures were computed from")
    if got.get("n_loads") != 1 or side["construction"]["factor_exponent"] != 1.0 or \
            side["spec"]["factor"] != FROZEN_FACTOR:
        raise RuntimeError("the build did not apply the combined factor once at w = 1")
    side["inputs"]["factor_note"] = FACTOR_NOTE
    side["combined_build"] = {
        "script": SCRIPT, "script_sha256": bs.sha256(SCRIPT_PATH), "form": f"{kind}-{base}",
        "base": {"name": base, "slot": slot, "x": x, "reference_sha256": base_sha,
                 "rng_streams": side["inputs"]["rng_streams"]},
        "tables": {k: {"file": v["file"], "sha256": v["sha256"], "recorded_in": v["recorded_in"],
                       "pinned_sha256": TABLE_SHA256[k]} for k, v in got["tables"].items()},
        "combination": got["combination"], "where": WHERE, "list_entry": LIST_ENTRY,
        "frozen_files_sha256": frozen, "tree_at_start": tree,
        "identity": f"the {base} base file byte for byte ({base_sha[:12]})" if kind == "identity" else None,
        "depth_step": None if depth is None else "combined_build.depth_line",
    }
    if depth is not None:
        drec, checks = dhook["record"], dhook["checks"]
        side["combined_build"]["depth_line"] = {
            "captain_ruling": depth, "line": DEPTH_LINE, "file": depth["file"], "arm": depth["arm"],
            "depth_from": depth["depth_from"], "exponent": drec["exponent"], "form": drec["form"],
            "global_ratio": drec["global_ratio"], "rule": drec["rule"], "decision_file": prec,
            "labels": drec["labels"], "labels_sha256": drec["labels_sha256"],
            "n_labels": drec["n_labels"], "n_labels_with_a_pooled_ratio": drec["n_labels_with_a_pooled_ratio"],
            "wt_share_of_labels_left_at_1": drec["wt_share_of_labels_left_at_1"],
            "ratios_frozen": drec["frozen"], "ratios_frozen_sha256": drec["frozen_sha256"],
            "ratios_applied": drec["applied"], "ratios_applied_sha256": drec["applied_sha256"],
            "sha256_rule": drec["sha256_rule"], "where": DEPTH_WHERE,
            "spec": {"replaced": dhook["replaced"], "added": {"depth": True}, "source": DEPTH_SOURCE},
            "checks": {k: v for k, v in checks.items() if k != "counts"},
            "emptied_by_thinning": {**checks["counts"], "handling": DEPTH_EMPTY_RULE},
            "c2_reference_sha256": C2_S2_SHA256,
            "unchanged": "C2's factor (combined_build.combination), rows and factor draws, composition, row number "
                         "(max_cells, read by code), seed, streams and coordinates rule",
            "taken": "the depth ratios only (never the allowed entry's factor, composition, stage or coverage)",
            "list_entry": DEPTH_LIST_ENTRY,
        }
    if ruled is not None:
        side["combined_build"]["gata4_line"] = {
            "captain_ruling": ruled, "line": G2_LINE, "file": G2_FILE, "gata4_weight": ruled["gata4_exponent"],
            "mab21l2_exponent": 1.0, "factor_rule": "phi_M ** 1 * phi_G ** w, nothing rescaled or clipped after the "
            "product", "r30": G2_R30, "figures": {k: v for k, v in figs.items() if k != "genes"},
            "factor_sha256_applied": got["factor_sha256"], "list_entry": G2_LIST_ENTRY,
        }
    if cd7 is not None:
        side["combined_build"]["composition_line"] = {
            "captain_ruling": cd7, "line": R36_LINE, "file": R36_FILE, "composition_from": composition_from,
            "decision_file": brec, "composition": block, "replaced": hook["replaced"],
            "log2_magnitude": R36_LOG2_MAGNITUDE, "ratio": R36_RATIO,
            "ratios_before_renormalising": hook["ratios"], "c2_reference_sha256": C2_S2_SHA256,
            "unchanged": "C2's factor (combined_build.combination), row number (max_cells, read by code), seed, "
                         "streams and coordinates rule",
            "list_entry": R36_LIST_ENTRY,
        }
    if c2g is not None:
        crec = got["combination"]
        side["combined_build"]["c2_gata4_line"] = {
            "captain_ruling": c2g, "line": C2G_LINE, "file": C2G_FILE, "gata4_weight_g": c2g["gata4_weight"],
            "c": crec["c"], "powers_applied": {TABLE_M: crec["mab21l2_exponent"], TABLE_G: crec["gata4_weight"]},
            "option_value": "on --form c2-s2 the value of --gata4-exponent is the gata4_pooled weight g, not a power",
            "factor_rule": "phi_M ** c * phi_G ** (g c), c = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), population SD "
                           "over the genes; nothing rescaled or clipped after the product",
            "guard": {k: v for k, v in guard.items() if k != "genes"},
            "factor_sha256_applied": got["factor_sha256"], "c2_reference_sha256": C2_S2_SHA256,
            "list_entry": C2G_LIST_ENTRY,
        }
    if gonly is not None:
        crec = got["combination"]
        side["combined_build"]["c2_gonly_line"] = {
            "captain_ruling": gonly, "line": GONLY_LINE, "file": gonly["file"], "g": crec["g"],
            "size_scale": crec["size_scale"], "c": crec["c"],
            "powers_applied": {TABLE_M: crec["mab21l2_exponent"], TABLE_G: crec["gata4_weight"]},
            "option_value": "on --form c2-s2 the value of --gata4-exponent is the gata4_pooled weight g, not a power; "
                            "g = inf keeps phi_G alone at C2's held size; --size-scale s multiplies ln phi after the "
                            "size-held combination (s = 1 when absent)",
            "factor_rule": "ln phi = s c_g (ln phi_M + g ln phi_G), c_g = SD(ln phi_M) / SD(ln phi_M + g ln phi_G); at "
                           "g = inf ln phi = s (SD(ln phi_M) / SD(ln phi_G)) ln phi_G; population SD over the genes; "
                           "nothing rescaled or clipped after",
            "guards": {k: v for k, v in guards.items() if k != "genes"}, "den_ruling": gonly["den_ruling"],
            "factor_sha256_applied": got["factor_sha256"], "c2_reference_sha256": C2_S2_SHA256,
            "list_entry": GONLY_LIST_ENTRY,
        }
    bs.dump_json(side, Path(out).with_suffix(".json"))
    return side


@app.callback()
def main() -> None:
    """The barred file with its Mab21l2 factor combined with the Gata4 table, in one of the fixed forms."""


@app.command()
def build(
    decision: Path = typer.Argument(..., help="decision file (experiments/t3_decomp/results/gata4/"
                                    "gata4_decisions.json)"),
    form: str = typer.Option(..., help=f"one of {', '.join(FORMS)}; nothing else"),
    wt: Path = typer.Option(..., help="WT reference (WT E8.75), as for the frozen builder"),
    board: str = typer.Option(..., help="T3 board key in boards/"),
    out: Path = typer.Option(..., help="output .h5ad (sidecar: same name with .json)"),
    run_dir: Path = typer.Option(None, help="run directory holding factors/"),
    seed: int = typer.Option(1_000_001, help="1000001, the list's seed; any other is refused"),
    n_cells: str = typer.Option("max", help="'max', the list's count: the board's max_cells, read in code, never "
                                "printed; any other is refused"),
    boards_dir: Path = typer.Option(None, help="boards directory (default: the repository's boards/)"),
    overwrite: bool = typer.Option(False, help="replace an existing output"),
    threads: int = typer.Option(16, help="BLAS threads while preparing the reference (the frozen default)"),
    gata4_exponent: float = typer.Option(None, help="2.0 only, with --captain-ruling R22: on --form c1-s2 row R22's "
                                         "second T3 line's file (phi_M * phi_G ** 2); on --form c2-s2 its third "
                                         "line's file, the gata4_pooled weight g = 2 in C2's size-held product (a "
                                         "weight, not a power); each while row R22 of TRIES.md is the captain's order "
                                         "naming it. On --form c2-s2 also inf (F1; with --size-scale 0.67, F3): "
                                         "phi_G alone at C2's held size, while row R22 names the file's whole option "
                                         "string. Under --captain-ruling R76 only row R76's: 1.5, 0.5, and 2.0 with "
                                         "--size-scale 0.67"),
    captain_ruling: str = typer.Option(None, help="R22, with --gata4-exponent 2.0 --form c1-s2 or c2-s2; or R36, "
                                       "with --composition-from allowed-s1 --form c2-s2; or R53, with --depth-from "
                                       "allowed-s1 --form c2-s2; or R76, with one of row R76's five option strings "
                                       "(checked first); nothing else"),
    composition_from: str = typer.Option(None, help="allowed-s1 only, with --form c2-s2 --captain-ruling R36: row "
                                         "R36's file cD7 (C2 with the decision file's submissions.allowed.s1."
                                         "composition), while row R36 of TRIES.md is the captain's ruling naming it"),
    size_scale: float = typer.Option(None, help="0.67 only, with --form c2-s2 --captain-ruling R22, alone (F2, C2's "
                                     "ln factor times 0.67) or with --gata4-exponent inf (F3), while row R22 of "
                                     "TRIES.md is the captain's order naming the file's whole option string. Under "
                                     "--captain-ruling R76 only row R76's: 0.85 and 0.5 alone, 0.67 with "
                                     "--gata4-exponent 2.0"),
    den_ruling: str = typer.Option(None, help="a ruling row's id (R and digits), with the options of F1, F2 or F3 "
                                   "only: the 5% den check gates nothing for the file while that row of TRIES.md names "
                                   "the file's whole option string and says \"the 5% check gates nothing for it\". "
                                   "Under --captain-ruling R76, R76 only, with the three strings that row lifts the "
                                   "check for (and never with the two weight-only strings)"),
    depth_from: str = typer.Option(None, help="allowed-s1 only, with --form c2-s2 --captain-ruling R53: the depth "
                                   "line, C2 with the frozen allowed branch's per-label depth step (T3_c2_s2_dd.h5ad), "
                                   "while row R53 of TRIES.md is the captain's order naming the file's whole option "
                                   "string"),
    depth_exponent: float = typer.Option(None, help="0.5 only, with --depth-from allowed-s1: the depth ratios to the "
                                         "power 0.5 (T3_c2_s2_ddh.h5ad)"),
    depth_form: str = typer.Option(None, help="global only, with --depth-from allowed-s1: one depth ratio for every "
                                   "label, the frozen ratios' geometric mean weighted by WT E8.75's cells per label "
                                   "(T3_c2_s2_dg.h5ad)"),
):
    """Build one fixed form."""
    try:
        if seed != LIST_SEED or n_cells != LIST_N_CELLS:
            raise Refused(LIST_INPUTS_REFUSAL)
        parse_form(form)
        if captain_ruling == R79_ROW:  # row R79's check first, before R76's, and only with --captain-ruling R79
            check_r79(form, gata4_exponent, size_scale, captain_ruling, composition_from=composition_from,
                      den_ruling=den_ruling, depth_from=depth_from, depth_exponent=depth_exponent,
                      depth_form=depth_form)
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, gata4_exponent=gata4_exponent,
                           captain_ruling=captain_ruling, composition_from=composition_from, size_scale=size_scale,
                           den_ruling=den_ruling, depth_from=depth_from, depth_exponent=depth_exponent,
                           depth_form=depth_form)
        elif captain_ruling == R76_ROW:  # row R76's check first, and only with --captain-ruling R76
            check_r76(form, gata4_exponent, size_scale, captain_ruling, composition_from=composition_from,
                      den_ruling=den_ruling, depth_from=depth_from, depth_exponent=depth_exponent,
                      depth_form=depth_form)
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, gata4_exponent=gata4_exponent,
                           captain_ruling=captain_ruling, composition_from=composition_from, size_scale=size_scale,
                           den_ruling=den_ruling, depth_from=depth_from, depth_exponent=depth_exponent,
                           depth_form=depth_form)
        elif is_depth_command(depth_from, depth_exponent, depth_form):
            check_depth_line(form, depth_from, depth_exponent, depth_form, captain_ruling,
                             composition_from=composition_from, gata4_exponent=gata4_exponent, size_scale=size_scale,
                             den_ruling=den_ruling)
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, captain_ruling=captain_ruling,
                           depth_from=depth_from, depth_exponent=depth_exponent, depth_form=depth_form)
        elif is_gonly_command(form, gata4_exponent, size_scale) or den_ruling is not None:
            check_c2_gonly(form, gata4_exponent, size_scale, captain_ruling, composition_from=composition_from,
                           den_ruling=den_ruling)
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, gata4_exponent=gata4_exponent,
                           captain_ruling=captain_ruling, composition_from=composition_from, size_scale=size_scale,
                           den_ruling=den_ruling)
        elif composition_from is not None or captain_ruling == R36_ROW:
            check_composition_from(form, composition_from, captain_ruling, gata4_exponent)
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, captain_ruling=captain_ruling,
                           composition_from=composition_from)
        elif form == C2G_FORM and gata4_exponent is not None:
            check_c2_gata4(form, gata4_exponent, captain_ruling)
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, gata4_exponent=gata4_exponent,
                           captain_ruling=captain_ruling)
        elif gata4_exponent is None and captain_ruling is None:
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads)
        else:
            check_gata4(form, gata4_exponent, captain_ruling)
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, gata4_exponent=gata4_exponent,
                           captain_ruling=captain_ruling)
    except DepthCheckFailed as e:
        bs.say(f"refused (the depth line's checks inside the build, before anything is written): {e}", "red")
        raise typer.Exit(EXIT_DEPTH_CHECK) from None
    except GonlyGuardFailed as e:
        bs.say(f"refused (the range guard or the 5% den check, before any cell is drawn): {e}", "red")
        raise typer.Exit(EXIT_GUARD_FAILED) from None
    except GuardFailed as e:
        bs.say(f"refused (the guard: the file's den-weighted total against C2's): {e}", "red")
        raise typer.Exit(EXIT_GUARD_FAILED) from None
    except FactorMismatch as e:
        bs.say(f"refused (the factor is not the one the line's figures describe): {e}", "red")
        raise typer.Exit(EXIT_FACTOR_MISMATCH) from None
    except Refused as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(2) from None
    chk = s["output"]["check"]
    cb = s["combined_build"]
    bs.say(f"form {cb['form']}: base {cb['base']['name']} (x = {cb['base']['x']}), Gata4 weight "
           f"{cb['combination']['gata4_weight']:.4f}, c = {cb['combination']['c']:.4f}")
    if cb.get("c2_gata4_line"):
        bs.say(f"row {G2_ROW}'s third T3 line: gata4_pooled weight g = {cb['c2_gata4_line']['gata4_weight_g']:g} in "
               f"C2's size-held product (a weight, not a power); the guard passed")
    if cb.get("c2_gonly_line"):
        gl = cb["c2_gonly_line"]
        lift = gl["den_ruling"]
        den = "passed" if lift is None else f"gates nothing for it (row {lift['row']})"
        bs.say(f"row {G2_ROW}'s options on C2's command: {gl['file']} (g = {gl['g']}, s = {gl['size_scale']:g}); the "
               f"range guard passed; the 5% den check {den}")
    if cb.get("composition_line"):
        bs.say(f"composition from {cb['composition_line']['composition_from']} (row {R36_ROW}): the "
               f"{R36_N_STATES} scenario states at {R36_RATIO:.4f} of their wild-type share before renormalising")
    if cb.get("depth_line"):
        dl = cb["depth_line"]
        how = (f"one ratio {dl['global_ratio']:.6g} for every label" if dl["form"] == "global" else
               f"the frozen per-label ratios to the power {dl['exponent']:g}")
        bs.say(f"the depth line (row {dl['captain_ruling']['row']}): {dl['file']}, {how}; the depth step received "
               "C2's counts row by row")
    if cb.get("r76_line"):
        rl = cb["r76_line"]
        lift = rl["den_ruling"]
        den = "passed" if lift is None else f"gates nothing for it (row {lift['row']})"
        bs.say(f"row {R76_ROW}'s refinement round on C2's command: {rl['file']} (g = {rl['g']:g}, s = "
               f"{rl['size_scale']:g}); the range guard passed; the 5% den check {den}")
    if cb.get("r79_line"):  # g is the string "inf" for T3_c2_s2_gonly_s05: r79_g_text prints it as it is
        rl = cb["r79_line"]
        lift = rl["den_ruling"]
        den = "passed" if lift is None else f"gates nothing for it (row {lift['row']})"
        bs.say(f"row {R79_ROW}'s refinement round 4 on C2's command: {rl['file']} (g = {r79_g_text(rl['g'])}, s = "
               f"{rl['size_scale']:g}); the range guard passed; the 5% den check {den}")
    for w in s["warnings"]:
        bs.say(f"  warning: {w}", "yellow")
    for e in chk["errors"]:
        bs.say(f"  check error: {e}", "red")
    bs.say(f"{'OK' if chk['ok'] else 'FAILED the format check'}: {s['output']['h5ad']}",
           "green" if chk["ok"] else "red")
    raise typer.Exit(0 if chk["ok"] else 1)


@app.command()
def g1(
    decision: Path = typer.Argument(..., help="decision file (experiments/t3_decomp/results/gata4/"
                                    "gata4_decisions.json)"),
    wt: Path = typer.Option(..., help="WT reference (WT E8.75), as for the build"),
    composition_from: str = typer.Option(..., help="allowed-s1 only"),
    threads: int = typer.Option(16, help="BLAS threads while preparing the reference (the frozen default)"),
):
    """G1, row R36's gate before any build: prints the five scenario states' shares over their wild-type shares
    before renormalising, one per line in the block's order (Python's shortest round-trip form), and nothing else;
    exit 0 when each is 0.3835131214774091 within a relative 1e-12, 1 otherwise, 2 refused, 3 not computed."""
    try:
        with redirect_stdout(io.StringIO()):  # anything the loading code prints is dropped: the ratios only
            ratios = run_g1(decision, wt, composition_from, threads)
    except Refused as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(EXIT_REFUSED) from None
    except Exception as e:  # the type only: a message or traceback could carry a figure
        bs.say(f"failed: {type(e).__name__}", "red")
        raise typer.Exit(EXIT_GATE_CRASHED) from None
    for r in ratios:
        print(repr(float(r)))
    raise typer.Exit(0 if ratios_pass(ratios) else EXIT_GATE_FAILED)


@app.command()
def g2(
    decision: Path = typer.Argument(..., help="decision file (experiments/t3_decomp/results/gata4/"
                                    "gata4_decisions.json)"),
    wt: Path = typer.Option(..., help="WT reference (WT E8.75), as for the build"),
    board: str = typer.Option(..., help="T3 board key in boards/"),
    run_dir: Path = typer.Option(None, help="run directory holding factors/"),
    seed: int = typer.Option(1_000_001, help="1000001, the list's seed; any other is refused"),
    n_cells: str = typer.Option("max", help="'max', the build's count: the board's max_cells, read in code, never "
                                "printed; any other is refused"),
    boards_dir: Path = typer.Option(None, help="boards directory (default: the repository's boards/)"),
    composition_from: str = typer.Option(..., help="allowed-s1 only"),
    threads: int = typer.Option(16, help="BLAS threads while preparing the reference (the frozen default)"),
):
    """G2, row R36's gate before upload: prints the six gate bars, one per line, and nothing else: C2's on N4's three
    gene sets (all genes, feasible down or up, the top fifth by |dp|), then cD7's on the same three; exit 0 when on
    each set cD7's bar is at most C2's plus 0.01, 1 otherwise, 2 refused, 3 not computed."""
    try:
        if seed != LIST_SEED or n_cells != LIST_N_CELLS:
            raise Refused(LIST_INPUTS_REFUSAL)
        if not isinstance(composition_from, str) or composition_from not in COMPOSITION_FROM:
            raise Refused(f"--composition-from {composition_from!r} is refused: only allowed-s1")
        b = bs.board_spec(board, boards_dir)
        with redirect_stdout(io.StringIO()):  # anything the loading code prints is dropped: the six bars only
            res = run_g2(decision, wt, run_dir, b, bs.resolve_n_cells(n_cells, b), seed, composition_from, threads)
    except Refused as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(EXIT_REFUSED) from None
    except Exception as e:  # the type only: a message or traceback could carry a count of the test build
        bs.say(f"failed: {type(e).__name__}", "red")
        raise typer.Exit(EXIT_GATE_CRASHED) from None
    for name in ("c2", "cd7"):
        for s in N4_SETS:
            print(repr(float(res[name][s])))
    raise typer.Exit(0 if res["passes"] else EXIT_GATE_FAILED)


if __name__ == "__main__":
    app()
