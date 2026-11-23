# Task 1 outside heart atlases (stage 1), strict form

Candidate MX2s of the outside-data round of 2026-11-23: a within-study sign-agreement consensus of three heart atlases
for the CM, Endocardium and Endothelium directions of the T1 direction step, at the ruled step's size. The captain's
word (ruling row R27 in `experiments/harness/TRIES.md`) reopened this parked line in its strict form, practice read
only, with the stage-parser exception below; register rows X019 (GSE244783) and X020 (GSE230531) and X002's MX2s
amendment are in `EXTERNAL_DATA.md`. This folder is stage 1: the code that plans, pins, screens, fetches, loads,
labels and tables the atlases (the specification's METHOD steps 1 to 3 and its SCREENS). The consensus tables of the
specification (its steps 4 to 6) and the practice read (steps 7 and 8) are other engineers' work after this one.
Nothing here changes a submission file. Task 1 only: E8.5 lies inside the T2 heart's interpolation window, so every
manifest says `tasks: ["T1"]` and no code outside this folder names it (`tests/test_t1_ext_mx2s.py`).

Written by the AI assistant (an engineer agent of the main session) on 2026-11-20 from a planner's design and a
challenger's review; third version 2026-11-21 (pinned samples as the first line); strict form 2026-11-24 from the cold
challenger's final specification of MX2s. No data file was fetched or read to write it, and no GEO record was opened;
the tests use synthetic records and synthetic files only.

**Who may read records.** Only two independent readers (two screening agents, neither the engineer who wrote this
code) read the GEO records, the processing text and the screens' output; the main session alone commits their pins,
the processing clearance and the screens record.

## Atlases (`atlases.yaml`)

| Accession | Status | Ends used (samples named) | Tissue | Groups formed | Cell call (UMIs, genes, mito share) |
|---|---|---|---|---|---|
| GSE193346 (X002) | in use | E9.5, E14.5 | chambers | CM, Endocardium, Endothelium, Epicardium (its ruled file, row N54) | as N54 |
| GSE244783 (X019) | use | E8.5 (GSM7829699), E16.5 (GSM7829703) | heart, Multiome RNA part | CM, Endocardium, Endothelium | 500, 300, 0.05 |
| GSE230531 (X020) | use | E8.5 (GSM7226268, GSM7226269), E16.5 (GSM7226274, GSM7226276) | E8.5 embryo, E16.5 heart, 10x 3' v2 unfiltered | CM | 1000, 500, 0.20 |
| GSE176063 + GSE153562 | held | | Microwell-seq | | |
| GSE216970 | held | | 10x, MACS-enriched at E9.5 | | |

Never fetched, by name (`never_fetch`; the fetch refuses them, and a pin naming one is refused): GSE244783's
GSM7829701 (E10.5), its ATAC samples GSM7829700, GSM7829702 and GSM7829704, and GSE244783_RAW.tar; GSE230531's
GSM7226279 and GSM7226280 (E10.5), GSM7226270 and GSM7226271 (E12.5), GSM7226272 and GSM7226273 (E14.5, dropped
under the strict form) and GSE230531_RAW.tar. GSE230531's P3 samples are not listed by accession (none was read here):
they are not named, so they can never be pinned, and only pinned samples are fetched. Every series archive
(`<GSE>_RAW.tar`) of an atlas is refused by name whether listed or not, and every data file must lie in its own
sample's GEO folder.

E8.5 and E16.5 lie outside (E9.5, E13.5] under any plug-day convention in [-1.0, +0.5] day (E8.5 reads at most 9.0,
E16.5 at least 15.5). The dissection differs between GSE230531's ends (whole embryo, then heart): recorded, and guarded
by the competing classes and the purity gate. No composition is computed: the composition block (level and change
forms) is deleted from `tables.py` and `atlases.yaml`, and no share is written.

## Two lines of defence

Two rounds of review showed that a text parser cannot be the only guard of the window rule: each round found more ways
a real record words a stage that the parser read as permitted. The same problem was closed in `experiments/t3_ext` by
a pinned list, and it is closed the same way here.

1. **The pins (first line).** A sample is used only when `atlases.yaml` names it and pins it with the sha256 of its
   full verbatim record, after two readers have read that record line by line (the checklist below). Each series of
   the atlas is pinned by the sha256 of its series record. The plan, the fetch (on the live record), the loader and
   the tables' read of a reduction refuse a record that no longer hashes to its pin.
2. **The stage rule (second line).** `t1x_common.record_stage` reads the whole record and can only refuse: a pinned
   sample it reads as WINDOW or UNKNOWN, or at a stage other than its pin's, refuses the atlas. It never makes an
   unpinned sample usable.

A sample is **USED** only when it is pinned AND the rule does not refuse it. Nothing is fetched or loaded for an atlas
whose pinned samples do not cover every named sample at both of its ends.

## The stage-parser exception (ruling row R27)

In the ruling's words: "a protocol or processing field that is byte-identical across every sample record of a series
is read as series text: it gives no stage reading, and it still refuses on pooling words; the sample's stage must then
come from its own title and characteristics, which must agree". As coded (`t1x_common.series_wide_fields`,
`record_stage`, `stage_reading`), and no wider:

- A field is one key of the record with all its values (GEO's repeated `!Sample_<key>` lines, in order); only the
  protocol and processing fields (`*protocol*`, `data_processing`) can be series text. It is series text when every
  sample record of the series (window and other-assay samples included) holds that key with the same values in the
  same order; a field missing from one sample or differing by one character is not. A series with one sample has
  none. The plan computes it from the whole family SOFT; the fetch's live re-plan must give the same; it is stored in
  the manifest and printed in the processing text (`series_text <key>[<i>]: <value>`, or `series_text: (none)`), so
  the readers' clearance covers it.
- Series text gives no stage reading, no leftover and no companion integer. It still feeds the plug-day convention
  (as the series' own summary and design do). A pooling word in it (pool, pools, pooled, pooling, merge, merges,
  merged, merging) makes the sample UNKNOWN.
- When a sample has a field read as series text, its stage must come from its own title AND its own characteristics:
  each must give a reading and the two must be the same nominal interval, else UNKNOWN. Every other field (source
  name, description, file names, the sample's own non-series protocol text) is still read and can still refuse.
- The exception can only turn a refusal into a reading from the title and characteristics; a sample whose own title
  or characteristics read in the window or unknown stays WINDOW or UNKNOWN. The rule stays the second line: it never
  makes an unpinned sample usable.

## Screens (`screens.py`, `screens.yaml`)

The specification's SCREENS, run on the fetched barcode and feature files before any matrix is fetched or opened; each
returns yes or no with counts and never a cell. The checks per atlas are its `screens` list; their constants are
`screen_rules` in `atlases.yaml`.

- (a) `suffix`: every barcode of every sample ends in -1 (no aggr).
- (b) `own_peaks` (GSE244783): the two samples' features files have different sha256 and different sets of Peaks rows
  (cellranger-arc count per library, not aggr), and each has Peaks and Gene Expression rows; the Peaks ids are hashed
  in memory only to compare the sets. Only Gene Expression rows are ever read as expression (by the loader).
- (c) `no_hashing`: one column per barcode line (no per-cell sample, stage or tag column); features files of 2 or 3
  columns (Multiome: 3 or 6); no feature type outside the assay's (Gene Expression; Multiome also Peaks); no hashtag,
  CMO or antibody-tag-like feature name.
- (d) `whitelist_sized` (GSE230531): every barcode file has a 10x whitelist's number of lines (737,280 or 6,794,880;
  general knowledge, not checked against a source): unfiltered matrices, and the team's own cell call applies.

A failing check marks the member dropped in `screens.yaml` (committed by the main session); the fetch of its matrices,
its load and its tables then refuse.

## Steps

In order; the main session commits between steps, and every output lies under `D` (`~/v/vec-ext/t1` on the DGX),
outside every checkout. `run_dgx.sh CHECKOUT MODE` runs each step detached on the DGX (modes plan, derive, heads,
screens, fetch, load, tables, rest).

1. `fetch.py plan ACC --dest D`: family SOFT metadata only. Every sample of every listed series is marked USED,
   WINDOW (never fetched; only the fields that gave a reading are kept), NOT_USED (permitted but not pinned, not
   named, on the never-fetch list, or not needed) or UNKNOWN. The plan computes each series' series-wide fields (the
   exception below) and writes `plan_table.tsv`, `records.txt` (every series record with its series-wide fields
   listed, and the full record of every sample the rule lets through, each under a header with its sha256),
   `records/<ACC>.txt` (each record's exact bytes) and `processing_text.txt`. Before any pin the plan refuses ("not
   pinned", "do not cover both ends", "named samples not pinned").
2. Two independent readers read every record of the named samples, every series record and the processing text, by
   the checklist below. `fetch.py pin-template ACC --dest D --gsm GSMa --gsm GSMb ...` prints, from the stored plan,
   the ready-to-paste `pinned:` block for exactly the samples the atlas names (no bulk selector; any other set is
   refused) and every series, and the `processing_cleared` line for those samples. `read_by` and `date` are
   placeholders (`READER_1`, `YYYY-MM-DD`) that the atlas check refuses until replaced.
3. The main session commits the pins and the clearance, syncs, and plans again: the plan must now pass.
4. `groups.py derive`: the extra competing classes for GSE230531's whole-embryo early end, from the released E8.5 by
   a fixed rule. The main session commits the output, pins its sha256 and freezes `atlases.yaml`.
5. `fetch.py heads ACC --dest D --plan-sha256 S`: every check of the fetch (the plan's sha256, its atlas entry with
   the pins, the clearance, a live re-plan that must give the same processing text, USED samples, records, files and
   series-wide fields), then each USED sample's barcode and feature files only, from its own GEO folder; never a
   matrix. Writes `D/raw/<ACC>/HEADS.json`.
6. `screens.py ACC --dest D`: the SCREENS below, on those files, before any matrix is fetched or opened. It writes
   `D/screens/<ACC>/screen.json` and prints yes or no per check with counts (never a cell, barcode or gene) and the
   line to paste into `screens.yaml`: `kept` when every check said yes, `dropped` otherwise.
7. The main session commits `screens.yaml`. A dropped member is dropped before any table is made; no exception is
   granted on the spot.
8. `fetch.py fetch ACC --dest D --plan-sha256 S`: the same checks as heads, then only for a member `screens.yaml`
   records kept, with a screen file that hashes to the record, says every check passed and ran on the heads on disk:
   the heads are re-hashed (not fetched again) and only the matrices are fetched. Writes `D/raw/<ACC>/MANIFEST.json`
   (with the series-wide fields and the screen).
9. `load.py ACC --dest D`: window removal FIRST (each stored record re-read by the rule, with its series-wide fields,
   before any file is opened); then every manifest sample pinned, its record hashing to its pin, both ends covered;
   the processing text rebuilt from the records (series-wide fields included) and checked against the clearance; the
   screens record kept for this member; then the file gates (sha256, integer counts, every 10x barcode ending in -1,
   Multiome feature_type with only Gene Expression rows read), the whitelisted columns only, the fixed cell call and
   the QC record.
10. `tables.py`: each reduction is read by `load.read_reduced` under the atlas entry, pins, clearance and screens
    record in force (the parked branch's known blocker, closed below); the labeller on a derived config, the checks,
    the directions, the rescale to 5 days, the consensus variants. No composition, no share.

## Pinned samples

Per atlas in `atlases.yaml` (written from `pin-template` output, never by hand):

```yaml
    pinned:
      samples:
        - accession: GSMn
          end: early            # or late
          stage: E8.5           # one of the end's labels; all pins of one end name one stage
          record_sha256: "<sha256 of the sample's canonical record text>"
          read_by: [<reader 1>, <reader 2>]   # two different readers; placeholders are refused
          date: "YYYY-MM-DD"
      series:
        - accession: GSEn       # every series of the atlas, once
          record_sha256: "<sha256 of the series' canonical record text>"
          read_by: [<reader 1>, <reader 2>]
          date: "YYYY-MM-DD"
    processing_cleared: {sha256: "<hex>", read_by: [...], date: "YYYY-MM-DD"}
```

The pinned block is part of the atlas entry's sha256 (a plan made under other pins is not fetched);
`processing_cleared` is not (it is checked on its own against the text). An atlas in use or held holds no pins.

## Canonical record text

Form `t1x-record-v1`, pinned by `tests/test_t1_ext_stage.py::test_record_text_canonical_form`:

- line 1 `t1x-record-v1 SAMPLE <GSM>` (or `t1x-record-v1 SERIES <GSE>`);
- then EVERY field of the record as `geo.parse_soft` gives it (the `Sample_` or `Series_` prefix removed, values
  stripped of surrounding whitespace and otherwise verbatim), keys in the order they first appear in the SOFT block, a
  key's values together in file order, one line `<key>[<i>]: <value>` per value (i from 1), or `<key>: (none)`.
  Nothing is left out: protocol and processing fields, relations, every supplementary file URL, contacts, dates and
  the lines the parser kept under `_other_lines` (lines of the block that are not `!Sample_` fields, such as a data
  table);
- under any value holding a character outside printable ASCII (32-126; a tab counts), one more line
  `<key>[<i>] escaped: <the value with Python's unicode_escape>`;
- lines joined by `\n`, one `\n` at the end, UTF-8. `record_sha256` is the sha256 of those bytes.

Stored records (plan.json, the fetch manifest) are lists of `[key, [values]]` pairs, because the JSON writer sorts
object keys; the order is part of the text.

## Readers' checklist (before any pin)

Each of two readers reads, line by line, the full record of every sample to be pinned (`records/GSMn.txt`), every
series record (`records/GSEn.txt`) and `processing_text_pinned.txt`, and confirms for each sample:

1. **One stage, permitted, stated by the record.** The record states exactly one stage, it is the pin's stage, it is
   outside (E9.5, E13.5] on the plug = E0.5 scale, and it is stated by the record itself (not inferred from an
   accession range, a file order or the series summary).
2. **No other stage anywhere.** No other stage, range, age, somite count or Theiler stage, no pooling across stages and
   no littermate or "collected at" wording that names another day, in ANY field: title, source name, every
   characteristics line, description, the protocol and processing text, relations, labels, file names and
   `_other_lines`. The one exception is ruling row R27's: a field the plan lists as series-wide (byte-identical in
   every sample record of the series, `records.txt` and `series_text` lines of the processing text) may name the
   series' other stages, provided it names no pooling across stages; the sample's own title and characteristics must
   then both state its stage and agree.
3. **The plug convention.** The series or the sample states the plug-day convention; for E9.5 and E14.5 ends it must
   be plug = E0.5 (E14.5 also E0). No second sentence about mating, conception or fertilisation gives another day.
4. **One library of one sample.** The record is one library from one sample (one animal or one stated pool of one
   stage), not a merge of libraries or a re-split of an aggregated matrix; its files are its own, named after it.
5. **Per-sample processing.** The matrix was made from that sample's own library (Cell Ranger count or the
   Microwell-seq pipeline per sample): no aggr, merge, integration, batch correction, joint normalisation, joint or
   all-stage quality filter, ambient or doublet removal learned across samples, no downsampling. The series'
   supplementary file names (in the processing text) give no sign of joint processing that the per-sample files came
   from.
6. **No enrichment that differs between the ends.** Sorting, MACS, FACS, dissociation or enrichment is the same at
   both ends.
7. **Hidden characters.** Every `escaped:` line matches its plain line; no look-alike letters.

A reader who doubts any point does not pin. Both readers must confirm every point for every pinned sample and series;
the pin names them both.

## Plug-day conventions

| Stated | Shift onto the plug = E0.5 scale | Why |
|---|---|---|
| E0.5 (or 0.5 dpc) | 0 | the scale the window is written in |
| E0 (or day 0) | [0, +0.5] | a label may count from the plug day's midnight or its noon |
| E1 (day 1, first day) | [-1.0, -0.5] | one day later than E0's count; E14.5 reads [13.5, 14.0], in the window |
| none | [-1.0, +0.5] | any of the three; E9.5 and E14.5 then touch the window, E8.5 and E16.5 do not |
| another such sentence with a day, or two conventions | UNKNOWN | fail closed |

A convention sentence is one that names a plug, mating, coitus (not "post coitus", a dpc unit), copulation,
conception or fertilisation (the verifier's lesser point 2).

## The stage rule (second line)

`t1x_common.record_stage` reads every field of the sample's record except `supplementary_file*` (read as file names),
`contact_*`, `platform_id`, `submission_date`, `last_update_date`, `status` (it holds the release date),
`data_row_count` and `_other_lines` (the verifier's lesser point 3, with `status` added). Text is NFKC-normalised
first (a full-width E reads as E). Descriptive fields are read value by value. The sample's own protocol and
processing fields (`*protocol*`, `data_processing`) are read sentence by sentence, with the convention's own tokens
(E0.5, day 1, 0.5 dpc) blanked in convention sentences; a reading there counts as a reading anywhere else (blocker 3).
Prefixes read: E, ED, GD, e.d., embryonic day (also glued), dpc, pc; ranges and lists closed by a dpc unit ("13.5-14.5
dpc", "12.5 and 14.5 days post-coitum") read as their span (blocker 1). Leftovers that give UNKNOWN: E, TS, dpc, ED,
GD followed by digits; a number with a d suffix; embryonic or stage before a number; somite words; the words day,
week, month and gestation (midgestation too) in descriptive fields; a half-day number no pattern read; and, in a value
or sentence that gave a reading, an integer from 5 to 21 (blocker 2).

## The known blocker of the parked branch, closed

The parked branch's tables step read a reduction left on disk (`load.read_reduced`) with its file hashes and window
checks but without a pin check, so a reduction made under earlier pins, another atlas entry, a lapsed clearance or a
later drop could be tabled. `read_reduced(acc, dest, cfg)` now refuses unless the reduction was made under the atlas
entry in force (pins included), every sample it holds is pinned now with its record hashing to its pin and its stage
its pin's, the pins cover every named sample at both ends, the clearance holds for its processing text, and a screened
member is still recorded kept (`tests/test_t1_ext_mx2s.py::test_the_tables_step_reads_no_reduction_made_under_other_pins`).
The launcher no longer deletes an old reduction before a load (it deletes nothing); the check above is the guard.

## Processing clearance

The plan prints `D/meta/<ACC>/processing_text.txt` (also at the end of `plan_table.tsv`), and `pin-template` prints the
same text for the samples it pins (`processing_text_pinned.txt`). Canonical form `t1x-processing-v2`, pinned by
`tests/test_t1_ext_stage.py::test_processing_text_canonical_form`:

- line 1 `t1x-processing-v2 <ACC>`;
- per series in the atlas's order, `SERIES <GSE>`, then for `summary` and `overall_design` one line
  `<key>[<i>]: <value>` per value (i from 1), or `<key>: (none)`; then `series_supplementary_file[<i>]: <file name>`
  per supplementary file of the series record, or `series_supplementary_file: (none)` (the verifier's lesser point 4;
  the version was raised from v1 so every older clearance lapses);
- (since the strict form, when the plan, fetch or load pass the series-wide fields) after each series' file names,
  `series_text <key>[<i>]: <value>` per series-wide field, keys sorted as text, or `series_text: (none)`;
- per USED sample, sorted by accession as text, `SAMPLE <GSM> (<GSE>)`, then the same for `data_processing`,
  `extract_protocol_ch1`, `description`, `growth_protocol_ch1`, `treatment_protocol_ch1`;
- under any value holding a character outside printable ASCII (a tab counts), `<key>[<i>] escaped: <value with
  Python's unicode_escape>`;
- lines joined by `\n`, one `\n` at the end, UTF-8. The sha256 is of those bytes (`sha256sum processing_text.txt`).

**Kept as its own clearance, not folded into the pins.** Every line of the processing text also lies inside a pinned
record (the series records and the USED samples' records), so folding would lose no line. It is kept because it puts
the processing of all USED samples and the series' file names side by side, which no single record shows, and because
the brief asked that it stay as it was unless folding was argued. `pin-template` prints its hash with the pins, so the
pins and the clearance are written and committed together.

## Fixed rules

- Consensus: members with a finite non-zero value count (an exact zero does not); the majority sign must be held by
  at least two members and more than half; c = (1 / k) times the sum of the agreeing values, else 0. Equal weights.
- Variants: `deciding` (each group rescaled to the ruled direction's L2 norm), `raw`, `noX002`, `noE85start`. A group
  with fewer than two members in a variant is the ruled vector unchanged. (The specification's own rescale, over the
  non-held genes where the ruled direction is finite, is its step 6: other engineers' work after this one.)
- Composition: deleted in the strict form. No share, level or change form is computed or written, and `groups.py` has
  no share check.

## Where this differs from the plan, the reviews and the verifier

- The challenger's changes are all taken except where noted here.
- A conflicting reading makes the sample UNKNOWN even when one reading lies in the window; a window mark never hides a
  conflict. An UNKNOWN sample refuses the atlas only if the atlas would use it.
- The verifier's half-day leftover is applied as written in descriptive fields. In protocol and processing prose it
  counts only half-day numbers from 9.5 to 14.5 (any number that could touch the window under some convention): 0.5,
  1.5 and 2.5 are everywhere in protocols (0.5 mg/ml, 1.5 ml tubes, clustering resolution 0.5). Both forms skip a
  number followed by a unit (%, ml, mg, min, h, x, ...) or preceded by v (a version).
- The verifier's companion-integer rule is applied as written (5 to 21, in a value that gave a reading), with a
  protocol sentence counting as a value, and numbers glued to letters (rep12, 10x) or followed by a unit skipped.
- In protocol and processing prose the rule reads no weeks, months, adult, neonatal or postnatal forms: none can meet
  the window, and in protocols they describe the dams ("8-12 weeks old females") or reagents ("primer P3"). Embryonic
  forms, Theiler stages, somites and days are read there as anywhere else.
- "stage" (or "embryonic" without "day") before a number is a leftover (UNKNOWN), not read as an E-day as the
  verifier offered: a bare "stage 12" may be a Theiler, Carnegie or day count. "day" before a number, glued or not
  (Heartday12.5), is read as days, either reading, as the verifier asked.
- Blocker 3 is closed by the brief's wording rather than the verifier's narrower fix: a stage read in the sample's own
  protocol text is a reading like any other (it can give WINDOW, UNKNOWN by conflict, or a leftover), not only an
  UNKNOWN after a USED decision. The default for shared protocol text listing every stage of a series is therefore
  that each sample reads UNKNOWN; ruling row R27 then set the one exception, coded above.
- Dense-table row names are matched whole, or token by token for names with separators, not as substrings: the
  labeller's pattern as a substring would refuse every real table (Epcam holds "pca"). Snn (stannin) is the one gene
  exempted from a whole-name match.
- The deciding variant uses the second review's norm rule: the ruled norm over genes where the consensus is finite and
  non-zero. Its record holds the factor, the share of genes at zero and the largest |deciding| / largest member.
- The extra competing classes are added only to whole-embryo atlases' configs, not to heart-dissected ones.
- The purity gate scores released cells with the atlas's derived config, whose controls come from the released file's
  own gene means; the atlas's cells use controls from the atlas. The midpoint comparison is across two control draws.

## Not checked

- No GSM or GSE record was opened. The stage field names, the plug-day statements, the per-sample chemistry and
  enrichment text, the file names and whether GSE176063's fetal heart record states a day are all unknown until the
  plan runs. The pins guard against what the rule misses only as well as the two readers read.
- False refusals. The rule fails closed on ordinary text that looks like a stage token, and it now reads protocol text
  and every field, so more real records may refuse: TIE2 in capitals, E2f1, a "D1" lane name, "the same day", "3D",
  "stage 2", a count such as "12 hearts" beside a stage, "timed-mated females" in a sentence that also names a stage
  (a convention sentence with a day and no recognised convention). A protocol shared by every sample that lists every
  stage is now read as series text (R27), but a series-wide field that also says "pooled" (as in "libraries were
  pooled and sequenced") refuses every sample of the series, and a protocol that differs by one character between
  samples is not series text. Its false-refusal rate on real GEO records is unknown. A refusal shows its reason in the
  plan table; no clearance or pin can override it.
- The screens' constants (the whitelist sizes, the hashtag-like names) are general knowledge, not checked against a
  source; the screens ran only on synthetic files.
- The verifier's lesser point 5: "genotype: Nkx2-5 E1 enhancer KO" reads E1 as a stage and conflicts (UNKNOWN), and a
  somite key such as "somite number: ..." is always a leftover. Both fail closed (usability costs only).
- The verifier's lesser point 1: the rule is not applied to the atlas in use. Without a stated E0.5 the unstated
  widening makes any E9.5 or E14.5 end WINDOW, but X002 (GSE193346, E9.5 and E14.5) is not planned here, so it is
  never checked the same way. Whether X002's own records state its plug convention is for the main session. The
  verifier's replacement text for the register amendment: "X002's ends (E9.5, E14.5) were not re-read under the T1-ext
  stage rule of bcd5832; under that rule they are permitted only if the plug convention is stated as E0.5 (E0 also for
  E14.5)."
- The Theiler and somite tables are general knowledge, widened, not checked against a source. Whether the Microwell-seq
  gene symbols hold every set gene (min_genes_found) is unknown.
- feng_member and released_inputs were tested only on stubs; groups and tables were not run end to end on real files.
