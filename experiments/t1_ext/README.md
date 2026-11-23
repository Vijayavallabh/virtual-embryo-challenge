# Task 1 outside heart atlases (stage 1)

Several public heart atlases that sample both sides of the protected window in one protocol give (1) a consensus
per-group direction, with a sign rule, in place of one atlas's estimate, and (2) a cardiomyocyte-share forecast for
E12.5. This folder is stage 1: the code that plans, pins, fetches, loads, labels and tables the atlases. Nothing here
changes a submission file. Task 1 only: E8.5 lies inside the T2 heart's interpolation window, so every manifest says
`tasks: ["T1"]`.

Written by the AI assistant (an engineer agent of the main session) on 2026-11-20, from a planner's design and a
challenger's review; third version 2026-11-21 (pinned samples as the first line, after the second verifier's three
blockers). No data file was fetched or read to write it, and no GEO record was opened; the tests use synthetic data
only.

## Atlases (`atlases.yaml`)

| Accession | Status | Ends used | Tissue | Groups formed |
|---|---|---|---|---|
| GSE193346 (X002) | in use | E9.5, E14.5 | chambers | CM, Endocardium, Endothelium, Epicardium (its ruled file, row N54) |
| GSE230531 | wave 1 | E8.5 embryo, E14.5 heart (E16.5 fallback) | 10x 3' | CM |
| GSE176063 + GSE153562 | wave 1 | E9.5 embryo, E14.5 fetal heart (only if its own record states the day) | Microwell-seq | CM |
| GSE244783 | wave 2 | E8.5, E16.5 heart | Multiome, RNA part | CM, Endocardium, Endothelium |
| GSE216970 | held | E9.5, E15.5 | 10x, MACS-enriched at E9.5 | Endocardium |

Wave 1 gives a CM-only consensus with k = 3; every other group keeps the ruled file's direction byte for byte.

## Two lines of defence

Two rounds of review showed that a text parser cannot be the only guard of the window rule: each round found more ways
a real record words a stage that the parser read as permitted. The same problem was closed in `experiments/t3_ext` by
a pinned list, and it is closed the same way here.

1. **The pins (first line).** A sample is used only when `atlases.yaml` pins it with the sha256 of its full verbatim
   record, after two readers have read that record line by line (the checklist below). Each series of the atlas is
   pinned by the sha256 of its series record. The plan, the fetch (on the live record) and the loader refuse a record
   that no longer hashes to its pin.
2. **The stage rule (second line).** `t1x_common.record_stage` reads the whole record and can only refuse: a pinned
   sample it reads as WINDOW or UNKNOWN, or at a stage other than its pin's, refuses the atlas. It never makes an
   unpinned sample usable.

A sample is **USED** only when it is pinned AND the rule does not refuse it. Nothing is fetched or loaded for an atlas
whose pinned samples do not cover both of its ends.

## Steps

1. `fetch.py plan ACC --dest D`: family SOFT metadata only. Every sample of every listed series is marked USED,
   WINDOW (never fetched; only the fields that gave a reading are kept), NOT_USED (permitted but not pinned, or not
   needed) or UNKNOWN. The plan writes `plan_table.tsv`, `records.txt` (every series record, and the full record of
   every sample the rule lets through, each under a header with its sha256), `records/<ACC>.txt` (each record's exact
   bytes, so `sha256sum records/GSMn.txt` is its hash) and `processing_text.txt`. Before any pin the plan refuses
   ("not pinned", "do not cover both ends"); candidates show `candidate end label` in `records.txt`.
2. Two readers read every record of the samples they mean to pin, every series record and the processing text, by the
   checklist below. `fetch.py pin-template ACC --dest D --gsm GSMa --gsm GSMb ...` prints, from the stored plan, the
   ready-to-paste `pinned:` block for the NAMED samples (no bulk selector) and every series, and the
   `processing_cleared` line for exactly those samples (it writes their processing text to
   `processing_text_pinned.txt`). It refuses a sample the rule does not let through, one that is not a candidate at a
   bracket end, a set that does not cover both ends, and a stored record that no longer hashes to the plan's value.
   `read_by` and `date` are placeholders (`READER_1`, `YYYY-MM-DD`) that the atlas check refuses until replaced.
3. The main session commits the pins and the clearance, syncs, and plans again: the plan must now pass.
4. `groups.py derive`: the extra competing classes for whole-embryo early ends, from the released E8.5 by a fixed
   rule. The main session commits the output, pins its sha256 and freezes `atlases.yaml`.
5. `fetch.py fetch ACC --dest D --plan-sha256 S`: the plan's sha256, its atlas entry (pins included) and its clearance
   are checked, and each USED sample against its pin; the metadata is downloaded and planned again, and must give the
   same processing text, USED samples, full records (each hashing to its pin), series records and files; then only the
   USED samples' files are fetched, from each sample's own GEO folder; never a RAW.tar. The manifest stores each
   sample's full record, its record_sha256 and the pins.
6. `load.py ACC --dest D`: window removal FIRST (the same stored record re-read by the rule before any file is opened);
   then every manifest sample must be pinned, its record hashing to its pin, its end and stage its pin's, the series
   records hashing to their pins, both ends covered; then the processing text rebuilt from the records and checked
   against the clearance; then the file gates (sha256, integer counts, every 10x barcode ending in -1, Multiome
   feature_type, dense-table rows), the whitelisted columns only, the fixed cell call, and a QC record whose window
   count is recomputed per cell from its sample's record. `read_reduced` re-checks every row of `cells.csv` the same
   way.
7. `tables.py`: per atlas, the labeller (`feng_relabel.relabel`, imported) on a derived config, the checks, the
   directions (`feng_dirs.directions_from_stats`, `direction_halves`, `eb_shrink`), the rescale to 5 days, the CM
   shares; then the consensus variants and the composition table. Outputs under `D/tables/<name>/` with a manifest.

`run_dgx.sh CHECKOUT MODE` runs these detached on the DGX (modes plan, derive, fetch, load, tables, rest).

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
   `_other_lines`.
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

## Processing clearance

The plan prints `D/meta/<ACC>/processing_text.txt` (also at the end of `plan_table.tsv`), and `pin-template` prints the
same text for the samples it pins (`processing_text_pinned.txt`). Canonical form `t1x-processing-v2`, pinned by
`tests/test_t1_ext_stage.py::test_processing_text_canonical_form`:

- line 1 `t1x-processing-v2 <ACC>`;
- per series in the atlas's order, `SERIES <GSE>`, then for `summary` and `overall_design` one line
  `<key>[<i>]: <value>` per value (i from 1), or `<key>: (none)`; then `series_supplementary_file[<i>]: <file name>`
  per supplementary file of the series record, or `series_supplementary_file: (none)` (the verifier's lesser point 4;
  the version was raised from v1 so every older clearance lapses);
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
  with fewer than two members in a variant is the ruled vector unchanged.
- Composition: the forecast is the within-atlas change form, logit(x) = logit(s_base) + days x median per-day slope
  (days 2 at the test, 1 on practice). The level form is a sensitivity: logit(x) = logit(s_base) + f (L - logit(s_base))
  (f 0.5 at the test, 0.2 on practice), L the median E14.5 level. Offered only with a change atlas, at least two levels
  within a 1.0-logit half-range, and both steps of the same sign. r_mesenchyme counts as heart only from E14.

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
  that each sample reads UNKNOWN; any exception is the main session's rule, not code.
- Dense-table row names are matched whole, or token by token for names with separators, not as substrings: the
  labeller's pattern as a substring would refuse every real table (Epcam holds "pca"). Snn (stannin) is the one gene
  exempted from a whole-name match.
- The change form uses one label set at both ends (never r_mesenchyme); the level form adds r_mesenchyme at a late end
  from E14 (the second review's point 3).
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
  (a convention sentence with a day and no recognised convention), a protocol shared by every sample that lists every
  stage. Its false-refusal rate on real GEO records is unknown. A refusal shows its reason in the plan table; no
  clearance or pin can override it.
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
