# Task 1 outside heart atlases, row T72 (candidate MX2s)

Candidate MX2s of the outside-data round of 2026-11-23: a within-study sign-agreement consensus of three heart atlases
for the CM, Endocardium and Endothelium directions of the T1 direction step, at the ruled step's size; practice read
only, under the captain's word of ruling row R27 (`experiments/harness/TRIES.md`; its third cell words the stage-parser
exception). Register rows X002 (its MX2s amendment), X019 (GSE244783) and X020 (GSE230531) are in `EXTERNAL_DATA.md`.
This folder is the row's DATA SIDE: it plans, pins, screens, fetches, loads, labels and writes the MEMBER TABLES of the
two new atlases. The consensus, the rescale, the two controls (`mx_mask`, `mx_perm`), the pre-gate, the runs and the
read are the practice-read script's (`experiments/t1_m1/mx2s_read.py`, launcher `experiments/t1_m1/mx2s.sh`), which
reads the tables in section B's format below. Nothing here changes a submission file. Task 1 only: E8.5 lies inside the
T2 heart's interpolation window, so every manifest says `tasks: ["T1"]` and no shared-core, T2, T3 or builder path
names this folder or its data (`tests/test_t1_ext_imports.py`).

Written by the AI assistant (engineer agents of the main session): the parked branch 2026-11-20 and 2026-11-21, the
strict form's port 2026-11-23, the row's data side 2026-11-24 from row T72's cells as revised after the cold attack
and the rule writer's engineer notes. No data file was fetched or read to write it, and no GEO record was opened; the
tests use synthetic records and synthetic files only.

**Who may read records.** Only two independent readers (two screening agents, neither the main session nor the
engineer who wrote this code) read the GEO records, the per-field hashes, the processing text and the screens' facts;
the main session alone commits their pins, the processing clearance and the screen blocks.

## Atlases (`atlases.yaml`)

| Accession | Register | Status | Ends used (samples named) | Tissue | Groups tabled | Cell call (UMIs, genes, mito share) |
|---|---|---|---|---|---|---|
| GSE193346 | X002 | in use | E9.5, E14.5 | chambers | (its ruled file, row N54; read by the practice script) | as N54 |
| GSE244783 | X019 | use | E8.5 (GSM7829699), E16.5 (GSM7829703) | heart, Multiome RNA part | CM, Endocardium, Endothelium | 500, 300, 0.05 |
| GSE230531 | X020 | use | E8.5 (GSM7226268, GSM7226269), E16.5 (GSM7226274, GSM7226276) | E8.5 embryo, E16.5 heart, 10x 3' v2 unfiltered | CM | 1000, 500, 0.20 |
| GSE176063 + GSE153562 | | held | | Microwell-seq | | |
| GSE216970 | | held | | 10x, MACS-enriched at E9.5 | | |

A member the readers' screens drop gets `status: dropped` with `screen: {verdict: dropped, failed: [points]}`: it is
never planned, fetched, loaded or tabled again, and no exception is granted on the spot.

NEVER FETCHED, by name, from a literal list in code (`t1x_common.NEVER_FETCHED`, row T72's cells 3 and 8; the yaml's
`never_fetch` lists must name exactly the same): GSM7829701 (E10.5), the ATAC samples GSM7829700, GSM7829702 and
GSM7829704, and GSE244783_RAW.tar; GSM7226279 and GSM7226280 (E10.5), GSM7226270 and GSM7226271 (E12.5), GSM7226272
and GSM7226273 (E14.5, dropped under the strict form) and GSE230531_RAW.tar. GSE230531's P3 samples are not listed by
accession (none was read here): they are not named, so they can never be pinned, and only pinned samples are fetched.
Every `_RAW.tar` is refused by name whether listed or not, and every data file must lie in its own sample's folder.

E8.5 and E16.5 lie outside (E9.5, E13.5] under any plug-day convention in [-1.0, +0.5] day (E8.5 reads at most 9.0,
E16.5 at least 15.5). GSE230531's ends differ in dissection (whole embryo, then heart): recorded (check (e)), guarded at
the E8.5 end by the derived competing classes and at the E16.5 end by the skm class and the purity clause. No
composition, level or share is computed or written.

## Two lines of defence

Two rounds of review showed that a text parser cannot be the only guard of the window rule. The same problem was closed
in `experiments/t3_ext` by a pinned list, and it is closed the same way here.

1. **The pins (first line).** A sample is used only when `atlases.yaml` names it and pins it with the sha256 of its
   full verbatim record, after two readers have read that record line by line. Each series of the atlas is pinned by
   the sha256 of its series record. The plan, the fetch (on the live record), the loader and the tables' read of a
   reduction refuse a record that no longer hashes to its pin.
2. **The stage rule (second line).** `t1x_common.record_stage` reads the whole record and can only refuse: a pinned
   sample it reads as WINDOW or UNKNOWN, or at a stage other than its pin's, refuses the atlas. It never makes an
   unpinned sample usable.

## The stage-parser exception (ruling row R27; row T72, cell 3)

"A protocol or processing field that is byte-identical across every sample record of the series is read as series
text: it gives no stage reading and still refuses on pooling words; the sample's stage must then be stated both by its
own title and by its own characteristics, and the two must be the same stage; a sample whose title or characteristics
state no stage is refused; every other part of the stage rule still refuses." As coded (`t1x_common.series_wide_fields`,
`record_stage`, `stage_reading`, `pooling_in_series_text`), and no wider:

- A field is one key of the record with all its values; only the protocol and processing fields (`*protocol*`,
  `data_processing`) can be series text. It is series text when every sample record of the series in the family SOFT
  (window and other-assay samples included) holds that key with the same values in the same order. A series with one
  sample has none. The plan computes it in code; the fetch's live re-plan must give the same; it is printed in the
  processing text (`series_text <key>[<i>]: <value>`), so the readers' clearance covers it.
- THE PER-FIELD HASHES. The plan prints, for every sample record of the series (window samples included), the sha256
  of each protocol and processing field and nothing else of a window sample: one line `field_sha256 <GSE> <GSM>
  <field> <sha256 or (absent)>` per field on stdout and in `meta/<ACC>/protocol_field_sha256.tsv`. A field's sha256
  is that of its canonical lines (`<key>[<i>]: <value>`, with the `escaped:` line under a non-ASCII value), joined by
  newlines with one at the end, as in the canonical record text. The readers compare these hashes and open no record
  of a never-fetched sample. The exception applies only where the code's test and both readers' comparison agree.
- Series text gives no stage reading, no leftover and no companion integer; it still feeds the plug-day convention.
- THE POOLING CLAUSE (the reviews of 290da3f and 75356ec). Each value of series text is read as written and with its
  Unicode format characters (category Cf) removed; a match in either refuses (UNKNOWN). (i) A form of pool or merge
  anywhere refuses. (ii) A staged pooling word (`POOLING_WORDS_STAGED`: combine, combination(s), in combination, mix,
  mixes, mixed, mixing, mixture(s), multiplex..., superload..., aggregat..., also with a glued or hyphenated prefix
  pre, inter, ad, co or re: premixed, intermixed, admixed, co-mixed, recombined) in ANY sentence of a value refuses
  when the value AS A WHOLE names stages: any of its sentences gives a stage reading or a leftover in protocol mode,
  or the value names a stage, an age, a time point (time point, time-point, timepoint) or a time course (time
  course, time-course, timecourse). The U+2010 and U+2011 hyphens count as hyphens.
- Benign prose that stays PERMITTED: "cDNA was combined with RT Master Mix.", "Cells were mixed with Trypan blue and
  counted.", "Libraries were prepared with the Single Cell Multiome kit.", and the same reagent prose in a value of
  several sentences none of which names a stage. Benign prose the whole-value rule now REFUSES (fails closed; tested):
  "Hearts were dissected at E16.5. cDNA was combined with RT Master Mix.", "cDNA was combined with RT Master Mix and
  incubated for 2 days." and "Samples were mixed by pipetting at each stage of the protocol.".
- When a sample has a field read as series text, its stage must be stated by its own title AND its own
  characteristics: each must give a reading and the two must be the same nominal interval, else UNKNOWN. Every other
  field is still read and can still refuse; a sample whose own title or characteristics read in the window stays
  WINDOW or UNKNOWN.

## The screens (before the row's commit, before any matrix)

`fetch.py fetch-meta ACC --dest ~/v/scratch/t72_screen --plan-sha256 S` requests only the barcodes and features files
of the six named samples (`t1x_common.SCREEN_SAMPLES`: GSM7829699, GSM7829703, GSM7226268, GSM7226269, GSM7226274,
GSM7226276), each from the sample's own folder; it refuses by name every never-fetched accession, every matrix and
every RAW.tar before any request; it writes only under `~/v/scratch/t72_screen` (exactly that folder: never the row's
data folders), and prints each file's name, bytes and sha256. It then runs `screens.py`, which prints, counts and
yes/no facts only (never a cell, a barcode or a gene):

- (a) `suffix`: per barcodes file its line count (as `zcat | wc -l`) and the count of barcodes not ending in -1;
- (b) `own_peaks` (GSE244783): whether the two features files have different sha256; whether at least one Peaks
  interval of each file overlaps an interval of the other without being identical to it (per-library peak calls;
  subsets of one merged peak set never partly overlap; the relation is symmetric, so one answer serves both files);
  whether the two files' Gene Expression rows are the same features in the same order;
- (c) `no_hashing`: the feature_type counts per features file, rows of a type outside the assay, hashtag-like names,
  column counts, barcode lines with more than one column;
- (d) `whitelist_sized` (GSE230531): whether each barcodes file has exactly 737,280 lines (the 10x 3' v2 whitelist,
  `screens.WHITELIST_V2`; the 3' v3 size of the port is dropped: the row names v2 only);
- and the yaml `screen:` block to fill (`verdict: READERS_VERDICT`, the files' names and sha256, `READER_1`,
  `READER_2`, the date: placeholders the yaml check refuses).

The two readers decide by their checklist (points 1 to 7 and (a) to (e)); a reader who doubts does not pin. The main
session writes, per atlas, either `screen: {verdict: kept, files: {GSM: {barcodes: {name, sha256}, features: {name,
sha256}}}, read_by: [two readers], date}` or `status: dropped` with `screen: {verdict: dropped, failed: [points]}`.

## Steps

BEFORE THE ROW'S COMMIT, on the DGX in a campaign checkout (`run_dgx.sh CHECKOUT MODE`, modes plan, fetch-meta,
screens, derive; outputs under `~/v/scratch/t72_screen`):

1. `fetch.py plan ACC --dest ~/v/scratch/t72_screen`: the family SOFT only (a bounded retry with a pause on HTTP 403,
   429, 5xx and timeouts; one host over https; no redirect followed). Writes `meta/<ACC>/plan.json`,
   `plan_table.tsv`, `records.txt`, `records/<ACC>.txt`, `processing_text.txt`, `protocol_field_sha256.tsv`; prints
   the summary and the per-field hashes. Exit 2 before pins (expected: "not pinned", "named samples not pinned").
2. `fetch.py fetch-meta ACC --dest ~/v/scratch/t72_screen --plan-sha256 <sha256 of meta/ACC/plan.json>`, then (or
   again) `screens.py ACC --dest ~/v/scratch/t72_screen`.
3. Two readers read and hand back; `fetch.py pin-template ACC --dest ~/v/scratch/t72_screen --gsm ...` prints the
   `pinned:` block and the `processing_cleared` line for the named samples.
4. `groups.py derive --released-e85 data/raw/T1/E8.5_RNA.h5ad --out experiments/t1_ext/results/extra_competing.json`
   (the released E8.5, a training stage; target-free; `--out` must be the path the yaml names).
5. THE PLACEMENT COMMIT (the main session; the row's commit): the pins, the clearances, the screen blocks, the derived
   file and its sha256 in the yaml, `status: frozen`, `row_t72.FROZEN_YAML_SHA256` set to the frozen yaml's sha256,
   the register amendments naming row T72, and row T72 in the ledger, in one commit.

THE ROW'S OWN DATA STEPS, at the row's commit, started by `experiments/t1_m1/mx2s.sh` (each only after the one before
exited 0; `--commit` is the row's commit, full sha):

6. `fetch.py fetch --commit C`: the guards; creates `data/external/t1_ext_t72/` exclusively; per kept member, the plan
   made again from the live metadata (into `meta/` there) must pass under the frozen yaml; then the pinned samples'
   barcodes, features and matrix files, each from its own folder (`raw/<ACC>/<GSM>/`, `raw/<ACC>/MANIFEST.json`,
   each barcodes and features file marked `screened_copy` true or false); writes `FETCH.json`.
7. `load.py --commit C`: the guards; creates `data/external/t1_ext_t72/reduced/` exclusively; per fetched member:
   window removal first; only pinned samples whose records still hash to their pins; each barcodes and features file
   with the sha256 of the screened copy the yaml records; checks (a), (b) and (d) again on the fetched files; Gene
   Expression rows only; the cell call on each cell's Gene Expression rows; every pinned sample loaded. Any failure is
   a failed (p1) for that member, recorded in `reduced/LOAD.json` (the step still exits 0; the pre-gate reads it).
8. `tables.py --commit C --released-dir data/raw/T1`: the guards; creates `data/derived/t1_ext_t72/` exclusively;
   labels and member tables (below) for each member whose load passed; `MANIFEST.json`.

THE GUARDS (`row_t72.guards`, row T72's cell 8), before any download or file open, each step: row T72 in the ledger;
the register strings (X002 'consensus member of T1 candidate MX2s' and 'row T72'; X019 and X020 'Candidate MX2s' and
'row T72'; X010 'skeletal-muscle competing class of row T71'); the yaml's sha256 equal to `FROZEN_YAML_SHA256`, status
frozen, equal to its version at HEAD; the derived file's sha256 the yaml's; HEAD the given commit, which must be the
commit that adds row T72 to the ledger (its parent's ledger has no row T72), on a clean tree; the step's output folder
absent and inside the checkout's own `data/` (a real directory, not a link); and the literal never-fetched list on
every request. A folder left by a stopped data step blocks a rerun (a rerun is a new row).

## The tables (section B's format)

Labels: `feng_relabel.relabel` (imported, never copied) on row N54's frozen yaml unchanged (log1p CP10k; control sets
of 25 bins and 10 genes per set gene, seed 20261119; min_genes_found 2; min_score 0.3, min_margin 0.25, split_margin
0.1; its label-to-group map), plus the competing class skm (Tnnt3, Tnni2, Mylpf, Myl1, Tnnc2; label r_skm; a set gene
is never a control gene) on both new atlases, plus, on X020's E8.5 samples only, the derived neural, paraxial and EXEM
classes. X019 is labelled in one run over both ends (as N54); X020 in two runs (E8.5 with skm and the derived classes,
E16.5 with skm only), since the derived classes are the E8.5 samples' only.

Per atlas a and group g (X019: CM, Endocardium, Endothelium; X020: CM) with at least 10 cells at each end and at least
one in every sample of the end: m_end = the equal-weight mean over the end's samples of each sample's mean over the
group's cells of the cell's CP10k (counts x 10^4 over the cell's total on ALL its Gene Expression rows, `total_ge` in
`cells.csv`), no downsampling; D = log((m_late + 0.5) / (m_early + 0.5)) x 5/8; mapped to T1's genes (the released
E9.5 file's gene names only) by unique symbol (a symbol T1 or the atlas holds twice gets no value); ^Hb[ab]- and ^mt-
genes NaN in D.

- `data/derived/t1_ext_t72/tables/<ACC>/<group>.npz`: `genes` (T1 order), `m_early`, `m_late`, `D`, `n_cells_early`,
  `n_cells_late`, `samples_early`, `samples_late`, `per_sample_early`, `per_sample_late` (samples x genes).
- `data/derived/t1_ext_t72/MANIFEST.json`: the guards' record (row, commit, clean, the yaml's and the derived file's
  sha256, tasks ["T1"]), (p1) per atlas and the screens' drops, the constants, the T1 gene file's sha256; per atlas its
  register row, ends, factor 5/8, the labeller runs (each config's sha256, the sets found), the load QC per sample
  (barcodes, cells called, cells kept), the skm purity figure of (m2) for a CM member (`skm_purity_late`: the late
  end's r_cm cells, those with skm score at least 0.3, the share); per group the table's path and sha256 and the
  cells per end and per sample, or why no table was written.

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

## The known blocker of the parked branch, closed

The parked branch's tables step read a reduction left on disk (`load.read_reduced`) with its file hashes and window
checks but without a pin check, so a reduction made under earlier pins, another atlas entry, a lapsed clearance or a
later drop could be tabled. The port (c34d8a6) closed it: `read_reduced(acc, dest, cfg)` refuses unless the reduction
was made under the atlas entry in force (pins included), every sample it holds is pinned now with its record hashing
to its pin and its stage its pin's, the pins cover every named sample at both ends, the clearance holds for its
processing text, and the member is still recorded kept with the same screened files
(`tests/test_t1_ext_mx2s.py::test_the_tables_step_reads_no_reduction_made_under_other_pins`). At the row's run the
reduction folder is also created exclusively, so no earlier reduction can sit there.

## Where this differs from the texts (each refuses more, or is forced by another text)

- The guards' "output folders absent" is read per step: fetch needs both folders absent, load needs `reduced/` and the
  derived folder absent, tables needs the derived folder absent (a literal reading would let no step after fetch run).
- "HEAD is this row's commit": the commit is given on the command line and must be HEAD and the commit that adds row
  T72 to the ledger (its parent's ledger lacks it).
- "follows no redirect off it" is read as: no redirect is followed at all (a same-host redirect refuses too).
- The standard slash field names accepted are a literal list of one, `contact_zip/postal_code` (under Series_ or
  Sample_); any other key with a slash still refuses the series.
- The retry: at most 4 requests per file, 20 s apart, on HTTP 403, 429, 5xx and timeouts only.
- A group's table needs at least 10 cells at each end AND at least one in every sample of the end (a sample without a
  cell of the group has no mean to weigh equally).
- The cell call and each cell's CP10k use the cell's counts on ALL its Gene Expression rows, before duplicated
  symbols are dropped (the port called cells after the drop).
- X020 is labelled in two runs (the derived classes on its E8.5 samples only), so its controls are drawn per end;
  X019 in one run over both ends, as row N54.
- T1's gene names are read from the released E9.5 file (a reference stage of every practice task), not E10.5.
- (d) is the 3' v2 whitelist's 737,280 lines only; the port also allowed the 3' v3 size.
- A fetched barcodes or features file that is not the screened copy is recorded by the fetch and refused by the load
  (a failed (p1)), as the row's LOAD clause places it.

## Not checked

- No GSM or GSE record was opened. The stage field names, the plug-day statements, the per-sample text, the file names
  and the files' contents are unknown until the plan and the screens run. The pins guard against what the rule misses
  only as well as the two readers read.
- False refusals. The rule fails closed on ordinary text that looks like a stage token; the pooling clause refuses
  benign values that name a stage in one sentence and mix a reagent in another (listed above). Its false-refusal rate
  on real GEO records is unknown. A refusal shows its reason in the plan table; no clearance or pin can override it.
- Whether GSE244783's RNA samples' features files carry Peaks rows at all is unknown: if not, (b) says no and X019 is
  dropped by the readers.
- The hashtag-like names and the Theiler and somite tables are general knowledge, not checked against a source; the
  whitelist size is the row's.
- The labeller, the tables and the three steps ran only on synthetic data (`tests/test_t1_ext_tables.py`,
  `tests/test_t1_ext_row.py`); the labeller's set sizes on the real atlases (min_genes_found) are unknown.
