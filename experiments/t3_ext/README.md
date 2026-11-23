# T3 outside-knockout factor tables (stage 9 round 2: bulk whole-embryo units, phases 1 and 2, on stage 8)

Campaign branches `campaign/t3-outside-knockouts` (stage 1), `-2` (stage 2), `-3` (stage 3), `-4` (stage 4: series
membership, files bound to their pin, the series text shown and pinned), `-5` (stage 5: file urls
derived from the hashed fields only, links only from the sample's own folder, no stale output, a parser that splits
on newlines only, plain field names that are screened) and `campaign/t3-ext-archives` (stage 6: a series' own
gene-named archives or loose matrices, pinned to the samples they hold with a split rule, and controls from a sister
series of the same study; sections "Series-level files" and "A control from a sister series" below) and
`campaign/t3-ext-3` (stage 7: the three blockers and the lesser points of the two reviews of stage 6 closed, section
"Stage 7" below, and an intake for per-sample unfiltered matrices, section "Unfiltered per-sample matrices") and `campaign/t3-ext-4` (stage 8:
the three blockers and the lesser points of the two reviews of stage 7 closed, section "Stage 8" below) and
`campaign/t3-ext-bulk` (stage 9: bulk whole-embryo knockout units, section "Bulk whole-embryo units") and
`campaign/t3-ext-bulk-2` (stage 9 round 2: the six blockers and the lesser points of the two reviews of stage 9
closed, section "Stage 9 round 2" below). Everything of stage 5 stands. Under the captain's ruling of 2026-11-20 (`docs/PLAN.md` decision log
2026-11-20 (j), "All three tiers"), outside single-cell data on knockouts of OTHER genes may shape the Task 3 file.
This folder turns such data into per-gene factor tables in the frozen factor's form, so that a later thin builder can
take one in place of the Mab21l2 factor. Stage 1 builds the tables only. Nothing here enters a model, a check or an
upload until (1) each dataset has its row in `EXTERNAL_DATA.md`, (2) try K21's row is committed before any table is
scored, and (3) the captain rules on any change to an upload, with a `docs/S2_RULES.md` amendment.

## Files

| File | What it does |
|---|---|
| `datasets.yaml` | The datasets, their tier, study unit, design, control, stages, sample selectors (suggestions only), status, and per cleared dataset the `pinned` block. Its sha256 is pinned in `t3x_common.DATASETS_SHA256`. |
| `t3x_screen.py` | The held-out-gene screen, the second line (Ctnnb1 / beta-catenin, Gata4, their drugs; review terms such as Gata6, Apc, Axin, Gsk3). It has no switch. |
| `t3x_common.py` | Fixed constants, `NEVER_FETCH`, `DATASETS_SHA256`, the dataset loader and the pinned-block check, the canonical metadata serialisation, the stage cap (E9.5 inclusive for every cell), shrinkage, the outside-the-repository guard. |
| `t3x_tables.py` | Finishing a log table into factors (rescale, clip, absent = 1), writing the pooled and per-group forms, the reference statistics. |
| `fetch.py` | `plan` writes a plan table of every sample of every screened series from GEO's family SOFT; `pin-template` prints the entries of named samples (with `--from-series-file` for samples whose cells are in a series file), or with `--series` the series entries, to paste; `series-file-template` prints one pinned series file entry (stage 6); `check-pins` is the gate; `datasets-sha` prints the hash; `fetch` downloads exactly the pinned files of the pinned samples and the pinned series files. |
| `process.py` | `reference` writes the WT E8.75 statistics. `run` turns the pinned samples of one fetched dataset (their own files, or their pinned series file split by its rule; since stage 7 a declared unfiltered matrix through the one fixed cell call) into one table per knockout gene plus `QC.json`. |
| `combine.py` | Tier tables (wnt, cardiac, unrelated), the all-knockouts table, one table per study, and leave-one-study-out tables. |
| `run_dgx.sh` | The detached DGX launcher, with exit files: mode `plan` (the default), `header` (stage 9: the first line of a bulk count table) or `all` (behind the pinned gate). |
| `bulk.py` | Stage 9: bulk whole-embryo units from GEO (per-sample text tables, one series text table): the file rules, the header step and the exact tie rule, the bulk pins, fetch, the readers, the id routes, the math and the guards G1 to G12, `run_bulk`. `fetch.py`, `process.py` and `combine.py` dispatch to it for a dataset with `modality: bulk`. |
| `ae.py` | Stage 9, phase 2: the many-line ArrayExpress resource E-MTAB-7322 (IDF, SDRF, methods text, the accession-keyed count matrix): withheld lines, records per Source Name, the views, the mechanical comparisons, the committed line-set rule, the run ties, its pins. |
| `make_panel_ids.py` | Stage 9: writes `panel_ids.tsv` (the gene-id map of the Ensembl route) from the features file of an already pinned and fetched 10x sample; pinned by `t3x_common.PANEL_IDS_SHA256` (None until the main session commits the file). |
| `chipq_fetch.py` | T3 candidate CHIPQ-S (register row X018, standard library only): `lists`, `rows`, `targets`, `verify` against a pinned generic ChIP database (two hosts, no redirect followed, mm10, the two lists and the +/-5 kb target-gene files only), a manifest with URL, size, sha256 and fetch time per file, nothing deleted or overwritten (a partial download is set aside as `.failed.<n>`; a failed download is retried up to three times; a short pause before every request); no message names a target-gene file by its path; the committed keyword filter for placebo antigens (`DROP_PATTERNS`). Data root `~/v/vec-ext/t3/chipatlas`. |
| `chipq.py` | CHIPQ-S's database part: the kept list (`chipq_kept.yaml` with `CHIPQ_KEPT.sha256`, both committed by the main session after the two readers' screen; refused with exit 2 until then), the one reader of a target-gene file (`open_target`), s_g, the bound set, the weights, the raw table exp(-w) finished through `t3x_tables.finish`, the exponent at S* and 4 S*, the placebo pool (at least 10 antigens) and the per-seed draw without replacement; steps `p1` (two integers, exit 0 or 3) and `tables` (`chipq_Gata4` and `chipq_rand_seed0`-`9`, never Ctnnb1's, into the folder `VEC_T3_EXT_TABLES_K27` names, `~/v/vec-ext/t3/tables_k27` on the DGX). Every reading is one named constant in its READINGS block, as the rule writer's settled points fix it. |

## The pinned sample list (the first line)

A text blacklist can never be shown complete, so it is not the first line of defence. The data in scope were five
cleared series (GSE205117, GSE168768, GSE186934, GSE137337, GSE204906), a few dozen samples in all, and the first line
is an explicit list of the samples two people have read. GSE205117 is held since stage 8 (its review term was not
cleared by the two readers, decision log 2026-11-20; it is on the team's held list, `t3x_common.HELD_EXTRA`).

1. **Plan** (metadata only). `fetch.py plan ACC` writes `meta/<ACC>/plan.json` and one table per screened series,
   `meta/<ACC>/plan_table_<SERIES>.txt`. Each table starts with a **SERIES block for every screened series**: its
   accession, its `series_metadata_sha256` and every field of its record verbatim (title, summary, overall_design,
   type, relations, sample list, contact fields: all shown, hashed and, since stage 5, screened). Then every
   sample of the series: its accession, the verbatim text of every field of its GEO record (all are screened, and
   since stage 5 every field NAME is screened too), the sha256 of those fields, the screen's verdict and hits, its
   **series membership** (`NOT PINNABLE` when its own `series_id` lists a series outside `screen_series` or none, or
   when a supplementary link names a held accession), the per-sample files allowed and every other supplementary link
   with the reason it is never fetched, the parsed stage (with the field(s) it was read from, or why it is `none`)
   and chemistry, and the selectors' suggested arm, which chooses nothing. **Stage parsing:** by default the stage
   is read from the characteristics (keys naming stage, age, time or day), else the title, as before. A dataset may
   list `stage_fields`, an ordered list of SAMPLE field names (plain keys: ASCII letters, digits, underscores; only
   GSE204906 has one, `[characteristics_ch1, source_name_ch1, description]`, because its samples state the stage only
   in the source name and the description). Then each listed field is read by the same parser (the characteristics
   through their stage keys, any other field through every value; the title only if listed); a field the record lacks
   states no stage; every field that states one must give the same single value; no field stating one, a range, or
   two values give `parsed stage none` with the reason, and such a sample is not pinnable. The E9.5 cap applies to
   the agreed value. plan, pin-template and verify_pins (so fetch and process) all use this one rule
   (`fetch.sample_stage` under `t3x_common.stage_fields`); the screen does not depend on `stage_fields`. **A file is allowed only when its link is the sample's own GEO folder**, exactly
   `ftp.ncbi.nlm.nih.gov/geo/samples/<stub>/<GSM>/suppl/<GSM>_<name>` (GEO's stub: GSM1234567 -> GSM1234nnn; no port,
   query or fragment; a plain file name), names no held accession, is not a denied type and matches the dataset's
   `sample_file_patterns`. A link into a series folder, another sample's folder or another host is listed as never
   fetched with its reason. Under any value that is not plain printable ASCII the table prints a second line, `escaped:`, its
   Python `ascii()` form, so a hidden, reordering or look-alike character shows. Lines of a sample or series block
   that are not `key = value` fields are kept under `_other_lines` (shown, hashed, screened); the header counts the
   SOFT lines outside the series and sample blocks (PLATFORM and DATABASE blocks), which are not shown. A dataset
   with a held-out hit anywhere, or with a sample of a screened series whose `series_id` also names a held series, is
   refused and none of its text is written: the table holds accessions, verdicts and the names of the fields hit, and
   nothing can be pinned from it. SOFT bytes that are not UTF-8 refuse the dataset. Since stage 5 the SOFT text is
   split on the newline character only (one trailing carriage return stripped per line), and the dataset is refused,
   its text withheld, when a line still holds a line or paragraph separator (CR, VT, FF, FS, GS, RS, NEL, U+2028,
   U+2029), when a `!Series_` or `!Sample_` field name is not a plain GEO key (ASCII letters, digits, underscore; the
   one exception is GEO's `contact_zip/postal_code`), or when a SAMPLE or SERIES block's accession is not `GSM` or
   `GSE` digits. `plan.json` says `in_progress` from the start of a run, so a run that dies leaves no passing plan
   behind.
2. **Read and pin.** The main session and a second, independent reader each read the plan tables (what they check is
   below). For each sample both agree on, `fetch.py pin-template ACC --dest D --arm ko|control --gene G --gsm GSM ...`
   prints the entry (one command per arm; every sample is named, there is no bulk selector), with every field of the
   sample as a comment above it, and it is pasted into the dataset's `pinned` list in `datasets.yaml`:
   `{accession: GSMnnn, arm: ko|control, gene: <knockout gene>, stage: E8.5, files: [...], metadata_sha256: "..."}`.
   `fetch.py pin-template ACC --dest D --series` prints the `pinned_series` block, one entry per screened series,
   `{accession: GSEnnn, series_metadata_sha256: "..."}`; it is pasted beside `pinned`, and is required whenever
   `pinned` is not empty. A control names the knockout gene it controls; a control that serves several knockouts is
   pinned once per gene. pin-template never writes `datasets.yaml` and refuses a held-out sample, a NOT PINNABLE
   sample, a sample whose stage does not parse to one value (or whose stored stage differs from the one its fields
   give under the committed `stage_fields`), and a sample with no allowed per-sample file; it marks a
   sample with review terms (pin it only together with the `review_cleared` entries both readers agree).
3. **Commit.** Every edit of `datasets.yaml` changes its sha256. The same commit puts the new value
   (`fetch.py datasets-sha`) into `t3x_common.DATASETS_SHA256`; otherwise every step refuses the committed file, and
   `tests/test_t3_ext_screen.py::test_the_committed_dataset_list_hashes_to_its_pin_in_t3x_common` fails with the
   reason.
4. **Fetch and process** touch only pinned accessions. They refuse a dataset whose pinned block is empty or has no
   `pinned_series` for every screened series; a plan that did not pass or was made from another `datasets.yaml`; a
   pinned series whose live text (when plan runs again) or stored text (fetch, process) no longer hashes to its pin;
   a pinned accession that is not in its own series' plan table or belongs to another series; a pinned sample whose
   own `series_id` field is missing, or lists a series outside the dataset's `screen_series`, a `NEVER_FETCH` series,
   or the accession or a screened series of any dataset that is not cleared; a pinned sample whose live or stored
   metadata no longer hashes to the pinned value; a pinned sample that no longer passes the screen; a pinned stage
   that differs from the record's; a pinned file not in the sample's plan entry; knockouts and controls whose stages
   do not pair, or whose chemistry differs where the dataset reads chemistry; a manifest fetched under other pins;
   a processed QC made under other pins (combine). **Files are bound to the pin they were fetched under:** each
   manifest file record holds the sample's pinned `metadata_sha256`, the url and the `pins_sha256` (samples and series)
   in force at download time. fetch reuses a file only when all three equal the current values, otherwise it deletes
   it and downloads again, and it deletes any other file in a pinned sample's folder. process refuses a file whose
   record holds another `metadata_sha256`, url or `pins_sha256`, and any file in a pinned sample's folder that is not
   one of its pinned files with a record. Arms, genes and stages come from the pins only. No count file of an
   unpinned sample is opened (a test holds one beside the pinned ones and fails on any open of it). **Since stage 5
   the urls come from the hashed fields only:** verify_pins derives each pinned sample's file list (names and urls)
   again from its pinned fields, refuses a stored list that differs and two files of one name, and fetch and process
   use only the derived list. fetch and process refuse a symlinked raw, dataset or sample folder, manifest or file.
   **`raw/` and `MANIFEST.json` are written by `fetch` only.** Nobody copies, edits or moves a file in them; to redo a
   sample, delete its folder and fetch again. **No stale output:** process first deletes `processed/<ACC>/`, and
   writes into `QC.json` the sha256 of the plan it ran under (`plan_sha256`), of the manifest and of each count file;
   combine takes `--dest` and refuses a QC unless `meta/<ACC>/plan.json` exists, passes (status ok, the committed
   `datasets.yaml`, its pins verified) and hashes to that `plan_sha256`; the launcher removes `processed/<ACC>/`
   before the accession's re-plan, so a refused re-plan, fetch or process leaves nothing for combine.
5. **Per-cell genotype designs** (cell-level labels inside one sample: MULTI-seq, gem.group, chimera sorts) are
   excluded in this stage: a dataset with `cell_labels` cannot hold pins, and `process.py` refuses it.

## Series-level files (stage 6)

Some studies put their matrices only in the SERIES' supplementary list, named by gene (GSE137337: archives such as
`GSE137337_RNF2_E8.5_1ab.tar.gz`, and three genes as loose barcodes / features / matrix files). Stage 6 lets such a
file be used without loosening any rule above.

1. **Declared.** A dataset's `series_files` lists, per named series (one it screens: its own, a superseries or a
   named control series), the patterns a file must match to be OFFERED. A file is offered only when its link is that
   series' own GEO folder, exactly `ftp.ncbi.nlm.nih.gov/geo/series/<stub>/<GSE>/suppl/<GSE>_<name>` (GEO's stub:
   GSE137337 -> GSE137nnn; no port, query or fragment; a plain file name that starts with the series accession),
   names no held accession, is not a denied type (the series' `_RAW.tar` is denied) and matches a pattern declared for
   that series. The series' supplementary list is a field of its record, so it is already shown, screened and hashed
   by the series pin. Each SERIES block of the plan table prints `series files allowed` and every other series file
   with the reason it is never fetched.
2. **Pinned.** `pinned_series_files` maps each file the readers take to the pinned samples whose cells it holds:
   `{series: GSEnnn, files: [one archive, or the three files of one loose 10x triple], samples: [GSM...], split:
   {by: none} | {by: barcode_suffix, map: {"-1": GSMa, "-2": GSMb}} | {by: folder, map: {"RNF2_1a": GSMa, ...}}}`.
   `fetch.py series-file-template ACC --dest D --series GSE --file F [--file F2 F3] --gsm GSM ... --split S [--key
   K=GSM ...]` prints it from the stored plan (the series' offered files derived again from its hashed fields); the
   samples are pinned as usual with `pin-template ... --from-series-file`, which prints `files: []`. The loader
   refuses: a file the series_files block does not declare for that series; anything but one archive or one triple;
   a file or a sample in two entries; a sample of an entry pinned with its own files, or a sample with `files: []`
   that no entry names; a split map that does not name every sample of the entry exactly; a suffix that is not one
   separator and letters or digits (`-1`); a folder that is not one plain path component; a folder split of a loose
   triple. **When a file cannot be split by a rule stated in the pin** (split none), every sample of the entry must be
   pinned with ONE arm and ONE stage, and for the same genes (a wild-type file may serve several knockouts); the file
   must hold exactly as many libraries as the entry names samples (stage 7, checked by process before any count is
   read); the file is then used as one unit and counts as ONE library (m), and its confidence share is recorded
   under one label `series:<GSE>/<file>` (with one sample, under that sample's accession).
3. **Checked again.** `verify_pins` (plan with pins, fetch, process, combine) derives the series' offered files again
   from the pinned series fields, refuses a stored list that differs, two links of one file name, a pinned file that
   is not offered, and a sample of an entry that is not in that series' plan table.
4. **Fetched.** `fetch` downloads exactly the pinned series files into `raw/<ACC>/<GSE>/`, each manifest record bound
   to the series' pinned text hash, the url and the pins (`pins_sha256`, which takes the series file entries and a
   named control series in; a dataset without either keeps its old hash). Any other file in that folder is deleted.
5. **Read.** `process` reads only the pinned files, after checking each manifest record and that the series folder
   holds nothing else. An archive is never extracted: it is read member by member in memory. Before any count is
   read, every member is checked, and the file is refused for an absolute path (or a drive letter or backslash), a
   `..` component, an empty or `.` component, a name part outside letters, digits and `. _ + -`, a symlink, a hard
   link, a device or FIFO, a duplicate name, a member name the held-out screen hits (held out or review), a member
   larger than 8 GiB, and **a file member that is not part of a 10x triple (a member it has no rule for)**. Split
   folder: every file member lies in a top-level folder the map names, and each folder holds exactly one triple.
   Split barcode_suffix: every barcode ends with exactly one key (with something before it), and every key's sample
   gets cells; a barcode no key fits refuses the file. Split none: each folder of the file holds exactly one triple,
   the file holds exactly as many libraries (folders with a triple, or the one loose or flat triple) as the entry names
   samples, else it is refused ("the file holds N libraries but the entry names M samples"), and all its libraries
   are pooled under one label. A member whose path is a denied type (an unfiltered
   `raw_feature_bc_matrix`, fragments, peaks, ATAC) is refused too. More than 60,000
   barcodes in one library (each folder, or the one loose triple, under every split) is refused as unfiltered. A
   gzipped member is decompressed in chunks and refused once it passes 8 GiB (stage 7). `QC.json` records each file's libraries (label, folder, barcodes)
   and its sha256 under `count_files_sha256` as `<GSE>/<file>`.

## A control from a sister series (stage 6)

GSE137337 has no wild-type single-cell sample: the study's wild-type embryos are in another series of the same lab.
A dataset may name that series as `control_series`:

```yaml
control_series:
  accession: GSEnnnnnn             # null: a placeholder, inert (screens and pins nothing)
  screen_series: [GSEnnnnnn]       # its own accession first, then any superseries; none shared with screen_series
  relation: {same_lab: "...", same_protocol: "...", same_publication: "..."}   # each a statement the readers checked
  reader_note: "..."               # required: what both readers read and agree
```

The sister series and its screen_series are planned, screened, shown, hashed and pinned exactly as the dataset's own:
`plan` downloads their family SOFT, the screen runs on every sample and series field, every plan table starts with
every SERIES block (marked `(control series)`) and the control series declaration, the sister's own table says
`THIS IS THE CONTROL SERIES` and lists every sample with its series membership, and `pinned_series` must name every
screened series of both. Its samples are pinned by accession like any other, as controls only (`pin-template --arm
control` offers them; a knockout pin of a sister sample is refused). A placeholder or any field still holding `TO FILL`
cannot be used; a sister series that is held (NEVER_FETCH or a held dataset's series) refuses the dataset list, and a
held dataset's sister series joins the held set. The plan stores the declaration; a plan made under another
declaration is refused. **d**: a knockout unit with any control from the sister series takes d = 0.5
(`t3x_common.SISTER_D`, the weaker value the pipeline gives controls that are not littermates or same-capture), whatever
the dataset's design; the unit's QC entry and unit file say `d` and `control_source: sister series <GSE>
(control_series): [GSM...]`, and `QC.json` lists `sister_series_controls`.

**Checked by code, not on typed text alone (stage 7).** (1) `load_datasets` refuses a control series, or any series of
its screen_series, that is another dataset's accession, one of another dataset's screen_series, or a series of another
dataset's control_series: a sister series belongs to no other entry. (2) `plan` refuses the dataset (status refused,
the reason stated) unless the dataset's own series record and the sister series' record share a `contact_institute`
value AND at least one value of `contact_name`, `contact_laboratory` or `pubmed_id`. Since stage 8 (the first
reviewer of ec86a1b, blocker 3) it also refuses when both records carry `pubmed_id` values and share none (two
publications: not one study), and otherwise requires a shared `pubmed_id` OR a repository link between the two series:
one record's `relation` names the other series (`SuperSeries of: GSEnnn` or `SubSeries of: GSEnnn`, GEO's exact form),
or both name the same SuperSeries (both `SubSeries of` one series). A shared institute and contact alone no longer
suffice; an unpublished sister series with no repository link is refused, by design. Values are compared exactly as
the SOFT parser gives them (whitespace stripped, case kept), the stricter of the readings the reviewer's text allows.
The values matched and the link (`link = ...`, which of the links matched) are printed on the `CONTROL SERIES` lines
of every plan table (`the code matched ...`), stored in the plan, and derived again from the stored series fields by
`verify_pins` (fetch, process, combine). (3) A series of
control_series.screen_series that holds a sample of the dataset's own screened series refuses the plan: a superseries
that holds samples of the dataset's own series is listed once, in the dataset's screen_series, never under
control_series; list under control_series.screen_series only series that hold the sister's samples alone. The match
says the two records share an institute, a contact or a publication, and a publication or a repository link; it
does not say the protocol is the same or
that the wild-type embryos fit the knockouts: the readers still check that (check 11).

## Unfiltered per-sample matrices (stage 7)

Some studies' only per-sample files are the vendor pipeline's UNFILTERED output (GSE270969: nine single-embryo
libraries, each sample with one `raw_feature_bc_matrix` archive of its own). Such a file holds every barcode the
pipeline saw, empty droplets included. Stage 7 takes it through one fixed cell call, without loosening any rule above.

1. **Declared.** A dataset may declare

   ```yaml
   raw_matrices:
     patterns: ['^GSM\d+_[A-Za-z0-9_.+-]*raw_feature_bc_matrix\.tar\.gz$']   # a sample's OWN file name (case kept)
     reason: "..."                                                             # why the raw files are taken
   ```

   Only for a dataset with that declaration, and only for a sample file of the sample's own folder that a declared
   pattern matches, the deny term for unfiltered matrices (`raw_feature_bc_matrix`) is lifted (`fetch.denied`); every
   other deny term stays (fragments, ATAC, peaks, BAM, the series' `_RAW.tar` ...), the file must still match
   `sample_file_patterns`, and a name a pattern misses stays denied. A raw file is taken only when it is the sample's only file of its
   own; a sample that lists any other file of its own never has its raw file fetched (stage 8; "of its own" is any
   link into the sample's own GEO folder, whatever the reason that link is or is not offered). A sample's pinned
   files are all raw or all filtered, and since stage 8 so are all the samples of one unit: `verify_pins` (plan with
   pins, fetch, process, combine) refuses a unit whose knockout and control samples are not all unfiltered matrices
   through the fixed cell call or all filtered (a sample of a pinned series file counts as filtered), and `process`
   refuses the same as a second guard. Pins, the screen, `NEVER_FETCH`, the derived held set and the dataset-file hash work
   as before; the plan stores the declaration and a plan made under another is refused; `pins_sha256` takes the
   declaration in when it is present (a dataset without it keeps its hash). Only GSE270969 declares it, with no pins.
   The plan table prints the declaration and the cell call at its top, and under each sample with a raw file an
   `UNFILTERED ...` line; pin-template marks such a sample.
2. **One fixed cell call**, the same for every sample of every dataset, as constants in `t3x_common` (`RAW_CALL_*`,
   `RAW_MIN_CELLS`; not configurable per dataset, from the command line or the environment). Per library: each
   barcode's total UMIs over gene-expression features; N = 3000; m = the 99th percentile (numpy's default linear
   interpolation) of the N largest totals (all of them when the file holds fewer); a barcode is called a cell when
   its total is at least max(500, m / 10). Then the pipeline's ordinary per-cell rule (200 detected genes) and its
   cap on barcodes (60,000) apply to the CALLED cells; a library with fewer than 500 called cells, or more called
   cells than the cap, is refused. `QC.json` holds `raw_cell_call`: the rule's constants and, per raw sample, the
   barcodes in the file, the barcodes with any count, m, the cutoff, the cells called and the share of all
   gene-expression UMIs that the called cells hold (with the feature and entry counts); since stage 8 also the 5th,
   25th and 50th percentiles of the called cells' gene-expression totals and the number of called barcodes whose
   total is under twice the cutoff. Each unit's QC says `raw_called` and `raw_called_samples`, and since stage 8 a
   raw-called unit's QC holds `raw_call_arms`: the per-arm medians of the cutoff, the cells called and the UMI share in
   called cells, and their knockout-to-control ratios. These are records only: no guard and no threshold reads them,
   and a call difference between arms is for the readers of the QC to note, not a refusal. Nothing else downstream (unit files, tables, combine) tells a raw-called
   library from a filtered one.
3. **Read sparsely.** A raw matrix can hold millions of barcodes and of the order of a hundred million entries. It
   is read in two streamed passes, never densified and never unpacked: the MatrixMarket header and size line, then
   the entries in blocks of 64 MiB of text; pass one sums each barcode's gene-expression total, pass two keeps only
   the called barcodes' entries. Refused: a file above 16 GiB (`RAW_MAX_FILE_BYTES`), a header stating more than 400
   million entries (`RAW_MAX_ENTRIES`), more than 20 million barcodes (`RAW_MAX_BARCODES`), a header that does not
   match the features and barcodes files, an entry out of range, negative or not finite, a line over 1 MiB, more or
   fewer entries than the header states, a 10x h5 file (the intake reads MatrixMarket triples only), and since stage
   8 a matrix line that is not exactly three fields and more than 1,000 comment or empty lines before the size line. A per-sample
   archive is read with the member rules of the archive routes (no absolute path, no parent step, no link or device
   member, no name the screen hits, exactly one barcodes, features and matrix triple in one folder, nothing else;
   for the members alone the unfiltered-matrix deny term is lifted). A leading `./` is accepted, and the triple's
   folder may have any plain name: a `filtered_feature_bc_matrix/` folder inside a raw-named archive is read through
   the fixed call like any other (harmless: the call keeps its cells). A loose raw triple of three files is read the
   same way.
4. **What it does not do.** Without the declaration, an unfiltered matrix whose name lacks `raw_feature_bc_matrix`
   (for example `GSMx_raw_matrix.mtx.gz`) goes through the filtered route: the barcode cap (60,000 per library) is the
   only mechanical guard against an unfiltered matrix under another name. It does not tell nuclei from cells, remove ambient RNA, detect doublets or check that a
   raw file is one library of one sample: the readers check that from the record (check 12), and the per-sample
   confidence floor and the guards act downstream as for any library. The cell call was fixed by the design brief
   before any raw count was read; it was not tuned on any data.

**The canonical serialisation.** A pin's `metadata_sha256` is the sha256 of the UTF-8 bytes of

```python
json.dumps({"accession": GSM, "fields": F}, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
```

where `F` maps every field name of the sample's GEO family SOFT record (the `!Sample_` prefix removed, so `title`,
`characteristics_ch1`, `supplementary_file_1` ...) to the list of its values, each a string exactly as
`fetch.parse_soft` gives it (surrounding whitespace stripped, nothing else changed), in the order of the SOFT file
(a key's values together, keys in the order they first appear), plus the key `_other_lines` holding the record's
other non-empty lines when it has any (`t3x_common.canonical_metadata`; `tests/test_t3_ext_pipeline.py` pins the exact
bytes for one example). The hash is taken from the sample's record in its own series' SOFT file. A series'
`series_metadata_sha256` is the same serialisation with the series accession and `F` its fields with the `!Series_`
prefix removed, contact fields included (`t3x_common.series_metadata_sha256`; the exact bytes are pinned by a test).

**Held datasets.** `t3x_common.NEVER_FETCH` holds GSE165300, GSE126128, GSE280543, GSE280587 (GSE280543's second
screened series, added in stage 4), GSE259342 and E-MTAB-7325. Since stage 8 `t3x_common.HELD_EXTRA` holds exactly the
team's held list (82 accessions, `NEVER_FETCH` among them; its hash is pinned by `tests/test_t3_ext_held.py`). The held
set the code checks is `NEVER_FETCH` and `HELD_EXTRA` together with the accession and every screened series of every
dataset that is not cleared (`t3x_common.held_series`), so a held series missing from the constants is still held; a
test pins that they agree for the committed file. A HELD_EXTRA accession is refused as a dataset (`dataset()` and, for
a cleared one, `load_datasets()`), as a screened series and as a control series, and a sample whose record lists one
is not pinnable. GSE205117 is on that list and registered: its status is `held_phenocopy` (the status the pipeline
uses for held datasets), its reason "the two readers did not clear its review term (decision log 2026-11-20)".
`dataset()` refuses them whatever a YAML says, `load_datasets()` refuses a YAML that gives any of them a status other
than `held_phenocopy` or `excluded_design` and a cleared dataset whose `screen_series` holds any held series, and
combine refuses any processed output of them. No CLI takes a
`--datasets` path: every step reads the committed file, checked against `DATASETS_SHA256`. Tests reach other YAML
content through a function argument only.

**macOS companion members (2026-11-20).** The real archives of the pinned study were made on macOS and hold, beside
each real member, a 176-byte AppleDouble companion (`._name`) and a Finder `.DS_Store`. `process.archive_members` never
reads such a member and never classes it as part of a triple (so `._matrix.mtx.gz` is not a matrix). It skips a `._X`
member only when it is at most 4,096 bytes and the member `X` (a file or a folder) it names is in the archive beside
it, and a `.DS_Store` member only when it is at most 65,536 bytes; all skipped members of one archive together hold at
most 1 MiB; anything else without a rule still refuses the archive. Every skipped member first passes the path, link,
character, screen and deny checks like any other. On the unfiltered route (`read_raw_sample`) the number skipped is in
the sample's cell-call record (`mac_companions_skipped`); the per-sample and series archive routes skip the same
members without a record.

## What the two readers check in a plan table before pinning a sample

Each reader reads the whole plan table of the dataset's own series and of every other screened series (the
superseries), the SERIES blocks at its top included, not only the rows the selectors suggest, and pins a sample (and
the series) only if both agree on all of these:

1. **Genotype and arm.** The sample's own text (title, source, characteristics, description, protocols, file names)
   states its genotype, and it is the arm being pinned: a knockout of the named gene, or that knockout's control as
   `control_definition` describes it. A label that names no gene (a bare 'F/F', 'cKO', 'mutant') is not pinned unless
   another field of the same record names the gene.
2. **No held-out or phenocopy label anywhere in the series.** Nothing in any SERIES block (title, summary,
   overall_design, type, relations, `_other_lines`) and no sample of any screened series mentions beta-catenin,
   Ctnnb1, Gata4, another allele of either, a drug or construct acting on them, or a perturbation of their pathway
   that could phenocopy them (Wnt, Wls, Apc, Axin, Gsk3, Tcf/Lef, Fzd, Dvl, Csnk1a1, Hand2, Mef2c, Zfpm2, any GATA
   factor; Tbx5, Baf60c/Smarcd3, Sox17, Foxa2 (Gata4 partners); Tbx6, Sp5/8, T/Brachyury, Mixl1, Cdx2 (Wnt targets);
   Evi, Notum, ICAT, Chibby, Jup/plakoglobin), in any spelling, in a field's value or in a field's NAME. The
   screen's verdicts help; they do not decide. If in doubt, nothing in the dataset is
   pinned and the question goes to the main session.
3. **No pooled design.** The sample is one embryo pool of one genotype: no MULTI-seq, cell hashing, CRISPR guide
   pool, chimera, mixed-genotype capture or per-cell genotype label (those are excluded in this stage). Read
   `extract_protocol`, `data_processing` and the series' overall_design in full: pooling, hashing or demultiplexing
   (HTO, CMO, MULTI-seq, demuxlet, souporcell, SNP, genotype calling) anywhere means do not pin.
4. **Same-study control.** The control comes from the same series and experiment, at the same stage, with the same
   chemistry where the dataset reads chemistry (GSE137337), and the relation is the one `control_relation` records.
   A control whose text, or whose series' text, cites another GEO series, another paper or "previously published"
   data is not pinned. The one exception (stage 6) is a control of the series that `control_series` names, under
   check 11.
5. **Series membership.** The sample's `series membership` line says `ok`: its own `series_id` lists only this
   dataset's screened series.
6. **Stage.** The parsed stage in the table matches what the text says, is a single stage at or before E9.5, and is
   in `stages_used`.
7. **Assay.** Single-cell or single-nucleus RNA (the GEX part of a multiome), not ATAC, bisulfite, methylation, bulk
   or spatial; the pinned files are that sample's own per-sample matrices (no fragments), or, under check 10 only, a
   pinned series file of the series' own folder.
   The `per-sample files allowed` line names only files of the sample's own folder; read the `never fetched` line's
   reasons too: a link that points into a series folder or another sample's folder is a reason to ask, not to pin.
8. **Hidden characters.** Every `escaped:` line is read: a control, format or bidirectional character, a look-alike
   letter or a non-ASCII digit is explained or the sample (or series) is not pinned. The header's count of lines not
   shown (PLATFORM and DATABASE blocks) is noted; those lines are platform descriptions, not sample labels.
9. **The hashes** in the pasted entries equal the ones printed in the table the reader read: each sample's
   `metadata_sha256` and each series' `series_metadata_sha256`. A pinned sample's row in every other screened table
   (the superseries) must read the same as its own-series row; if not, do not pin (only the own-series row is hashed).
10. **A series file (stage 6).** For each `pinned_series_files` entry, both readers check: (a) the file is listed under
    `series files allowed` in the SERIES block of its series, its link is that series' own folder, and they read the
    `other series files` line's reasons; (b) the file's name, and the titles and descriptions of the samples named,
    tie the file to exactly those samples: the gene in the file name is the gene those samples are pinned for, the
    stage in the name (if any) is their stage, and the replicate tags (for example `1ab`) match the samples' titles;
    (c) no other sample of any screened series could also be in the file (a series text or sample text that says the
    file holds more libraries than the entry names means do not pin), and the number of libraries the file is
    expected to hold equals the number of samples named (process refuses a file whose library count differs); (d) the split rule is stated by the records
    themselves: a barcode suffix only where a record (the sample's or the series' data_processing, description or
    overall_design) says which suffix is which sample, a folder only where a record names the folders; otherwise
    `split: {by: none}`, which needs every sample of the entry to be ONE arm and ONE stage, and pinned for the same
    genes (a wild-type file may serve several knockouts), and which counts the file as ONE library; (e) every
    sample of the entry is pinned with `files: []`, by `pin-template --from-series-file`; (f) the file holds expression matrices of single-cell RNA, not bisulfite, ATAC or methylation
    data (the name and the series' overall_design and data_processing say so). What the archive really holds is
    checked only by process, after the download, and a file that does not fit its rule is refused, never repaired.
11. **A control from a sister series (stage 6).** Both readers read the sister series' SERIES block and EVERY sample of
    its plan table in full, with checks 1 to 9 applied to it as to the own series, and check: (a) the three relation
    statements of `control_series` hold on the two series' records (same lab: contact name, institute, laboratory;
    same protocol: extract protocol, library construction and chemistry text, 10x version; same publication: the same
    citation or PubMed field, or the same publication named in both summaries); a relation that the records do not
    show is not written, and then nothing of the sister series is pinned; (b) the sister's samples pinned are
    wild-type (or the control `control_definition` names) at the knockouts' stage, of the same strain background and
    chemistry where the records say them; (c) nothing in the sister series' text names a held-out gene, its pathway
    or a phenocopy (check 2); (d) the reader note says what both read and agree, in plain words; (e) every unit with
    such a control will take d = 0.5, and the plan table's `CONTROL SERIES` lines match the declaration they agreed;
    (f) stage 7: the `the code matched` lines name the values the code found shared (a contact_institute and a
    contact name, laboratory or PubMed id, and since stage 8 the `link` line: a shared PubMed id or the super- or
    sub-series relation that ties the two series); both readers find those values in both SERIES blocks, and the match is a
    floor, not the relation: (a) to (d) still decide; (g) control_series.screen_series lists only series that hold
    the sister's samples alone (a shared superseries goes in the dataset's own screen_series).
12. **An unfiltered matrix (stage 7, `raw_matrices`).** For each sample whose `UNFILTERED` line names a raw file,
    both readers confirm from the record (title, description, data_processing, the file name, the series'
    overall_design) that the raw file is ONE library of ONE sample: one embryo (or one pool of one genotype) in one
    capture, not an aggregate of several captures or samples, and no other sample of any screened series names the
    same capture; a record that leaves this open means do not pin. They also check that the raw file is the sample's
    only file of its own (a raw file is taken only when it is the sample's only file of its own; a sample that lists
    any other file of its own never has its raw file fetched, and the plan table lists it as never fetched with that
    reason), that the raw file is single-cell or single-nucleus RNA (not ATAC, not a multiome's peaks), and that the
    knockout and its controls are read the same way (the code refuses a unit whose arms are read differently).

## The rules each step keeps

- **Held-out genes (the second line).** Before any file is chosen, `fetch.py plan` screens every sample of every listed series, plus
  every series-level field (title, summary, overall_design, type, file names, contact fields since stage 5, and the
  rest) and every field name. A held-out hit anywhere, a review term not cleared for that sample (or, for a series-level field, that
  series accession), or a record with no text refuses the whole dataset. Clearing one review term never passes a
  second term in the same text, and nothing clears a held-out term. The spelling families caught, the review terms
  and the known over-refusals are listed in `t3x_screen.py`'s docstring; each has a test string, and a property
  test inserts every held-out spelling into every field of a clean record. A Gata4 reporter knock-in counts as
  held out. A label that names neither gene nor drug (for example 'Mesp1-Cre; F/F' alone) cannot be caught by any
  text screen: that is why two people read every plan table and pin samples before any fetch. `fetch.py fetch` and
  `process.py run` screen each pinned record again before a file is requested or a count is read. Since stage 3 the
  screen over-refuses by design: no allowlist on either side of the beta short names, catalogue wording and the like
  as review terms, typo forms as review terms, and an EXTENDED form (a Ctnnb1 or Gata4 typo, beta + 'cat',
  'catenine', and since stage 4 '$\beta$-catenin' and 'β1-catenin') voids every review clearance of the sample it
  appears in. Stage 4 also deletes Unicode control and format characters (and U+034F) in an extra variant, folds
  Greek and Latin look-alike letters and non-ASCII digits, makes a text holding any control or format character
  (TAB excepted) a review term 'control-or-format character', and adds the spelling gaps the reviewers of 9c981e1
  listed as review terms (`t3x_screen.py` docstring). Stage 5 adds the aliases Evi, Notum, ICAT, Chibby / Cby1,
  dHand / eHand, Jup / plakoglobin / gamma-catenin and the Wnt targets Tbx6, Sp5 / Sp8, Tbxt / Brachyury / 'T/T',
  Mixl1 and Cdx2 as review terms. Per-cell label tables go through
  `t3x_screen.screen_cells`, but no loader uses it and a dataset that needs it is refused. Series-level files are
  never fetched unless the readers pin one under `pinned_series_files` (stage 6, "Series-level files"); a series'
  `_RAW.tar` never.
- **Phenocopy clause.** Wls (GSE165300), Hand2 (GSE126128) and Mef2c (GSE280543) have status `held_phenocopy`, and
  nothing fetches them. Each waits for one priced question to the captain. Tbx1 (GSE168768, GSE186934) is cleared
  on the engineer's reading that a Tbx1 loss is not a phenocopy of Gata4 loss; that reading goes into the same
  priced question (its `status_reason`).
- **Confidence floor.** The lineage transfer has never run on real 10x cells restricted to the panel. `process.py`
  writes each sample's confident and unconfident shares to `QC.json` and refuses a knockout gene's table (status
  `refused_low_confidence`, nothing written) when any of its samples has a confident share under
  `qc_floors.min_confident_share_per_sample` in `datasets.yaml` (0.5, fixed before any download).
- **Design.** Chimeras need a wild-type-chimera second difference. Without one they are `excluded_design`
  (Stat3 GSE259342, Tal1 E-MTAB-7325).
- **Stage window.** Every outside cell is capped at E9.5 inclusive, whole-embryo cells included. A range, an
  unparseable stage or a later stage fails closed. Outputs are T3 only, and a test asserts that no T1 or T2 builder
  reads this folder.
- **Scores.** No constant comes from a score. They are listed in `t3x_common.py`: confidence 0.6, half of an arm
  confident, call-share ratio 2, Spearman margin 0.05, 100 cells, 20 UMIs, pseudocount 10, n_h/(n_h+200) capped by
  m/(m+1), d = 1 or 0.5, clip [0.25, 4].
- **Repository.** Every write goes through `t3x_common.outside_repo`, which refuses this checkout and the main
  checkout. Downloads are never unpacked. A per-sample tar.gz and a pinned series archive are read member by member
  in memory, every member checked first (the per-sample route since stage 7).
- **Printing.** No step prints a gene-level value. Logs hold counts, statuses and sample titles of clean samples.

## How a table is computed

1. **Lineage groups.** `vec.harmonise.fit_transfer` is trained on the WT E8.75 reference's lineage groups directly
   (`obs['group']`), with the frozen settings: logreg, 50 PCs, C 1.0, clip 10, 5 depth draws, seed 0. The reference
   is refitted thinned to the query's depth when the query is shallower. So an outside call and a reference cell's
   own group mean the same thing in a builder (the challenger's change 2). A cell is confident when its
   `label_conf` is at least 0.6.
2. **Guards per stage and group.** A stage drops if under half of either arm's cells are confident. A (group, stage)
   drops if any arm has under 100 confident cells, or if the group's share of an arm moves more than 2x across arms.
   It also drops unless the control profile's best Spearman match among reference groups is its own group, ahead of
   the next by 0.05. Where author labels exist, under half agreement drops it too.
3. **Log ratio.** On summed UMIs per gene: `log((x_ko+10)/N_ko) - log((x_ctl+10)/N_ctl)`, clipped to log
   [0.25, 4]. A gene under 20 UMIs in any arm is NaN there (not measured): it holds no value in any stage, group or
   study mean, and a gene no group or study holds ends as exactly 1 in the finished table, not rescaled. Chimeras
   take [mut vs host] in knockout chimeras minus the same in wild-type chimeras.

   Two departures from the frozen factor: within a group the log ratio is taken on summed UMIs with a pseudocount
   of 10 (challenger change 3), not on props_type_means. Groups are pooled by the reference share-weighted mean of
   log values. Only the finishing step (rescale, clip, absent genes = 1) is gene_factor's.
4. **Shrinkage.** `L = k * ratio` with `k = min(n_h/(n_h+200), m/(m+1)) * m/(m+1) * d`. Here m is the smallest
   number of libraries of an arm (a series file used whole counts as one). d is 1 for germline or conditional
   knockouts against littermate or same-capture controls, and 0.5 otherwise; a unit with a control from a sister series
   takes 0.5 (stage 6). Stages are averaged with equal weight.
5. **Tables.** A study's genes are averaged with equal weight, then studies with equal weight. With 2 or more
   studies, a value is kept only if at least 2 studies share its sign and no opposite study exceeds the mean in
   magnitude; otherwise it is 0. The pooled form weights groups by the reference's group shares. Finishing then
   rescales so the reference's weighted panel total is unchanged, clips to [0.25, 4], and sets absent genes to
   exactly 1, last.

## Outputs (all under `~/v/vec-ext/t3/` on the DGX)

```
reference/refstats.npz                     WT E8.75 panel genes, groups, shares, mean proportions, WT sha256
meta/<ACC>/<series>_family.soft.gz         GEO metadata, sha256 in plan.json
meta/<ACC>/plan.json                       screen record, every sample's table row, pins verified (if any), units
meta/<ACC>/plan_table_<SERIES>.txt         the plan table the two readers read (verbatim fields, hash, verdict)
raw/<ACC>/<GSM>/<file>, raw/<ACC>/MANIFEST.json   pinned downloads with url, size, sha256; the pins and their sha256
processed/<ACC>/QC.json                    cells per stage and arm, call shares, guards, panel genes, status,
                                           plan_sha256, manifest_sha256, count_files_sha256
processed/<ACC>/<ACC>__<gene>.unit.json    log-space table for combine
processed/<ACC>/tables/<ACC>__<gene>/      finished per-dataset tables
tables/<name>/factors/<name>.csv           pooled form (gene, factor; reference gene order) + factors/factors.json
tables/<name>/<name>_per_group.csv         per-group form (group, gene, factor); a missing group = 1
tables/<name>/<name>_manifest.json         studies, genes, sign agreement, input sha256, 'checked' note
tables/COMBINE.json                        every single-cell table's files and sha256; its left-out and not-cleared datasets
tables/COMBINE_bulk.json                   stage 9: the same for the bulk family (only when bulk units are combined)
logs/<step>.log, logs/<step>.exit, logs/ALL.exit
```

`raw/` stays until the main session decides to delete it after the K21 reading. The units and tables do not need it.

## Running on the DGX

Nothing of this folder writes into the DGX's reference checkout `~/v/virtual-embryo-challenge`: no fetch into its
object store, no worktree of it, no output under it (`t3x_common.outside_repo` and `run_dgx.sh` both refuse it). The
code reaches the DGX through the main session's usual sync of `main` into the campaign checkout `~/v/vec-campaign-c`
(`docs/COMPUTE.md`, "The campaign checkout on the DGX"), after this branch is merged into `main`. The launcher refuses
unless `vec` is imported from that checkout. The WT E8.75 file is read, never written, through the checkout's
`data/` links (environment variable `WT` to point elsewhere). The bundle route of stage 2 is dropped.

The steps for the main session, in order (none was run by the engineer):

1. **Merge.** On this machine, merge `campaign/t3-ext-3` (stage 7, which holds stage 6) into `main` after its
   review; the merge changes `DATASETS_SHA256`, so every stored plan is stale and is made again by the launcher's
   re-plan (`run_dgx.sh ... all`) or a plan run, never by editing `plan.json` by hand. To use a
   sister series, first fill GSE137337's `control_series` placeholder (accession, screen_series, the three relation
   statements, the reader note) in a reviewed commit with the new `t3x_common.DATASETS_SHA256`; the relation
   statements and the note are written only after the two readers have read both series (step 4), so the first plan
   of the sister series runs with a declaration whose statements say what the hunt saw and are confirmed or corrected
   by the readers before any pin. If the readers change any word of control_series, commit it with the new
   DATASETS_SHA256 and run plan mode again before any pin-template; pins are made only from a plan made under the
   final declaration. The plan of a dataset with a sister series is refused unless the two series' records share an
   institute and a contact or publication value (section "A control from a sister series", stage 7).
2. **Sync.** Only when nothing runs in `~/v/vec-campaign-c`: the `docs/COMPUTE.md` rsync line with the target
   `PrakashDGX_H2:v/vec-campaign-c/`, then on the DGX `git -C ~/v/vec-campaign-c checkout --detach <main's commit>`
   and `git -C ~/v/vec-campaign-c status --porcelain --untracked-files=no` (empty). Check `df -h /home`.
3. **Plan run** (metadata only), detached:
   `cd ~/v/vec-campaign-c && mkdir -p ~/v/vec-ext/t3 && setsid nohup bash experiments/t3_ext/run_dgx.sh "$PWD" plan > ~/v/vec-ext/t3/launcher_plan.log 2>&1 < /dev/null &`
   and wait for `~/v/vec-ext/t3/logs/ALL.exit` with a background poller. `plan_<ACC>.exit` 0 means the screen passed;
   2 is a refusal whose log names it (a review term, or a held-out hit with the text withheld). GSE270969 is not in
   the launcher's default list: plan it with `run_dgx.sh "$PWD" plan "GSE270969"` (same detached form).
4. **Read.** Copy `~/v/vec-ext/t3/meta/` (metadata only, no counts) to this machine if wanted. The main session and a
   second, independent reader each read every `plan_table_<SERIES>.txt` of each dataset and apply the checks above.
   A review refusal is resolved only by `review_cleared` entries both readers agree (sample, term, reason); a
   held-out refusal is never resolved, and that dataset is dropped.
5. **Pin.** For each agreed arm, `fetch.py pin-template ACC --dest ~/v/vec-ext/t3 --arm ko|control --gene G
   --gsm GSM1 --gsm GSM2 ...` (every sample named; no network, it reads the stored plan) and paste its output under
   the dataset's `pinned:`. Then `fetch.py pin-template ACC --dest ~/v/vec-ext/t3 --series` and paste its output
   under a `pinned_series:` key beside `pinned:`. For samples whose cells are in a series file (stage 6): `pin-template
   ... --from-series-file` (it prints `files: []`), and per file `fetch.py series-file-template ACC --dest
   ~/v/vec-ext/t3 --series GSE --file F --gsm GSM ... --split none|barcode_suffix|folder [--key K=GSM ...]`, pasted
   under a `pinned_series_files:` key. Controls of the sister series are pinned with `--arm control` like any other.
   A sample of a `raw_matrices` dataset (GSE270969) is pinned with plain `pin-template` (its raw file is in
   `per-sample files allowed`) only after both readers have applied check 12 to it.
6. **Commit.** `fetch.py datasets-sha`, put the printed value into `t3x_common.DATASETS_SHA256`, run the t3_ext
   test files (tests/test_t3_ext_pipeline.py, _ratios.py, _screen.py, _stage_fields.py, _archives.py, _raw.py,
   _held.py, and since stage 9 _bulk.py and _bulk_ae.py),
   and commit `datasets.yaml` and `t3x_common.py` together on `main` (or on a branch merged into it), with the
   two readers named in the message. The `EXTERNAL_DATA.md` rows for the datasets that will be fetched go in before
   any table enters a model or a check.
7. **Sync again** as in step 2, so the DGX checkout holds the pinned `datasets.yaml`.
8. **Fetch, process, combine**, detached:
   `cd ~/v/vec-campaign-c && setsid nohup bash experiments/t3_ext/run_dgx.sh "$PWD" all "<pinned accessions>" > ~/v/vec-ext/t3/launcher_all.log 2>&1 < /dev/null &`.
   The launcher first runs the pinned gate (`fetch.py check-pins`; ALL.exit 6 if any accession has no pins or the
   YAML hash is not the pinned one), then the reference statistics, then per accession a fresh plan that checks each
   pinned series' and sample's LIVE metadata against its pinned hash (`replan_<ACC>`), and only if it passes the fetch of the
   pinned files and the processing of the pinned samples; combine last (`--dest`, so each QC is checked against its
   stored plan). Before each re-plan the launcher deletes that accession's `processed/<ACC>/`. Wait on
   `logs/ALL.exit` in the background. If fetch stops at the byte cap (`--max-gb 30`), raise `--max-gb` in a reviewed
   commit; do not split the run by hand.
9. **Afterwards**: read `processed/<ACC>/QC.json` (status, confidence_by_sample, the guard records; no gene-level
   value is printed anywhere) and `tables/COMBINE.json`. `raw/` stays until the main session decides to delete it.

Fix a selector, a pin or a clearance only in a reviewed commit, then plan again before any fetch of that dataset.
Never relax the screen.

## What stage 2 needs

- **A thin builder.** It is modelled on `experiments/t3_decomp/build_barred_exponent.py` and imports
  `build_submission` and `t3_decompose` unchanged, with their pinned sha256. Within its own process it replaces the
  barred s1 factor source by a named outside table. It refuses unless the frozen spec is
  `{'source': 'mab21l2', 'exponent': 1.0}`. The frozen `load_factor` takes the decision file's pinned
  `transfer_mab21l2_all` first, so no outside table goes through it under that name (the challenger's change 1). The
  builder checks the table's sha256 against its own pin and the table's `factors.json`. Its sidecar gets one more
  block naming the outside file, its sha256 and the form. A test must show that, given the Mab21l2 table, it
  reproduces the reference barred s1 bytes (3712c3f9...).
- **Forms, fixed before any score.** (1) The outside pooled table alone. (2) The product with the Mab21l2 factor:
  multiply, rescale, clip, then absent = 1. For the per-group form, build with phi = None, then apply the table per
  drawn cell's `obs['group']` through `apply_gene_factor` on a named stream `<stream>:outside`. For the per-group
  product, apply the Mab21l2 phi times the group table in one call per group.
- **Try K21 on the released Gata4 knockout.** Its row comes before any table is scored. Until a ruling clears Hand2
  and Mef2c, C1 is Tbx1 alone (GSE168768 plus GSE186934). It needs the permuted controls C1p, C3p and A1p, the
  leave-one-study-out reading, and the K18 check arm, without which the run is VOID. The row says what a pass does
  not show: Wnt transfer, cKO against germline, or blindness of g1 and g2.

## What stage 1 could not check

GSM titles beyond the scouts' records, so the selectors may refuse at plan time. File sizes and formats. Whether the
GSE205117 arms are chimeras (`design_refuse` checks this). Whether GSE137337's matrices are per sample (if they are
series-level only, its arms refuse). Whether the transfer from 10x cells on the panel genes to the MERFISH groups
holds up (the guards are there to drop groups where it does not). Per-study licences beyond GEO's default.

## What stage 3 could not check

Nothing ran against GEO: the plan tables, the hashes of real records and the launcher's `all` mode were tested on
synthetic SOFT text and a stub download only. Whether real GEO records parse to one stage per sample (pin-template
refuses a sample whose stage does not parse, so a sample may need a reviewed selector or parser change before it can
be pinned). Whether GEO edits a record between the reading and the fetch (the live hash check refuses then, by
design). How many review round trips the over-refusing screen costs on the five series (catalogue wording, 'Chi',
'2i', APC fluorochromes and G4 ES cells are all review terms). Cyrillic look-alikes are folded; other scripts'
look-alikes are not.

## What stage 4 could not check

Still nothing ran against GEO. That a real sample record lists every series it belongs to in `!Sample_series_id`
(one line per series) is GEO's documented SOFT format as I remember it, not checked on a real file: if a real record
lists only its own series, the membership check passes it, and the superseries is still screened through
`screen_series`; if it lists a superseries missing from `screen_series`, that sample is NOT PINNABLE until a reviewed
commit adds the superseries. Whether real family SOFT files carry sample table rows or other lines that would now be
kept under `_other_lines` (and so change the hash and enter the screen): I expect none for scRNA records. How a real
terminal shows bidirectional characters in the raw value line: the `escaped:` line below it is ASCII and shows them.
The Greek and Latin fold covers only the letters listed in `t3x_screen.py`; other scripts are not folded. Hand2, Mef2c
and Zfpm2 are review terms now, so a cardiac series that names them as markers needs a review clearance.

## Bulk whole-embryo units (stage 9, campaign/t3-ext-bulk)

Built from the reviewed design of 2026-11-21 as its challenger amended it, phases 1 and 2 only. Phase 3 (the reads
route, `quantify.py`) is not built: a dataset that would need reads has no count form here and cannot be declared.
Nothing of the single-cell path changes: a dataset without `modality` is single-cell and is checked, planned, pinned,
fetched, processed and combined exactly as before; its `pins_sha256`, its table names, its sign rule and its table
bytes are unchanged (tests in `tests/test_t3_ext_bulk.py`). `datasets.yaml` is unchanged in this stage (no bulk dataset
is registered and nothing is pinned): the main session writes the entries.

**The unit.** One knockout gene of one study at one pinned stage: its homozygous or null embryos against wild-type
embryos of the same series, line and stage. One library (one embryo, or one pool of one genotype) is one replicate.
Arms, genes and stages come from the pins only. Hets are never a knockout arm or a control; point, catalytic, ENU,
hypomorph and isoform alleles are refused by the readers (check B1); a double knockout is its own gene (not this
round). d is `design_factor` as now (1 only for littermates stated in the text and recorded as `control_relation:
littermate`). At least `BULK_MIN_LIBS` = 2 libraries per arm per stage, checked when the pins load and again after the
guards. Somite-only staging is refused (a sample whose only stage statement is a somite count parses to no E-day).

**The keys** (`t3x_common.check_bulk`): `modality: bulk`; `count_form` per_sample | series_table | arrayexpress;
`id_type` ensembl | symbol; `gene_columns` (one or two header names, the id column first; a second, symbol column is
recorded only); `count_column` (per_sample); `stage_from_series` (series_table: a plain series field whose one stage
every sample takes; a sample field stating another stage, a range or two values refuses the sample); `sex_field` and
`somite_field` (GEO: `characteristics_ch1:<key>` or a plain field); `source_name_pattern` and `withheld_lines`
(arrayexpress only). A bulk dataset takes none of the single-cell declarations (cell_labels, series_files,
control_series, raw_matrices, chemistry); a single-cell dataset takes none of the bulk keys.

**Pins** (`t3x_common.check_bulk_pins`): per_sample `{accession: GSM, arm, gene, stage, files: [one text table],
metadata_sha256}`; series_table the same with `column: "<exact header>"` in place of files; arrayexpress `{accession:
ERS..., source_name, arm, gene, stage, runs: [ERR...], metadata_sha256}`. The count table of a series_table or
arrayexpress dataset is pinned under `pinned_series_files` as `{series, name, header_sha256}` (exactly one entry).
`pins_sha256` takes it in only when present, as before.

### Intake A: per_sample

The single-cell own-folder rule, DENY and `sample_file_patterns` apply, and the file must be a plain text table
(`.txt`, `.tsv`, `.csv`, `.tab`, optionally gzipped). XLSX, RDS, h5 and other binary forms are listed as never fetched
("not a plain text table"). A sample listing more than one such table of its own offers none. The file's first line
must hold the declared gene columns (a headerless file, STAR's ReadsPerGene included, is refused); exactly one other
column is the count column, unless `count_column` names one, and then (or whenever the file has more than one other
column) that header must tie to this sample by the exact tie rule below; a second integral column (featureCounts'
Length) without count_column is refused. Two pinned samples whose files have one sha256, or identical count vectors,
refuse the dataset (G8). Round 2: the file's first line is screened before any value is parsed (every field and the
whole line through `t3x_screen`): a held-out hit fails the dataset's QC with the text withheld, and a review term
refuses unless a `review_cleared` entry names that sample. A single count column whose header ties by the exact rule
to ANOTHER sample (its accession, title, a description value or a characteristics value) is refused. count_column is
one header for the whole dataset, so it can tie to at most one sample: in practice a per_sample dataset takes only
files holding the gene column(s) and one count column; featureCounts (seven columns), HTSeq (headerless) and
multi-sample files are refused.

### Intake B: series_table, the header step and the exact tie rule

1. Only a file of a screened series' OWN folder is offered (`fetch.own_series_file`, a plain name, not denied, a plain
   text table, not NAMED as normalised or expected counts); its url is derived from the hashed
   `!Series_supplementary_file` value. A name matching `bulk.NORMALISED` (fpkm, rpkm, tpm, cpm, rsem, normalised,
   expected counts) is refused at the offer, again at the header step and again at the read (round 2), for this form
   and the ArrayExpress matrix alike.
2. **The header step** (refused while `t3x_common.GEO_BULK_OPEN` is False; `fetch.py header ACC --dest D --series GSE --name FILE`, or `run_dgx.sh ... header
   "ACC:GSE:FILE"`), after a passing plan and before any pin: the file is streamed (`fetch._open`), gzip is
   decompressed incrementally, reading stops at the first newline (at most `HEADER_MAX_BYTES` = 1 MiB of decompressed
   text, at most 16 MiB off the wire; a longer line refuses), the connection is closed, and only the line and its
   sha256 are written (`meta/<ACC>/header_<FILE>.txt` and `.json`). Every header field and the whole line go through
   `t3x_screen`: a held-out hit refuses with the line withheld (only the number of columns is written); a review term
   needs a `review_cleared` entry naming the series. The step prints the line, its `escaped:` form when it is not
   plain ASCII, and the tie table. The stream is read in blocks of 4,096 bytes (gzip decompressed 4,096 bytes at a
   time). Before the pins, at most one read block past the first newline passes through memory; it is never parsed or
   written.
3. **The tie rule** (mechanical, no tokens): a header ties to a sample of the screened series when it equals EXACTLY
   the sample's accession, its whole title value, a whole description value, or the value part of one characteristics
   entry (the text after the first `': '`). Every header is a declared gene column or ties to exactly one sample; no
   sample ties to two headers; each gene column appears once; no header appears twice; otherwise the whole table is
   refused. Readers confirm the printed tie table; they do not decide ties. Column order decides nothing.
4. `pin-template ACC --dest D --arm ko|control --gene G --gsm ... [--gsm-file F] [--column H]` prints each sample with
   the column the tie table ties to it (a `--column` given must equal it); `table-template ACC --dest D` prints the
   `pinned_series_files` entry with the header's sha256. The next plan carries the stored header records
   (`count_table_headers`), and `verify_pins` derives the tie table again from the stored line and the plan's rows.
5. `fetch` downloads the full file only after the pins are committed, into `raw/<ACC>/_series/`, its manifest record
   holding `series_metadata_sha256`, url, `pins_sha256` and `header_sha256`; a file whose first line does not hash to
   the pin is deleted and refused at once. `process` refuses it again unless its first line hashes to
   `header_sha256`, and reads it with `usecols` (the gene columns and the pinned columns only): the other columns'
   bytes pass through memory, never parsed (accepted by the design; a test holds a text column there). A
   tab-separated table is read without quote handling (`csv.QUOTE_NONE`), so a quote in an unpinned text column cannot
   merge rows.
6. The superseries rule: a series whose `relation` names a held accession refuses the plan (for a bulk dataset the
   own series' SOFT is fetched first, so the refusal comes before any other screened series' SOFT is fetched); a
   sample whose `series_id` lists a series outside `screen_series` is NOT PINNABLE, as before.
7. **Heterozygous wording** (round 2, both GEO forms): a sample whose own title, source name or a characteristics
   value matches `bulk.HET_WORDING` (`het` or `heterozygo...` as a word, `+/-`, `-/+`) is NOT PINNABLE in the plan
   table and refused by `pin-template` and by `verify_pins` in either arm. A false hit costs one sample; a reviewed
   exception can restore it later.

### The math (`bulk.bulk_log_ratio`, pure numpy)

Per library i: x_gi its counts (runs or lanes of one library summed), T_i its total over ALL gene rows of its file or
column (HTSeq `__*` and STAR `N_*` rows left out). y_gi = log((x_gi + 10) / T_i); L_g = mean over knockout libraries
of y_gi minus mean over control libraries, clipped to log [0.25, 4]. Because it is a mean of logs, every library has
equal weight and any per-library normalisation constant is a gene-independent offset that the finish's rescale
removes. Not measured (NaN; exactly 1 in the finished table): the median of x_gi over an arm's libraries under 20 in
either arm; a panel gene absent from any pinned library; the knocked-out gene's own row (bulk only; the single-cell
path is not changed); Xist, Tsix, the four Y genes of the sex call and every chrY panel gene of `panel_ids.tsv`,
always; the other chrX panel genes when the sex check (G3) finds an imbalance or cannot run. k = (m/(m+1))^2 * d
(`t3x_common.shrink_k_bulk`, equal to `shrink_k` with no cell term), m the smaller number of libraries left in an arm
after the guards; stages are averaged with equal weight. The unit file holds L under `__whole__` with `"modality":
"bulk"`, `present`, k and m. **Counts only**: every value finite, non-negative and integral within 1e-6; a table whose
library totals all lie within 1% of 1e6 (a CPM or TPM signature), a column or file named as FPKM, RPKM, TPM, CPM, RSEM,
normalised or expected counts, and a second integral column without count_column are refused (G2: the dataset's QC
fails).

### The id routes and the gene-id map (`panel_ids.tsv`)

`symbol`: process.to_panel's rule (exact symbol, then the committed ALIASES; duplicates summed). `ensembl`: the
version suffix is stripped; two table rows of one stripped id refuse the dataset; each panel symbol maps through
`panel_ids.tsv`; a panel gene whose id is NA is not measured; a dataset with both an id and a symbol column uses the
id and records the number of symbol disagreements. Under half the panel matched fails the QC, as before.

**The gene-id map** (the main session's decision: no annotation file is downloaded). `make_panel_ids.py --dest D
--accession ACC --gsm GSM --file FILE --refstats R [--extra SYMBOL ...]` reads the features file of an ALREADY PINNED
AND FETCHED 10x sample (`raw/<ACC>/<GSM>/<FILE>`, a pinned file of a cleared single-cell dataset whose manifest record
holds the current pins and the file's sha256; a loose features file or a per-sample archive, read through
`process.archive_members`) and the reference's gene list, and writes `experiments/t3_ext/panel_ids.tsv`: one row per
reference symbol in the reference's order (symbol, ensembl_gene_id or NA, chromosome, reason), NA where the source
gives the symbol no id or several, or the id names two panel symbols; header lines with the source file's name and
sha256, the reference list's sha256, and one `# extra` line per sex gene and per `--extra` symbol (the knocked-out
genes G9 reads). `t3x_common.PANEL_IDS_SHA256` pins the file; it was None until the main session committed the file in
a reviewed commit (set 2026-11-20: 'The committed map' below); while it is None or the file does not hash to it, every
dataset with `id_type: ensembl` is refused (the pinned gate, fetch, process). A 10x features file states no chromosome, so the chromosome column is NA: the
chrY rule then covers Xist, Tsix and the four Y genes only, and a unit whose sex check finds an imbalance (or cannot
run) is REFUSED (`refused_sex_chromosomes`) rather than read with unknown chrX genes.

**The committed map (2026-11-20).** `panel_ids.tsv` was made on the compute machine by `make_panel_ids.py` from the
features file of GSM8532456 (the pinned wild-type sample of GSE277839, fetched under its current pins), with the 25
pinned lines of E-MTAB-7322 as extras. The source was chosen on coverage alone, before any count was read: maps were
made from the features files of three pinned wild-type samples (GSM5169142 of GSE168768, GSM6202159 of GSE204906,
GSM8532456 of GSE277839); they give 494, 498 and 500 of the 500 panel genes one id, agree on every id that two of them
give, and give the same ids for the extras. Ccn5 and Cemip2 get an id from GSM8532456 alone (the two older features files
lack these symbols), so their ids have no second source. E-MTAB-7322 has no symbol column, so a symbol mismatch between
the reference and the map would not be counted as a disagreement; all three sources are 10x references. In all three, 23 of the 25 lines get an id, equal to the id the resource's
own sample table states for that line; Mir96 (a microRNA) and Zfyve20 (a symbol the 10x annotations no longer use) get
none, so the genotype guard (G9) refuses those two units. The pin is `t3x_common.PANEL_IDS_SHA256`;
`tests/test_t3_ext_panel_ids.py` holds the file to it and holds the extras to the ids the resource states. An
independent reviewer re-derived every figure of this paragraph and found the commit ready (decision log 2026-11-20
(ggg)).

### The guards (fail closed; design numbers fixed before any count; verdicts to QC.json; no gene-level print)

- **G1 library size**: T_i >= 1,000,000 and a present-panel total >= 10,000, else the library is dropped; the largest
  to smallest library total of the unit is recorded.
- **G2 counts only**: above; a failure fails the dataset's QC.
- **G3 sex**: per library a call from Xist and Ddx3y, Eif2s3y, Kdm5d, Uty (CPM over T_i): male when the Y sum >= 10
  and Xist < 10, female the reverse, else unclear; a clear call contradicting the recorded sex drops the library; a
  male share differing between the arms by more than 0.5, or a sex gene absent ('not possible'), sets the chrX panel
  genes not measured (or refuses the unit when the map carries no chromosome).
- **G4 outlier and grouping**: in an arm of 3 or more libraries a library is dropped when its Spearman with the
  leave-one-out arm mean (log CPM, present panel genes) has a robust z under -3.5 AND lies more than 0.05 under the
  arm's median; the unit is refused when the mean Spearman within the controls is under the mean across the arms
  minus 0.05.
- **G5 single-lineage check** (bulk; changed 2026-11-21, before any bulk table was read): pb = sum over groups of
  shares[g] * mu[g] (`t3x_tables.whole_den`). A control library whose Spearman with pb is not strictly greater than its
  Spearman with every single reference group (margin > 0, `t3x_common.BULK_WHOLE_MARGIN`) is dropped. A control arm mean
  that fails refuses the unit. Knockout libraries are recorded only. The margin, the threshold and the best group are
  recorded per library and for the control arm mean. G5 catches a library dominated by one or two lineages (a dissected
  single tissue). It does not show that a library is a whole embryo. In the challenger's synthetic check (random group
  profiles; not kept in the repository), partial embryos (head-off, trunk-heavy) and embryos mixed with a tissue the
  reference lacks scored about like whole embryos at any margin; no test here shows it. Those
  cases rest on the readers' record check (B6), the composition guard (G12) and the genotype guard (G9). G5 first used
  SPEARMAN_MARGIN (0.05), a number borrowed from the single-cell group check. On E-MTAB-7322 (217 single whole E9.5
  embryos) it cut one tight unimodal margin distribution at its median (controls: median 0.051, range 0.038 to 0.069)
  and dropped 46 of 111 controls. Across all
  libraries, the share failing was highest at the somite extremes (0.70 at 12 or fewer somites, 0.60 at 25 to 27, 0.78
  at 28 or more). The single-cell route keeps SPEARMAN_MARGIN = 0.05
  unchanged.
- **G6 somites**, where recorded: non-overlapping arm ranges refuse; the median difference, the ranges and the
  overlap count are recorded; nothing is removed.
- **G7 stage and design**: the E9.5 cap, the range refusal, same-stage arms (the pins).
- **G8 duplicates**: identical count vectors or file sha256 among pinned libraries refuse the dataset; one count file
  sha256 under two bulk datasets refuses combine.
- **G9 genotype in the counts**: the knocked-out gene's own clipped log ratio >= 0 (or the gene not found in the
  counts) refuses the unit; only the verdict is written.
- **G10**: every single-cell refusal is kept (the panel match 0.5, pins, plan and manifest binding, the screen again
  before any count, symlinks, stale output deleted first, NEVER_FETCH and the held set, outside_repo, no gene-level
  print).
- **G12 composition**: the NNLS fit of each arm's mean panel proportion profile onto the reference groups' mu[g],
  normalised to shares; a group with a control share >= 0.05 whose knockout-to-control share ratio is above 2 or
  under 0.5 refuses the unit. The shares are recorded (group-level).

### The bulk table family (combine)

`load_units` accepts the key `__whole__` only from a bulk unit (and a bulk unit holds no other key); a unit's
modality must be its dataset's; a study spanning modalities (or tiers) is refused. Single-cell units build the
single-cell family exactly as before (names, rule, bytes). Bulk units build `bulk_all`, `bulk_study_<s>` and, when
`bulk_all` has 2 or more studies, `bulk_all__loso_<s>` for every study: pooled form only (no per-group file; the
manifest says `per_group: none (bulk)`), finished against the whole-embryo pseudo-bulk. No `both` family and no bulk
tier tables this round. Inside a study equal weight per gene, across studies equal weight per study. **Sign rule**:
a bulk table of 2 or fewer studies keeps today's rule; of 3 or more, the mean is kept iff the studies sharing its
non-zero sign number at least max(2, ceil(2h/3)), h the studies holding the gene, else 0 (the design's median clause
is struck: it is vacuous once two thirds share the mean's sign). Single-cell tables keep today's rule for every number
of studies. Manifests add `modality_by_study`, `sign_rule`, `k_by_unit`, `libraries_by_unit`, `guards_by_unit`, the
within-study line agreement per gene (a study of several genes, for reading only) and `checked: by the practice read
named in its TRIES row only`. Round 2: each family's manifests carry a combine record limited to its own accessions
(inputs, left-out units, not-cleared datasets): the single-cell manifests and `COMBINE.json` are byte-identical with
and without bulk units in the same run (apart from each manifest's `created` time), and the bulk family's record goes
to `COMBINE_bulk.json`.

## The many-line resource (stage 9, phase 2: E-MTAB-7322)

`count_form: arrayexpress`, for E-MTAB-7322 only (`ae.AE_ACCESSIONS`; study key MGP-DMDD; E-ERAD-* and E-MTAB-5224
are not taken). Every url is derived from the accession on one host (`fetch.AE_HOSTS`, https, no port, no query):

```
IDF, SDRF      https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/322/E-MTAB-7322/Files/E-MTAB-7322.idf.txt (.sdrf.txt)
count matrix   https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-7322/E-MTAB-7322-raw-counts.tsv.undecorated
methods text   https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-7322/E-MTAB-7322-analysis-methods.tsv
```

1. **plan** downloads the IDF, the SDRF and the methods text only. Text that is not UTF-8, a line or paragraph
   separator inside a line, an IDF key or SDRF header outside the closed grammar (`ae.IDF_KEY`, `ae.SDRF_HEADER`), a
   row of another width, or a required column (Source Name, ENA_SAMPLE, mutant line, genotype, developmental stage,
   somite stage, sex, individual, allele, ENA_RUN) not appearing exactly once refuses the resource with its text
   withheld.
2. **withheld_lines**, applied BEFORE the screen, the views, the hashes and the tables: every SDRF row whose
   `Characteristics[mutant line]` equals a withheld value exactly is dropped, and so is every row of a Source Name
   one of whose rows shares a Source Name, ENA sample, individual or run (`Comment[ENA_RUN]`) with a withheld row
   (to a fixed point; an empty value links nothing). Only `<gene>: n rows withheld (reason)`, the number of linked
   rows and the number of withheld runs are printed. The withheld runs are kept in the plan as sha256 digests only;
   `pin-template` and `verify_pins` refuse a pin or a record that lists one (round 2). A line
   value equal to a withheld value after case folding and removing spaces, `-`, `_` and `.` (a near-match) refuses the
   resource. A held-out term or a held accession anywhere in what remains withholds the whole resource, as now.
3. **Records**: a sample is one Source Name; its fields map `<3-digit column number>|<header>` to its rows' values in
   file order (nothing stripped), hashed with `canonical_metadata(Source Name, fields)`. The series block is the IDF
   (`idf|<key>`), the SDRF header row and the methods text (shown verbatim, hashed, pinned as `pinned_series`). Each
   Source Name's runs (ENA_RUN, paired rows counted once) are summed into one library.
4. **The views**: per column the distinct values with their number of Source Names (past 60 values, their number and
   the pattern view, digits as `#`), and per line the genotype-by-stage counts, the sexes, the somite range and the
   pinnable E9.5 Source Names.
5. **The mechanical comparisons** (a) to (f), each printing its number of disagreements (module docstring of `ae.py`;
   (g), the pin hash against the plan hash, runs in `verify_pins`); a disagreeing Source Name is NOT PINNABLE. Only
   `embryonic day 9.5` is pinnable (days 10 and 10.5 fall under the stage cap; 8, 8.5 and 9 are not pinnable); hets and
   the line `baseline` are never pinned.
6. **The committed rule fixes the line set**: every line not withheld (not `baseline`) with at least 2 homozygous and
   2 wild-type pinnable E9.5 Source Names, comparison (h) included. Round 2: the header ties are part of the rule.
   Once the header step has passed and the plan is made again, (h) marks every Source Name with a run missing from
   the matrix header (`(h) a run missing from the matrix header`) or tied through a header that is not a run (`(h)
   tied through a non-run header`) as NOT PINNABLE, and the line set is computed over the remaining rows; a line that
   drops under 2 in an arm leaves the set. One function (`ae.committed`) gives the set to the plan, `pin-template`
   and `verify_pins`. Before the header step the plan table says `THE LINE SET IS PROVISIONAL`, and `pin-template`
   refuses until the plan carries the passing header record. `pin-template ACC --dest D --arm ko|control --gene LINE
   --source ... [--source-file F]` refuses unless the names are exactly the rule's eligible Source Names of that line
   and arm, and `verify_pins` refuses pins that are not exactly the rule's set (no hand-picking).
7. **The header step** (`fetch.py header E-MTAB-7322 --dest D`): the matrix's first line, streamed and closed at the
   first newline, screened, tied by exact equality with an ERR run (or ENA sample, or Source Name) of a KEPT record. A
   header tying to two records, a doubled header or a gene column not appearing once refuses the table; a header tying
   to no kept record (the 369 unannotated runs, the withheld lines' runs) is never read and only counted; a record with
   a run missing from the header, or tied through a non-run header, is not pinnable (comparison (h), item 6). The line
   is kept in the record for the re-derivation and never printed. The matrix name is checked against
   `bulk.NORMALISED` at the header step and at the read.
8. **fetch** downloads the matrix (about 380 MB) into `raw/<ACC>/_series/` only after the pins; **process** reads it
   with usecols = `Gene ID` and the pinned runs; the matrix's gene ids are Ensembl ids, so the resource needs the
   pinned gene-id map.

An example entry (the main session writes the real one; `source_name_pattern` follows the form the main session read
in the resource's public metadata, no count value: `<line>_<genotype><n>`, the genotype one of `wt`, `het`, `hom`, for
example `Pigf_wt1` or `4933434E20Rik_het1`, the same string as `Characteristics[individual]`; the withheld lines are
the two readers' result):

```yaml
  - accession: E-MTAB-7322
    status: cleared
    status_reason: "..."
    tier: unrelated
    study: MGP-DMDD
    kind: "DMDD whole-embryo bulk RNA-seq of knockout lines and wild-type littermates, E9.5"
    design: germline
    control_relation: littermate          # only if the IDF states littermates in words (check B4); else same_series
    repository: arrayexpress
    modality: bulk
    count_form: arrayexpress
    id_type: ensembl
    gene_columns: ["Gene ID"]
    withheld_lines:                       # applied before the screen, the views, the hashes and the tables
      - {line: Trps1, reason: held}       # a GATA-type zinc-finger factor ('any GATA factor')
      - {line: Anks6, reason: pending_phenocopy}
      - {line: Atxn10, reason: pending_phenocopy}
      - {line: B9d2, reason: pending_phenocopy}
      - {line: Cc2d2a, reason: pending_phenocopy}
      - {line: Ift140, reason: pending_phenocopy}
      - {line: Kif3b, reason: pending_phenocopy}
      - {line: Kifap3, reason: pending_phenocopy}
      - {line: Rpgrip1l, reason: pending_phenocopy}
      - {line: Mad2l2, reason: pending_phenocopy}
      - {line: Tmem98, reason: pending_phenocopy}
      - {line: Nek9, reason: pending_phenocopy}
      - {line: Mtor, reason: pending_phenocopy}
      - {line: Map2k7, reason: pending_phenocopy}
      - {line: Traf6, reason: pending_phenocopy}
      - {line: Traf2, reason: pending_phenocopy}
      - {line: Rala, reason: pending_phenocopy}
      - {line: Pigf, reason: pending_phenocopy}
      - {line: Pigl, reason: pending_phenocopy}
      - {line: Pgap2, reason: pending_phenocopy}
      - {line: Hira, reason: pending_phenocopy}
      - {line: Smc3, reason: pending_phenocopy}
      - {line: 1700007K13Rik, reason: pending_phenocopy}
      - {line: 4933434E20Rik, reason: pending_phenocopy}
      - {line: Hmgxb3, reason: pending_phenocopy}
      - {line: Zkscan14, reason: pending_phenocopy}
      - {line: Zkscan17, reason: pending_phenocopy}
      - {line: Dhx35, reason: pending_phenocopy}
      - {line: Gpr107, reason: pending_phenocopy}
    source_name_pattern: '(?:M\d+_)?(?P<gene>[A-Za-z0-9]+)_(?P<genotype>hom|het|wt)\d+'   # [M<number>_]<line>_<wt|het|hom><n>
    screen_series: [E-MTAB-7322]
    stages_used: [E9.5]
    licence: "..."
    why: "..."
    risks: "..."
    pinned: []
```

## What the two readers check for a bulk dataset (the design's final checklist)

Checks B1 to B10 replace check 7 for a bulk dataset; checks 1 to 6, 8 and 9 stay as written. Any doubt means nothing in
the dataset is pinned and the question goes to the main session.

- **B1 Genotype and arm.** The sample's own text states homozygous, null or -/- for the named gene (knockout arm), or
  wild type (control). A bare KO, Mut or CTRL label needs another field of the same record naming the gene. Hets,
  point, ENU, catalytic, hypomorph and isoform alleles are not pinned.
- **B2 Held and phenocopy.** Check 2 applies to every series block, sample, IDF line, the header line and the line
  names. Lines in withheld_lines (Trps1: held; the 28 others of the example entry: pending phenocopy) never appear in
  the views; confirm only the printed `<gene>: n rows withheld` lines. Trps1 is not a screen term; a Trps1 mention in
  any GEO bulk dataset is caught by readers only and means do not pin. Check the accession
  character by character: E-MTAB-7322, not E-MTAB-7325. If any series text links a gene to Wnt, beta-catenin or a GATA
  factor, do not pin, and send the question to the main session. Sox7 (SoxF) and Cables2: the screen and this check
  decide. A Sox17 mention is a reader's review item (check 2), not a screen term: the code does not catch it.
- **B3 Library design.** One embryo, or one single-genotype pool, per library. Read the IDF protocols or the
  extract_protocol, data_processing and overall_design in full: pooling of mixed genotypes, hashing or demultiplexing
  means do not pin. Where the text does not say embryo or pool, note it. GEO has no individual field; two GSMs of one
  embryo are a reader check (B3).
- **B4 Same-study control.** Same series and line, same stage. Littermates stated in words (quote them) sets d = 1,
  else 0.5. Not 'previously published', not a baseline staging series, not another strain or line.
- **B5 Stage.** One E-day at or before E9.5, from stage_fields or stage_from_series. Samples past E9.5 in the same
  series (E10, E10.5, adult, foetal liver) stay unpinned; ES and EB arms stay unpinned. Any stage value with leading or
  trailing space is a disagreement.
- **B6 Assay.** Bulk RNA-seq of a whole embryo, with yolk sac and placenta not included or stated. Not ribo-seq,
  bisulfite, ATAC, microarray, spatial, dissected tissue or ES/EB. The file must not be described as FPKM, TPM, RSEM
  or normalised. Values are checked by code, never by eye.
  Where the text says whole embryo or whole organism and is silent on yolk sac and placenta, write it as a note: silence
  is not a refusal (the main session's ruling of 2026-11-20 for E-MTAB-7322, docs/PLAN.md decision log 2026-11-20
  (fff)); a text that names a dissected part, or different handling of the two arms, is.
- **B7 Files and ties.** per_sample: the file is in the sample's own folder and is exactly one text table.
  series_table: the series file url is in the series' own folder; the header was screened; the printed tie table
  shows every non-gene header tied to exactly one sample by exact equality, with 0 problems; gene_columns are the right
  columns; header_sha256 equals the printed one. A table whose headers do not tie is not pinned; do not decide a tie
  by eye. Superseries: if a sample lists a superseries, confirm its sample list holds no held accession; if unsure, do
  not take the dataset.
- **B8 Mechanical comparisons (ArrayExpress).** (a) to (h), including one Source Name per individual and per ENA
  sample, printed 0 disagreements for the Source Names pinned, or the disagreeing ones stay unpinned. The line set
  equals the committed rule's output: non-withheld lines with >= 2 hom and >= 2 WT at E9.5 after (h). The plan table
  read for pinning says the line set is final (after the header step), not provisional. No hand-picking.
- **B9 Hidden characters, hashes and superseries rows**: as in checks 8 and 9.
- **B10 Sex and somite fields.** Read their distinct values and note an arm that is all one sex. Nothing is decided on
  this: G3 and G6 handle it.

Each reader writes per dataset: the units (gene, knockout and control library counts per stage); d, the relation and
its source words; count_form, id_type, count_column and gene_columns; any doubt. Never read a count value, a QC gene
value or a table value. Nothing is pinned or added after a practice read without a new TRIES row.

## Stage 9: the steps for the main session, in order (none was run by the engineer)

1. **Merge** `campaign/t3-ext-bulk-3` (it holds `campaign/t3-ext-bulk-2`, the round 3 changes, and main) into `main`
   after its review. The GEO bulk forms stay closed (`t3x_common.GEO_BULK_OPEN = False`); only the arrayexpress form
   is usable. This branch changes nothing in `datasets.yaml`, so `DATASETS_SHA256` stays main's; stored plans stay
   valid for single-cell datasets. Sync the DGX campaign checkout.
2. **The gene-id map** (E-MTAB-7322 and any `id_type: ensembl` dataset) comes AFTER the pins (step 9) and before step
   10, because its `--extra` list is the pinned lines: make the reference once (`process.py reference --wt WT --out
   D/reference/refstats.npz`, as `run_dgx.sh all` does), then run `make_panel_ids.py --dest D --accession <a cleared
   single-cell dataset already fetched> --gsm GSM --file <its pinned features file or per-sample archive> --refstats
   D/reference/refstats.npz --extra-file <the pinned lines, one a line>`; commit `panel_ids.tsv` with its sha256 in
   `t3x_common.PANEL_IDS_SHA256` in one reviewed commit. A pinned line missing from the extras is refused at G9.
3. **Dataset entries**: one reviewed commit per batch, each entry with `modality: bulk` and its keys, the dataset-file
   hash updated in the same commit (`fetch.py datasets-sha`). For E-MTAB-7322 the example above, withheld lines
   first. Every stored plan is then made again (the dataset-file hash changed).
4. **Plan** (metadata only): `run_dgx.sh "$PWD" plan "ACC ..."`. For the resource the plan table says the line set is
   provisional.
5. **The header step** (series_table and arrayexpress only): `run_dgx.sh "$PWD" header "ACC:GSE:FILE"`; for
   E-MTAB-7322 `run_dgx.sh "$PWD" header "E-MTAB-7322"`. Read the printed status, tie counts and (h) lists.
6. **Plan again**, so the plan carries the header records; for the resource the plan table now says the line set is
   final and lists (h).
7. **The readers**: two readers read the series block, the column and per-line views, the comparison counts, the
   line set and the full record of every Source Name in the line set (for the resource); every plan table, header
   line and tie table (for a GEO dataset). The held-out screen has run mechanically over every record; a reader who
   sees a problem in a view opens the records it names. Any doubt means nothing is pinned. They apply checks 1 to 6,
   8, 9 and B1 to B10.
8. **Register rows**: the `EXTERNAL_DATA.md` row of each dataset to be pinned, written by the main session (this
   branch does not edit `EXTERNAL_DATA.md`) before any table enters a model or a check.
9. **Pins**: `fetch.py pin-template ACC --dest D --series`, `pin-template` per arm (`--gsm`/`--gsm-file`, or
   `--source`/`--source-file` with exactly the final line set for the resource), `table-template` for a count table;
   paste, `fetch.py datasets-sha`, commit with the two readers named. Then plan again (the plan with pins checks them
   live).
10. **Fetch, process, combine**: `run_dgx.sh "$PWD" all "<accessions>"` (the pinned gate refuses an Ensembl-route
    dataset while the id map is not pinned). Raise `--max-gb` only in a reviewed commit if the matrix exceeds it.
11. **Afterwards**: read each `processed/<ACC>/QC.json` (status, the guard verdicts per unit and stage, the library
    records, each library's count column header under `column`, which check B2 applies to; no gene-level value),
    `tables/COMBINE.json` and `tables/COMBINE_bulk.json`, then the practice read named
    in its TRIES row (its forms fixed first: bulk_all alone; bulk_all times the frozen factor; the comparison with the
    single-cell all table; and bulk_all__loso_MGP-DMDD if the resource is in).

## What stage 9 could not check

Nothing ran against GEO, BioStudies or Expression Atlas, and nothing ran on real data: the tables, headers, IDF,
SDRF, methods text and matrix of the tests are synthetic. The resource's Source Name form comes from the main
session's read of its public metadata (round 2); the engineer did not see it. If a real name does not match
`source_name_pattern`, comparisons (b) and (c) make that Source Name NOT PINNABLE until a reviewed change. Whether wild-type littermates'
`Comment[allele]` names the line's allele was not seen (comparison (c) needs it). Whether `Characteristics[individual]`
values are unique across lines was not seen (comparison (e) refuses a repeated value). The SDRF's paired rows were
assumed to repeat ENA_RUN with one FASTQ_URI each. The matrix's real header was not read; whether every annotated
run is a column is checked by the header step. A 10x features file states no chromosome, so chrX and chrY panel genes
other than the named sex genes are unknown: an imbalanced unit is refused. Trps1 was not added to the screen (see the
hand-back). How the guards behave on real whole embryos (G5 against an E8.75 pseudo-bulk, G12 on a bulk profile) is
unknown. Note (2026-11-21): No library known not to be a whole embryo has been measured, so there is no positive
control. Whether an E7.5 or later-stage embryo, or a whole embryo carrying yolk sac, fails G5 at margin > 0 is unknown.
The earlier expectation that E7.5 studies fail G5 was made for 0.05 and is withdrawn.

## Stage 9 round 3 (campaign/t3-ext-bulk-3, by the main session): the real tables, and two closures

After the second verification (not ready: two blockers, both on the GEO forms) the main session ran the planner, in a
scratch copy, on the resource's own description and sample tables (metadata only; no count read) and changed four
things, each with a test in `tests/test_t3_ext_bulk_ae.py`:

1. **The description file's keys may hold a hyphen.** The real file has the key `MAGE-TAB Version`; the key grammar
   refused it and with it the whole plan. Any other punctuation in a key still withholds.
2. **An allele is written `line<allele>` or `line^allele`.** The real sample table uses the caret; with only `<` read,
   569 of 680 kept records failed comparison (c).
3. **The named resource must withhold its held line.** `check_bulk` refuses E-MTAB-7322 unless `withheld_lines` holds
   `{line: Trps1, reason: held}` (`t3x_common.REQUIRED_HELD_LINES`).
4. **The two GEO bulk forms are closed.** `t3x_common.GEO_BULK_OPEN` is False: a dataset with `count_form: per_sample`
   or `series_table` is refused at load until a reviewed change closes the two open blockers (heterozygous wording
   written with a Unicode minus or as an allele pair such as KO/+; a per-sample count column that names another sample
   by a token) and sets it to True. The arrayexpress form passed its end-to-end verification and is not affected.

Real forms the main session saw in the sample table: Source Names are `M<number>_<line>_<wt|het|hom><n>` for most
embryos and `<line>_<wt|het|hom><n>` for some (the example entry's pattern takes both); a few lines spell the line
differently in the Source Name and in the line column, which leaves them not pinnable; `Term Source REF`, `Term
Accession Number` and `Comment[technical replicate group]` carry stray whitespace, and no comparison reads them (they
are hashed verbatim); `Comment[allele]` also carries some, and comparison (c) reads only its prefix. Two README sentences are corrected: before
the pins the code reads at most one 4,096-byte block past the first newline, and the HTTP and TLS layers below it may
buffer a few more kilobytes, none of it parsed or written; after the header step the plan carries the matrix header
line verbatim for the re-derivation, withheld runs' accessions included, never printed or tied.

After the verification of round 3 (two small blockers, both closed here): comparison (c) takes the allele prefix from
every row of a Source Name (whitespace after the allele part differing between rows no longer fails it); a
description key may hold a hyphen only between two parts (`MAGE-TAB Version` yes, `Key--` no); a row is withheld also
when its Source Name's gene token or its allele prefix names a withheld line, whatever its line column says; the
many-line tests run with the GEO forms closed; and the operator steps name this branch and put the gene-id map after
the pins. **Expectation for the reader of the QC:** with no chromosome column in the gene-id map, a unit whose arms
differ in sex mix by more than half, or whose sex cannot be called, is refused (`refused_sex_chromosomes`); with two
to four embryos an arm this takes a noticeable share of units and is not a fault.

## Stage 9 round 2 (campaign/t3-ext-bulk-2): what changed after the two reviews of 120af2f

Each blocker has a test that replays the reviewer's own probe (`tests/test_t3_ext_bulk.py`,
`tests/test_t3_ext_bulk_ae.py`, names starting `test_r2_`). Every refusal of the earlier rounds stays.

1. **A withheld line's runs are withheld.** Every run of a withheld row joins the linked keys; a kept row that carries
   one (or shares a Source Name, ENA sample or individual) is withheld with every row of its Source Name, to a fixed
   point, and counted; empty values link nothing. The withheld runs stay in the plan as sha256 digests only, and
   `pin-template` and `verify_pins` refuse a pin or a record that lists one.
2. **Per-sample headers are screened** before any value is parsed (`bulk.screen_sample_header`): a held-out hit fails
   the dataset with the text withheld; a review term refuses unless a `review_cleared` entry names that sample.
3. **A single count column that ties by the exact rule to another sample** (its accession, title, description or a
   characteristics value) is refused.
4. **A count table named as normalised or expected counts** (`bulk.NORMALISED`, now with `expected counts`) is
   refused at the offer, at the header step and at the read, for the series form and the ArrayExpress matrix.
5. **Heterozygous wording** in a GEO sample's own title, source name or characteristics (`bulk.HET_WORDING`) makes
   the sample NOT PINNABLE and is refused by `pin-template` and by `verify_pins` in either arm. The pattern is
   stricter than the reviewer's (`het` is matched between non-letters, so `HET_rep1` and `Pigf_het1` hit).
6. **The header ties are part of the committed rule** (`ae.committed`, comparison (h)): see "The many-line
   resource", item 6.

Lesser points applied: the pre-pin read block is 4,096 bytes; each table family gets a combine record limited to its
own accessions (`COMBINE_bulk.json` for the bulk family); a tab-separated table is read with `csv.QUOTE_NONE`; the
readers' step for the resource reads the views and the full records of the line set, not every record; Trps1 and the
GEO individual sentences in B2 and B3; the count_column sentence in Intake A; `process.py run` prints a bulk unit's
status and reason; a four-field matrix line is refused on the filtered single-cell route too (`process._three_fields`,
before scipy's `mmread`). The example entry for E-MTAB-7322 carries the Source Name form and the 29 withheld lines.

Not changed: the screen (`t3x_screen.py`). Trps1, Sox17 and the other reader-only terms of check 2 are not screen
terms, so a per-sample count column headed by one of them is read; its header is recorded in QC.json under
`libraries` (`column`), and the readers read those headers after `process` (step 11 above).

## Stage 8 (campaign/t3-ext-4): what changed after the two reviews of stage 7

Each blocker of the first review of ec86a1b has a test that replays the reviewer's own reproduction
(`tests/test_t3_ext_raw.py`, `tests/test_t3_ext_archives.py`):

- **A raw file only when it is the sample's only file of its own.** A sample listing a raw archive beside any other
  file of its own folder (a filtered h5, a filtered archive, a file the patterns do not offer) never has the raw file
  fetched: it is listed as never fetched with the other names. Before, a filtered file the patterns did not offer left
  the raw file offered.
- **Both arms of a unit are read the same way.** `verify_pins` (`fetch.check_unit_forms`) refuses a unit whose
  knockout and control samples are not all raw-called or all filtered, and `process` refuses the same as a second
  guard ("an unfiltered and a filtered library in one unit").
- **A sister series needs a publication or a repository link.** Two different PubMed ids refuse; otherwise a shared
  PubMed id or a super- or sub-series link is required, and the plan prints which link matched.

The lesser points: a matrix line of four fields is refused (it was read with its first field as an index); more than
1,000 comment lines before the size line refuse; the team's held list is in the code (`HELD_EXTRA`) and GSE205117 is
held; per-sample and per-arm cell-call records in QC (records only); README lines on the barcode cap, archive folder
names and GSE270969's plan command.

What stage 8 could not check: nothing ran against GEO or on real data. Real GEO `relation` values were not seen: the
code takes only the exact forms `SuperSeries of: GSEnnn` and `SubSeries of: GSEnnn`, and a record that writes them
otherwise gives no link (the plan then refuses unless a PubMed id is shared). Whether GSE270969's samples list any file
besides the raw archive is unknown: if they do, none of its samples can be pinned. The GSE137337 sister path may now
be blocked if the sister series has no PubMed id and no repository link to GSE137337.

## Stage 7 (campaign/t3-ext-3): what changed after the two reviews of stage 6

Commit 1 of the branch, each with a test that replays the reviewer's reproduction (`tests/test_t3_ext_archives.py`):

- **Split none counts libraries.** A file used as one unit is refused unless it holds exactly as many libraries
  (folders with a triple, or the one loose or flat triple) as its entry names samples. An archive with folders
  `KO_1a`, `WT_littermate` and `unlisted_3` pinned to one knockout sample is refused and nothing is written. So a
  flat or loose file that aggregates two samples' libraries (one triple) cannot be used whole any more: it needs a
  barcode suffix rule from the records. The 60,000-barcode cap applies to each library.
- **A sister series is checked by code** (section "A control from a sister series"): not registered as or screened
  by another dataset, a shared institute and a shared contact or publication value, no superseries of the own series
  under the sister.
- **A one-library file may control several knockouts.** Split none needs one arm and one stage, each sample pinned
  for the same genes, not one gene.
- **The per-sample archive route** (stage 1) reads every member through the same checks as a series archive, and a
  per-sample archive holds exactly one triple in one folder and nothing else; a gzipped member is decompressed in
  chunks and refused past 8 GiB.

Commit 2: the intake for per-sample unfiltered matrices (section "Unfiltered per-sample matrices"), declared on
GSE270969 only, with no pins, and its tests (`tests/test_t3_ext_raw.py`).

What stage 7 could not check: nothing ran against GEO or on real data. GSE270969's real file names, archive layout
(a `raw_feature_bc_matrix/` folder or loose members, gzipped or not), feature types and sizes were not seen: a name
the declared pattern misses is listed as never fetched, and an archive that does not fit the member rules is refused
after the download. How the fixed cell call behaves on single-nucleus libraries of whole E9.5 embryos (how many
nuclei it calls, how much ambient RNA the called nuclei carry) is unknown until the QC records exist; the readers of
those records do not change the rule. The sister-series match was tested on synthetic records only: real GEO
contact fields may differ in spelling between two series of one lab, and the plan then refuses (by design) until a
reviewed change.

## What stage 6 could not check

Nothing ran against GEO or on real data: the archives, loose triples and sister series of the tests are synthetic.
The series folder rule assumes GEO's series links have the form
`ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137337/suppl/GSE137337_<name>`, as I remember it; a real record
whose links differ offers no file. GSE137337's archive names other than the brief's example, the loose files' names,
what the archives hold (one library or several; folders or barcode suffixes; whether each holds a 10x triple at all;
member names with other characters), and whether the records state a split rule were not seen: a file whose contents
do not fit its pin is refused by process after the download, and a member name outside `[A-Za-z0-9._+-]` (a space,
for example) refuses the whole file until a reviewed change. A per-sample tar.gz (stage 1 route) was still read
without the stage 6 member checks; stage 7 closes that. The sister series was not named (the placeholder
is inert), so no real record of it was screened.

## What stage 5 could not check

Still nothing ran against GEO. The own-folder rule assumes GEO's sample file links have the form
`ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1234nnn/GSM1234567/suppl/GSM1234567_<name>` (and `GSMnnn` for accessions
of three digits or fewer), as I remember it: a real record whose links differ gets no allowed file and cannot be
pinned until a reviewed change. The one exception to the plain-key rule, `contact_zip/postal_code`, is GEO's own
field name from memory; if real files hold other field names with punctuation, every such dataset is refused, text
withheld, until a reviewed change lists them. A real SOFT file holding a form feed or U+2028 inside a field is refused
the same way. The Wnt-target review terms (Tbx6, Sp5/8, T, Mixl1, Cdx2) will cost review round trips on series that
name them as markers.
