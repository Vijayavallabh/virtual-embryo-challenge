# T3 outside-knockout factor tables (stage 3: the pinned sample list)

Campaign branches `campaign/t3-outside-knockouts` (stage 1), `-2` (stage 2) and `-3` (stage 3, this text). Under the captain's ruling of 2026-11-20 (`docs/PLAN.md` decision log
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
| `fetch.py` | `plan` writes a plan table of every sample of every screened series from GEO's family SOFT; `pin-template` prints a pinned block to paste; `check-pins` is the gate; `datasets-sha` prints the hash; `fetch` downloads exactly the pinned files of the pinned samples. |
| `process.py` | `reference` writes the WT E8.75 statistics. `run` turns the pinned samples of one fetched dataset into one table per knockout gene plus `QC.json`. |
| `combine.py` | Tier tables (wnt, cardiac, unrelated), the all-knockouts table, one table per study, and leave-one-study-out tables. |
| `run_dgx.sh` | The detached DGX launcher, with exit files: mode `plan` (the default) or `all` (behind the pinned gate). |

## The pinned sample list (the first line)

A text blacklist can never be shown complete, so it is not the first line of defence. The data in scope are five
cleared series (GSE205117, GSE168768, GSE186934, GSE137337, GSE204906), a few dozen samples in all, and the first line
is an explicit list of the samples two people have read.

1. **Plan** (metadata only). `fetch.py plan ACC` writes `meta/<ACC>/plan.json` and one table per screened series,
   `meta/<ACC>/plan_table_<SERIES>.txt`: every sample of the series, its accession, the verbatim text of every field
   of its GEO record (all are screened), the sha256 of those fields, the screen's verdict and hits, the per-sample
   files the dataset's patterns allow (and the other file names, never fetched), the parsed stage and chemistry, and
   the selectors' suggested arm, which chooses nothing. A dataset with a held-out hit anywhere is refused and none of
   its text is written: the table holds accessions, verdicts and the names of the fields hit, and nothing can be
   pinned from it.
2. **Read and pin.** The main session and a second, independent reader each read the plan tables (what they check is
   below). For each sample both agree on, `fetch.py pin-template` prints the entry, and the block is pasted into the
   dataset's `pinned` list in `datasets.yaml`:
   `{accession: GSMnnn, arm: ko|control, gene: <knockout gene>, stage: E8.5, files: [...], metadata_sha256: "..."}`.
   A control names the knockout gene it controls; a control that serves several knockouts is pinned once per gene.
   pin-template never writes `datasets.yaml` and refuses a held-out sample, a sample whose stage does not parse to one
   value, and a sample with no allowed per-sample file; it marks a sample with review terms (pin it only together with
   the `review_cleared` entries both readers agree).
3. **Commit.** Every edit of `datasets.yaml` changes its sha256. The same commit puts the new value
   (`fetch.py datasets-sha`) into `t3x_common.DATASETS_SHA256`; otherwise every step refuses the committed file, and
   `tests/test_t3_ext_screen.py::test_the_committed_dataset_list_hashes_to_its_pin_in_t3x_common` fails with the
   reason.
4. **Fetch and process** touch only pinned accessions. They refuse a dataset whose pinned block is empty; a plan that
   did not pass or was made from another `datasets.yaml`; a pinned accession that is not in its own series' plan table
   or belongs to another series; a pinned sample whose live metadata (when plan runs again) or stored metadata (fetch,
   process) no longer hashes to the pinned value; a pinned sample that no longer passes the screen; a pinned stage
   that differs from the record's; a pinned file not in the sample's plan entry; knockouts and controls whose stages
   do not pair, or whose chemistry differs where the dataset reads chemistry; a manifest fetched under other pins;
   a processed QC made under other pins (combine). Arms, genes and stages come from the pins only. No count file of
   an unpinned sample is opened (a test holds one beside the pinned ones and fails on any open of it).
5. **Per-cell genotype designs** (cell-level labels inside one sample: MULTI-seq, gem.group, chimera sorts) are
   excluded in this stage: a dataset with `cell_labels` cannot hold pins, and `process.py` refuses it.

**The canonical serialisation.** A pin's `metadata_sha256` is the sha256 of the UTF-8 bytes of

```python
json.dumps({"accession": GSM, "fields": F}, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
```

where `F` maps every field name of the sample's GEO family SOFT record (the `!Sample_` prefix removed, so `title`,
`characteristics_ch1`, `supplementary_file_1` ...) to the list of its values, each a string exactly as
`fetch.parse_soft` gives it (surrounding whitespace stripped, nothing else changed), in the order of the SOFT file
(`t3x_common.canonical_metadata`; `tests/test_t3_ext_pipeline.py` pins the exact bytes for one example). The hash is
taken from the sample's record in its own series' SOFT file.

**Held datasets.** `t3x_common.NEVER_FETCH` holds GSE165300, GSE126128, GSE280543, GSE259342 and E-MTAB-7325.
`dataset()` refuses them whatever a YAML says, `load_datasets()` refuses a YAML that gives any of them a status other
than `held_phenocopy` or `excluded_design`, and combine refuses any processed output of them. No CLI takes a
`--datasets` path: every step reads the committed file, checked against `DATASETS_SHA256`. Tests reach other YAML
content through a function argument only.

## What the two readers check in a plan table before pinning a sample

Each reader reads the whole plan table of the dataset's own series and of every other screened series (the
superseries), not only the rows the selectors suggest, and pins a sample only if both agree on all of these:

1. **Genotype and arm.** The sample's own text (title, source, characteristics, description, protocols, file names)
   states its genotype, and it is the arm being pinned: a knockout of the named gene, or that knockout's control as
   `control_definition` describes it. A label that names no gene (a bare 'F/F', 'cKO', 'mutant') is not pinned unless
   another field of the same record names the gene.
2. **No held-out or phenocopy label anywhere in the series.** No sample, series title, summary or design of any
   screened series mentions beta-catenin, Ctnnb1, Gata4, another allele of either, a drug or construct acting on
   them, or a perturbation of their pathway that could phenocopy them (Wnt, Wls, Apc, Axin, Gsk3, Tcf/Lef, Hand2,
   Mef2c, any GATA factor), in any spelling. The screen's verdicts help; they do not decide. If in doubt, nothing in
   the dataset is pinned and the question goes to the main session.
3. **No pooled design.** The sample is one embryo pool of one genotype: no MULTI-seq, cell hashing, CRISPR guide
   pool, chimera, mixed-genotype capture or per-cell genotype label (those are excluded in this stage).
4. **Same-study control.** The control comes from the same series and experiment, at the same stage, with the same
   chemistry where the dataset reads chemistry (GSE137337), and the relation is the one `control_relation` records.
5. **Stage.** The parsed stage in the table matches what the text says, is a single stage at or before E9.5, and is
   in `stages_used`.
6. **Assay.** Single-cell or single-nucleus RNA (the GEX part of a multiome), not ATAC, bisulfite, methylation, bulk
   or spatial; the pinned files are that sample's own per-sample matrices (no series-level files, no fragments).
7. **The hash** in the pinned entry equals the one printed for that sample in the table the reader read.

## The rules each step keeps

- **Held-out genes (the second line).** Before any file is chosen, `fetch.py plan` screens every sample of every listed series, plus
  every series-level field (title, summary, overall_design, type, file names and the rest; contact fields are
  skipped). A held-out hit anywhere, a review term not cleared for that sample (or, for a series-level field, that
  series accession), or a record with no text refuses the whole dataset. Clearing one review term never passes a
  second term in the same text, and nothing clears a held-out term. The spelling families caught, the review terms
  and the known over-refusals are listed in `t3x_screen.py`'s docstring; each has a test string, and a property
  test inserts every held-out spelling into every field of a clean record. A Gata4 reporter knock-in counts as
  held out. A label that names neither gene nor drug (for example 'Mesp1-Cre; F/F' alone) cannot be caught by any
  text screen: that is why two people read every plan table and pin samples before any fetch. `fetch.py fetch` and
  `process.py run` screen each pinned record again before a file is requested or a count is read. Since stage 3 the
  screen over-refuses by design: no allowlist on either side of the beta short names, catalogue wording and the like
  as review terms, typo forms as review terms, and an EXTENDED form (a Ctnnb1 or Gata4 typo, beta + 'cat',
  'catenine') voids every review clearance of the sample it appears in. Per-cell label tables go through
  `t3x_screen.screen_cells`, but no loader uses it and a dataset that needs it is refused. Series-level files are
  never fetched.
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
  checkout. Downloads are never unpacked. A per-sample tar.gz is read member by member in memory.
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
   number of libraries of an arm. d is 1 for germline or conditional knockouts against littermate or same-capture
   controls, and 0.5 otherwise. Stages are averaged with equal weight.
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
processed/<ACC>/QC.json                    cells per stage and arm, call shares, guards, panel genes, status
processed/<ACC>/<ACC>__<gene>.unit.json    log-space table for combine
processed/<ACC>/tables/<ACC>__<gene>/      finished per-dataset tables
tables/<name>/factors/<name>.csv           pooled form (gene, factor; reference gene order) + factors/factors.json
tables/<name>/<name>_per_group.csv         per-group form (group, gene, factor); a missing group = 1
tables/<name>/<name>_manifest.json         studies, genes, sign agreement, input sha256, 'checked' note
tables/COMBINE.json                        every table's files and sha256; left-out and not-cleared datasets
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

1. **Merge.** On this machine, merge `campaign/t3-outside-knockouts-3` into `main` after its review.
2. **Sync.** Only when nothing runs in `~/v/vec-campaign-c`: the `docs/COMPUTE.md` rsync line with the target
   `PrakashDGX_H2:v/vec-campaign-c/`, then on the DGX `git -C ~/v/vec-campaign-c checkout --detach <main's commit>`
   and `git -C ~/v/vec-campaign-c status --porcelain --untracked-files=no` (empty). Check `df -h /home`.
3. **Plan run** (metadata only), detached:
   `cd ~/v/vec-campaign-c && mkdir -p ~/v/vec-ext/t3 && setsid nohup bash experiments/t3_ext/run_dgx.sh "$PWD" plan > ~/v/vec-ext/t3/launcher_plan.log 2>&1 < /dev/null &`
   and wait for `~/v/vec-ext/t3/logs/ALL.exit` with a background poller. `plan_<ACC>.exit` 0 means the screen passed;
   2 is a refusal whose log names it (a review term, or a held-out hit with the text withheld).
4. **Read.** Copy `~/v/vec-ext/t3/meta/` (metadata only, no counts) to this machine if wanted. The main session and a
   second, independent reader each read every `plan_table_<SERIES>.txt` of each dataset and apply the checks above.
   A review refusal is resolved only by `review_cleared` entries both readers agree (sample, term, reason); a
   held-out refusal is never resolved, and that dataset is dropped.
5. **Pin.** For each agreed sample, `fetch.py pin-template ACC --dest ~/v/vec-ext/t3 --arm ko|control --gene G
   --gsm GSM...` (no network; it reads the stored plan) and paste its output under the dataset's `pinned:`.
6. **Commit.** `fetch.py datasets-sha`, put the printed value into `t3x_common.DATASETS_SHA256`, run the three test
   files, and commit `datasets.yaml` and `t3x_common.py` together on `main` (or on a branch merged into it), with the
   two readers named in the message. The `EXTERNAL_DATA.md` rows for the datasets that will be fetched go in before
   any table enters a model or a check.
7. **Sync again** as in step 2, so the DGX checkout holds the pinned `datasets.yaml`.
8. **Fetch, process, combine**, detached:
   `cd ~/v/vec-campaign-c && setsid nohup bash experiments/t3_ext/run_dgx.sh "$PWD" all "<pinned accessions>" > ~/v/vec-ext/t3/launcher_all.log 2>&1 < /dev/null &`.
   The launcher first runs the pinned gate (`fetch.py check-pins`; ALL.exit 6 if any accession has no pins or the
   YAML hash is not the pinned one), then the reference statistics, then per accession a fresh plan that checks each
   pinned sample's LIVE metadata against its pinned hash (`replan_<ACC>`), and only if it passes the fetch of the
   pinned files and the processing of the pinned samples; combine last. Wait on `logs/ALL.exit` in the background.
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
