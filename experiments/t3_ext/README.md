# T3 outside-knockout factor tables (stage 1)

Campaign branch `campaign/t3-outside-knockouts`. Under the captain's ruling of 2026-11-20 (`docs/PLAN.md` decision log
2026-11-20 (j), "All three tiers"), outside single-cell data on knockouts of OTHER genes may shape the Task 3 file.
This folder turns such data into per-gene factor tables in the frozen factor's form, so that a later thin builder can
take one in place of the Mab21l2 factor. Stage 1 builds the tables only. Nothing here enters a model, a check or an
upload until (1) each dataset has its row in `EXTERNAL_DATA.md`, (2) try K21's row is committed before any table is
scored, and (3) the captain rules on any change to an upload, with a `docs/S2_RULES.md` amendment.

## Files

| File | What it does |
|---|---|
| `datasets.yaml` | The datasets, their tier, study unit, design, control, stages, sample selectors and status. Only `cleared` ones are fetched. |
| `t3x_screen.py` | The held-out-gene screen (Ctnnb1 / beta-catenin, Gata4, their drugs; review terms such as Gata6, Apc, Axin, Gsk3). It has no switch. |
| `t3x_common.py` | Fixed constants, the dataset loader, the stage cap (E9.5 inclusive for every cell), shrinkage, the outside-the-repository guard. |
| `t3x_tables.py` | Finishing a log table into factors (rescale, clip, absent = 1), writing the pooled and per-group forms, the reference statistics. |
| `fetch.py` | `plan` screens every GSM of the series (and superseries) from GEO's family SOFT and picks the samples and per-sample files. `fetch` downloads exactly those files with sha256 and size in a manifest. |
| `process.py` | `reference` writes the WT E8.75 statistics. `run` turns one fetched dataset into one table per knockout gene plus `QC.json`. |
| `combine.py` | Tier tables (wnt, cardiac, unrelated), the all-knockouts table, one table per study, and leave-one-study-out tables. |
| `run_dgx.sh` | The detached DGX launcher, with exit files. |

## The rules each step keeps

- **Held-out genes.** Before any file is chosen, `fetch.py plan` screens every sample of every listed series, plus
  the series title, design, summary and file names. A held-out hit anywhere, a review term not cleared for that
  sample, or a record with no text refuses the whole dataset. `process.py run` screens each chosen record again
  before it reads a count. Per-cell label tables (MULTI-seq, gem.group, chimera sorts) go through
  `t3x_screen.screen_cells`. That function drops held-out, unassigned, negative and doublet cells. It refuses a
  table with no genotype column, a review term, or guides without a gene map. No stage-1 loader uses it yet, and a
  dataset that needs it is refused. Series-level files are never fetched in stage 1.
- **Phenocopy clause.** Wls (GSE165300), Hand2 (GSE126128) and Mef2c (GSE280543) have status `held_phenocopy`, and
  nothing fetches them. Each waits for one priced question to the captain.
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
3. **Log ratio.** On summed UMIs per gene: `log((x_ko+10)/N_ko) - log((x_ctl+10)/N_ctl)`. It is 0 below 20 UMIs in
   any arm, then clipped to log [0.25, 4]. Chimeras take [mut vs host] in knockout chimeras minus the same in
   wild-type chimeras.
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
meta/<ACC>/plan.json                       screen record, chosen samples and files, units (knockout gene -> samples)
raw/<ACC>/<GSM>/<file>, raw/<ACC>/MANIFEST.json   downloads with url, size, sha256
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

Use a campaign checkout at this branch's commit, synced only when nothing runs in it (`docs/COMPUTE.md`), and check
`df -h /home` first. The plan's download estimate is under 40 GB, but nobody has checked it.

```bash
mkdir -p ~/v/vec-ext/t3
cd ~/v/vec-campaign-d
setsid nohup bash experiments/t3_ext/run_dgx.sh "$PWD" > ~/v/vec-ext/t3/launcher.log 2>&1 < /dev/null &
```

Wait for `~/v/vec-ext/t3/logs/ALL.exit` with a background poller. A `plan_<ACC>.exit` of 2 is a refusal: its log
names the reason (a screen hit, a selector that found the wrong number of samples, a stage). Fix a selector only in a
reviewed commit, before any fetch of that dataset. Never relax the screen.

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
