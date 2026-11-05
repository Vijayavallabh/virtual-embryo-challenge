---
task_id: c
role: Bioinformatics curator and research-integrity auditor (provenance of pathway, gene-set and motif resources)
status: complete
sources_found: 19
---
## Sources
[1] The Molecular Signatures Database (MSigDB) hallmark gene set collection (Liberzon et al., Cell Systems 2015) | https://pmc.ncbi.nlm.nih.gov/articles/PMC4707969 | Source-Type: academic | Accessibility: public | As Of: 2015-12 | Authority: 9/10
[2] MSigDB mouse collections: collection details (GSEA-MSigDB) | https://www.gsea-msigdb.org/gsea/msigdb/mouse/collection_details.jsp | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 9/10
[3] License Terms for MSigDB released after April 2017 | https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 9/10
[4] Guide to GO evidence codes (Gene Ontology Consortium) | https://geneontology.org/docs/guide-go-evidence-codes/ | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 9/10
[5] Evidence codes (GO wiki, Texas A&M) | https://gowiki.tamu.edu/wiki/index.php/evidence_codes | Source-Type: community | Accessibility: public | As Of: 2026-11 | Authority: 5/10
[6] Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities (CollecTRI; Müller-Dott et al., Nucleic Acids Research 51(20), 2023) | https://academic.oup.com/nar/article/51/20/10934/7318114 | Source-Type: academic | Accessibility: public | As Of: 2023-10 | Authority: 9/10
[7] DoRothEA README (saezlab; rdrr.io mirror) | https://rdrr.io/github/saezlab/dorothea/f/README.md | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 6/10
[8] JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles (Nucleic Acids Research 2024) | https://pmc.ncbi.nlm.nih.gov/articles/PMC10767809 | Source-Type: academic | Accessibility: public | As Of: 2023-11 | Authority: 9/10
[9] UCSC genome-announce: JASPAR 2024 TFBS tracks for hg19, hg38, mm10 and mm39 | https://groups.google.com/a/soe.ucsc.edu/g/genome-announce/c/evkVNtpj4LA | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 7/10
[10] HOCOMOCO v12 help and downloads (the server redirects it to the v12 downloads page of the current HOCOMOCO site) | https://hocomoco12.autosome.org/help | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 8/10
[11] HOCOMOCO in 2024: a rebuild of the curated collection of binding models for human and mouse transcription factors (Nucleic Acids Research 2024) | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10767914/ | Source-Type: academic | Accessibility: public | As Of: 2023-11 | Authority: 9/10
[12] CIS-BP FAQ (Catalog of Inferred Sequence Binding Preferences) | https://cisbp.ccbr.utoronto.ca/faq.html | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 7/10
[13] ChIP-Atlas licence (DBCLS Database Archive) | https://dbarchive.biosciencedbc.jp/en/chip-atlas/lic.html | Source-Type: official | Accessibility: public | As Of: 2020-05 | Authority: 8/10
[14] ChIP-Atlas README, LATEST data release | https://dbarchive.biosciencedbc.jp/data/chip-atlas/LATEST/README_e.html | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 8/10
[15] Cistrome Data Browser 3.0 (Nucleic Acids Research) | https://pmc.ncbi.nlm.nih.gov/articles/PMC10767960 | Source-Type: academic | Accessibility: public | As Of: 2023-11 | Authority: 8/10
[16] ENCODE Data Use, Software, and Analysis Release Policies | https://www.encodeproject.org/about/data-use-policy/ | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 9/10
[17] KEGG legal: copyright and licensing | https://www.kegg.jp/kegg/legal.html | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 9/10
[18] Reactome licence | https://reactome.org/license | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 9/10
[19] WikiPathways licence terms | https://classic.wikipathways.org/index.php/WikiPathways:Copyrights | Source-Type: official | Accessibility: public | As Of: 2026-11 | Authority: 8/10

## Findings
1. The hallmark collection was built from 8,380 founder sets of MSigDB v4.0 collections C1-C6 (C2 mixes curated pathways with experimental signatures; C6 holds oncogenic perturbation signatures), and each hallmark kept only genes that tracked a phenotype contrast in public expression datasets, so the Wnt/beta-catenin hallmark (49 founders) cannot be certified free of beta-catenin perturbation data, and the mouse MH set is an orthology copy of it [1][2].
2. In the mouse MSigDB, CGP sets are UP/DN signatures from perturbation studies, CP sets come from BioCarta, Reactome and WikiPathways, M3 TFT sets are GTRD-pipeline promoter binding sets (-1000/+100 bp), and M5 GO sets are not filtered by evidence code [2].
3. MSigDB is CC BY 4.0 except the KEGG_LEGACY and BioCarta sets (licensed to the Broad with qualified permission) and KEGG MEDICUS-derived sets (CC BY-SA 4.0); KEGG's FTP is subscription-only and its API academic-only; Reactome data are CC0 (illustrations CC BY 4.0) and WikiPathways content is CC0 [3][17][18][19].
4. GO separates mutant-phenotype (IMP/HMP), genetic-interaction (IGI/HGI) and expression-pattern (IEP/HEP) evidence from direct assay (IDA/HDA) and physical interaction (IPI), and IMP covers knockout and overexpression experiments, so dropping those six codes removes perturbation- and expression-response-based annotations by rule [4][5].
5. CollecTRI merges ExTRI text mining, NTNU curation, TRRUST, SIGNOR, GEREDB, IntAct, GOA, TFactS, HTRIdb, CytReg, Pavlidis 2021 and DoRothEA A; each edge keeps its sources and pooled PMIDs but no species, stage, cell type or method, and mouse and rat edges were translated to human symbols and treated as human [6].
6. CollecTRI signs come from the balance of activating against repressing PMIDs, then the TF's regulon majority, then a default of activation (about half of all edges are default-signed and 86% are activating); co-regulators are kept as regulators, and perturbation data (KnockTF) serve only as a benchmark [6].
7. DoRothEA is for academic use only because it includes TRED, its mouse version is an orthology transform of the human regulons, and its maintainers now recommend CollecTRI instead [7].
8. JASPAR 2024 CORE matrices are curated from binding assays (ChIP-seq, ChIP-exo and DAP-seq via GTRD, PBM, compiled matrices) and admitted only with orthogonal literature support; they are served by web, REST API, Bioconductor and pyJASPAR, with genome-wide TFBS tracks for mm39 and mm10 on UCSC, and the article states no data licence [8][9].
9. HOCOMOCO v12 is rebuilt from GTRD ChIP-seq and HT-SELEX, ships mouse subsets in MEME, JASPAR, HOMER and TRANSFAC formats with threshold-to-P-value tables, and is distributed under WTFPL (usable as CC BY), while CIS-BP is free but mixes direct and DBD-similarity-inferred motifs (flagged) and withholds licensed TRANSFAC motifs [10][11][12].
10. ChIP-Atlas is CC BY-SA 4.0, covers mouse only on mm9/mm10, calls peaks at three q-value levels with TSS +-1/5/10 kb target-gene tables, and records stage, genotype and treatment only in free-text sample metadata; Cistrome DB 3.0 (hg38/mm10, per-sample regulatory-potential scores) states no licence, and ENCODE data are free for unrestricted use on release [13][14][15][16].

## First-principles analysis
Assumptions are marked A1-A6. Scorer formulas are read from the pinned veckit code (git cc3859a, `common/core_metrics.py`); no docstring number is used. No data file was read; the only computation is a synthetic numpy bracket (section G).

**A. What T3 can reward from a gene-level prior.** With pb = per-gene mean of log1p CP10k over cells, dp = pb(pred) - pb(WT) and dt = pb(KO) - pb(WT):
- de_score = (raw - chance)/(1 - chance). raw = the share of the truth's n_up up-genes and n_dn down-genes found among dp's top n_up and bottom n_dn; chance = the better of the same overlap for rankings by +pb(WT) and -pb(WT). It is rank-only: a prior moves it only by pushing genes across the n_up/n_dn cut-offs.
- de_direction = partial Spearman correlation of dp and dt given pb(WT).
- severity = log(beta), beta = <dt,dp>/<dt,dt> over truth-DE genes, floored at log(10^-3) when beta <= 0 or r^2 <= 0.01 (r^2 of dp on beta*dt over those genes).
- skill = d_f/(d_f + d), d = distance of the value from the ceiling, d_f = floor-to-ceiling distance.

**B. Eligibility and lineage dilution (derived).** A gene is truth-DE only if |dt_g| >= 0.25 (code constant). log1p values are non-negative, so a down-call needs pb_g(WT) >= 0.25. If the knockout acts only in Mesp1-lineage cells, dt_g = sum over lineage types t of f_t * dm_gt, roughly l_g * delta_g, where l_g is the lineage share of g's WT pseudobulk and delta_g the within-lineage change, so a gene can be called only if l_g * delta_g >= 0.25. The useful prior is therefore

  u_g = s_g * l_g * 1[pb_g(WT) >= 0.25 + margin],

with l_g and pb_g computed by code from released WT E8.75 only. A raw genome-wide propensity s_g spends its rank budget on genes that can never be scored. The panel also lacks Ctnnb1 and the classic readouts (`docs/PLAN.md` 4.4), so a cell-level "Wnt activity proxy" cannot be read off the panel itself: without T1 scRNA (closed, A15) it collapses into a propensity-weighted panel score, the same object as u_g.

**C. Information bound (derived).** Let T_g mark the panel genes that really respond within the lineage (prevalence pi; A1: pi = 0.05-0.10) and let the score have AUC A for them (A2). For a Gaussian score, corr(s, T) = d*sqrt(pi(1-pi)) / sqrt(1 + d^2*pi(1-pi)), with d = sqrt(2)*PhiInv(A). At pi = 0.08 this is 0.10 (A = 0.6), 0.20 (0.7) and 0.31 (0.8). The prior's correlation with the truth delta is r_p = corr(s,T)*sqrt(S_w), where S_w is the share of dt's variance carried by those responses (A3: 0.1-0.4), so r_p lies between 0.03 and 0.19.

**D. Magnitude-signal trade-off (derived).** Add lambda*p to an entry dp0 whose rank correlation with dt is r0, with p uncorrelated with dp0, and set kappa = lambda*sd(p)/sd(dp0). Then

  r(kappa) = (r0 + kappa*r_p) / sqrt(1 + kappa^2).

- The prior helps only while kappa < 2*r_p*r0/(r0^2 - r_p^2), about 2*r_p/r0.
- The best gain is sqrt(r0^2 + r_p^2) - r0, about r_p^2/(2*r0), which is second order whenever the entry already carries direction signal.
- At kappa = 0.25 the change is Delta r = 0.24*r_p - 0.03*r0.
- Near the floor, de_direction's skill is about 1/2 + v/(4C) (C = its ceiling; A4: C about 0.6), worth about 10*Delta r board points.

So even the best case (r_p = 0.19, r0 = 0) moves de_direction by under +0.5 board points at kappa = 0.25. A larger kappa helps only if the entry truly has no signal, and it raises the costs in E.

**E. Severity cliff (derived).** On the truth-DE genes, a prior uncorrelated with dt adds variance to dp and scales the explained share r^2 by about 1/(1 + kappa_DE^2). An entry whose r^2 sits near the 0.01 gate drops to the floor, a loss of 25*(skill - 0.5) board points (up to about 10 from a well-scored severity). Down-pushes on genes that truly rise also lower beta.

**F. Motif statistics (derived).**
- TCF/LEF matrices are short HMG-box motifs. At p < 1e-4 per position and strand, a window of L bp carries about 2e-4*L chance hits: about 0.5 in a 2.5 kb promoter and about 20 in a 100 kb regulatory domain.
- Raw counts therefore track window length, GC/CpG content and repeat masking. Developmental regulators, with large conserved regulatory domains, collect hits for any motif.
- A specificity residual against shuffled matrices and unrelated motif families is mandatory.
- A5 (not sourced this session): binding of this family is enhancer-biased and gated by accessibility. A promoter-motif score's AUC for cell-type-specific responders is therefore likely 0.55-0.65; a binding-informed score from relevant contexts might reach 0.65-0.8.

**G. Synthetic bracket (my numpy calculation, no data).**
- Set-up: 500 genes and 250 truth draws per cell. Shared processing carries 45% of dt's variance and composition 20%; about 10% of genes are DE.
- Prior and entries: the prior is -max(z(s),0)*l*pb on eligible genes. The entry has either no transferred factor (rho_f = 0) or a half-transferred one (rho_f = 0.5), and its composition is half right in both cases.

Mean change in T3 board points at kappa = 0.25 (rho_f = 0 / rho_f = 0.5; the ceilings assumed at 0.6, 0.6 and 0.05):

| Prior AUC | S_w = 0.10 | S_w = 0.25 | S_w = 0.40 |
|---|---|---|---|
| 0.5 (noise) | -0.2 / -0.5 | -0.1 / -0.6 | not run |
| 0.6 (motif-like) | -0.1 / -0.3 | +0.2 / -0.3 to -0.4 | +0.2 / -0.3 |
| 0.7 | +0.2 / -0.1 | +0.3 to +0.4 / 0.0 | +0.5 / 0.0 |
| 0.8 (binding-like, optimistic) | not run | +0.6 / +0.2 | +0.8 / +0.4 |

- At kappa = 0.75, every cell with AUC <= 0.7 is negative when a factor is present (-1.6 to -3.3), and a noise prior costs about -1 even without one.
- A binary top-25 prior at kappa = 0.5 reaches +0.9 to +1.3, but only at AUC 0.8.
- The per-draw SD is 0.8-1.5 points, so a +0.4 mean clears +0.3 in only about half of the truth draws.

Message: a motif-only prior is worth about 0 (-0.4 to +0.2). A binding-informed one is worth +0.2 to +0.8 only under optimistic AUC and S_w. Either one loses points on an entry that already carries signal, or at any but a small magnitude.

**H. Provenance lattice (derived), from cleanest to excluded.**
1. Sequence and in-vitro binding (matrices from HT-SELEX or PBM, genome sequence, conservation): clean by construction. Matrices fitted to ChIP peaks keep only sequence preference, so treat them as clean as well.
2. In-vivo binding (ChIP of TCF/LEF factors, or of the held-out co-activator): not a perturbation, but tied to its context. A judgement call that needs coded filters for windows, genotypes and treatments.
3. Unperturbed co-expression (DoRothEA's inferred layer; WT co-variation in released data): clean of perturbation, but context-dependent.
4. Literature-curated or text-mined signed edges (CollecTRI, DoRothEA curated, TRRUST, SIGNOR, GO IMP/IGI/IEP and annotation extensions, transcriptional output nodes of KEGG, Reactome and WikiPathways): a sign almost always comes from a gain- or loss-of-function statement in the cited paper. For the held-out regulators these edges restate perturbation results: excluded.
5. Perturbation footprints and expression-refined collections (MSigDB CGP, C6, hallmarks; PROGENy-type footprints): excluded by construction.

How to tell a tainted regulon edge:
- An edge is tainted if (a) its regulator is the held-out co-activator, one of its DNA-binding TCF/LEF partners, Gata4 or Gata6 (the validation genotype carries Gata6 F/+), or any node whose loss phenocopies Wnt loss, and (b) its evidence is curated or text-mined.
- CollecTRI pools PMIDs per edge and records no species, stage, cell type or method, and it merges translated mouse and rat edges into human symbols [6]. Deleting one PMID therefore never removes an edge.
- The only clean action is to drop those regulators' regulons whole, by rule, without opening them. In DoRothEA, keep only their motif (TFBS) and co-expression edges, plus ChIP edges if the captain allows binding.
- A check that reads no content: CollecTRI's own columns show whether a sign is PMID-derived or defaulted. A default-signed edge carries no direction evidence, so even a clean one adds nothing to a direction prior.
- Extra W-class risk: because stage is stripped, any curated edge may rest on an in-window developmental finding, so for T1 and the T2 heart curated regulons are a window risk as well as a genotype risk.

**I. Decision table.**

| Resource | Provenance (how built) | Verdict | Licence | Access |
|---|---|---|---|---|
| JASPAR 2024 CORE vertebrates | Curated matrices from ChIP-seq/ChIP-exo/DAP-seq via GTRD, PBM and compiled matrices; orthogonal literature support required [8] | Allowed (in-vitro-derived matrices for the strictest reading) | Open access; no licence stated in the article [8]; a secondary catalogue says CC BY 4.0 (unconfirmed) | Web, REST API, JASPAR2024 (Bioconductor), pyJASPAR, flat files; UCSC TFBS tracks for mm39 and mm10 [8][9] |
| HOCOMOCO v12 (v14 current) | Rebuilt from GTRD ChIP-seq and HT-SELEX; mouse subset [10][11] | Allowed | WTFPL, usable as CC BY [10] | Downloads (PWM/PCM/PFM; MEME, JASPAR, HOMER, TRANSFAC; threshold tables); Zenodo [10] |
| CIS-BP | Direct (PBM, SELEX) and DBD-similarity-inferred motifs, flagged; licensed TRANSFAC motifs withheld [12] | Allowed, direct motifs first | Free; licence text not found [12] | Bulk download |
| Own scan of GRCm39 | Genome sequence + matrices (+ conservation) | Allowed (cleanest route) | Open genome and annotation (not verified this session) | Our committed code |
| ChIP-Atlas | Reprocessed public ChIP-seq; peaks at three q-value levels; TSS +-1/5/10 kb target tables; free-text metadata [14] | Judgement call (TCF/LEF binding; the held-out protein flagged separately, default out); only after coded window/genotype/treatment filters | CC BY-SA 4.0 (share-alike on derivatives) [13][14] | dbarchive files; mm9/mm10 only [14] |
| ENCODE | Consortium ChIP-seq, mostly human cell lines | Judgement call (binding; human to mouse by orthology) | Unrestricted on release; cite accessions [16] | Portal/API, rate-limited [16] |
| Cistrome DB 3.0 | Reprocessed GEO/ENCODE/Roadmap; per-sample RP scores; ontology-coded cell type and tissue [15] | Judgement call; avoid until terms are known | Not stated [15] | Web and REST; no bulk route stated [15] |
| MSigDB M3 TFT (GTRD) | Promoter (-1000/+100) binding sets from the GTRD pipeline [2] | Excluded as shipped (aggregated samples cannot be filtered); rebuild per experiment instead | CC BY 4.0 [3] | GMT/JSON |
| MSigDB H/MH Wnt/beta-catenin hallmark | Founders from C1-C6 incl. perturbation signatures; refined on phenotype-contrast data; MH = orthology copy [1][2] | Excluded | CC BY 4.0 [3] | - |
| MSigDB M2 CGP, C6 | Perturbation signatures (UP/DN) [2] | Excluded (any Wnt, TCF, beta-catenin or Gata set; the rest irrelevant) | CC BY 4.0 [3] | - |
| MSigDB M2 CP | BioCarta, Reactome, WikiPathways [2] | As its source resource | CC BY 4.0; BioCarta and KEGG_LEGACY under Broad licence; KEGG MEDICUS CC BY-SA 4.0 [3] | GMT |
| KEGG Wnt map | Hand-drawn literature map: components plus a few transcriptional output nodes of unrecorded provenance | Components allowed but carry no target signal; output nodes excluded | Website free to academics; FTP subscription-only; API academic-only [17] | KEGG REST (academic), KGML |
| Reactome | Curated reactions with references (mouse by orthology inference; not verified this session) | Signalling machinery allowed; outputs of transcription events regulated by the TCF/LEF complex excluded | CC0 data, CC BY 4.0 illustrations [18] | Downloads, Content Service |
| WikiPathways | Community-curated GPML diagrams [19] | Components allowed; target nodes excluded | CC0 [19] | GPML and GMT releases [19] |
| GO (mouse annotations) | Annotations with evidence codes [4] | Allowed after dropping IMP, IGI, IEP, HMP, HGI, HEP and annotation extensions; IBA/ISO/ISS a judgement call (propagated); role terms, not targets, so no use for a propensity score | CC BY 4.0 per GO's stated policy (not verified this session) | GAF releases |
| DoRothEA | Curated, ChIP, TFBS and co-expression evidence, levels A-E; mouse by orthology [7] | Partly: TFBS and co-expression allowed; ChIP a judgement call; curated edges of the held-out regulators and partners excluded | Academic only (TRED) [7] | R package, OmniPath (check the per-evidence flags exist before filtering) |
| CollecTRI | Literature-curated and text-mined edges with pooled PMIDs; signs from PMIDs, regulon majority or default; mouse edges merged [6] | Excluded for the held-out co-activator, the TCF/LEF family and Gata4/Gata6 (whole regulons, unopened); other TFs allowed but irrelevant here | Not stated (article CC BY) [6] | OmniPath, decoupler |
| PROGENy-type footprints | Built from perturbation experiments by design (not researched this session) | Excluded | - | - |

**J. Recipe for a provenance-clean TCF/LEF propensity, computed in code.** Nothing is printed as a gene list; outputs are hashed files read by code, as with test-board counts.
0. Before any download: `EXTERNAL_DATA.md` rows (resource, version, licence, what is removed and how, T3 only), and the admission filter of lever L0 as committed code with a unit test.
1. Gene map: panel symbols to MGI accessions to Ensembl GRCm39 gene ids (pinned MGI and Ensembl releases). Unmapped genes stay NA; no hand fixes beyond MGI synonyms.
2. Windows: the union of [TSS - 2,000, TSS + 500] over annotated transcripts. A second arm uses a GREAT-style domain (basal -5/+1 kb, extension up to 500 kb). Pre-register one arm as primary.
3. Motifs: select the TCF/LEF family by its TF-family annotation in code, never by hand: JASPAR CORE vertebrates as primary, the HOCOMOCO mouse subset as replication.
   - Controls: 50 matrices from unrelated families, matched on information content and GC.
   - Shuffles: 20 dinucleotide shuffles of each family matrix.
4. Scan: FIMO or MOODS on both strands at p < 1e-4, with a local GC background. Exclude soft-masked repeats and merge overlapping family hits. Optionally weight hits by mm39 phyloP clipped at 0.
5. Score:
   - R_g = the number of hits, or the sum of -log10 p over hits.
   - Z_g = R_g against the shuffles in the same window.
   - epsilon_g = the residual of Z_g(TCF/LEF) on the median Z_g of the control families, log window length and CpG o/e, fitted genome-wide over protein-coding genes.
   - s0_g = the genome-wide percentile of epsilon_g.
6. Target-free sanity checks, pre-registered; failing any one stops the work:
   - Spearman(JASPAR, HOCOMOCO) over the genome >= 0.7;
   - the panel's share of genome-wide top-decile genes lies within 2-fold of 10%;
   - |Spearman(s0, WT pb)| < 0.2 over the panel.
7. Binding layer, only on the captain's yes: ChIP-Atlas mm10 experiments for TCF/LEF antigens; the held-out protein is a separate flag, out by default.
   - Drop embryonic samples whose stage or Theiler range touches (E7.25, E8.0), (E8.25, E8.75) or (E9.5, E13.5], and embryonic samples whose stage cannot be parsed.
   - Drop samples carrying any genotype or treatment token for knockout, knockdown, Cre/flox, CRISPR, dominant-negative or mutant, or for a pathway agonist or antagonist (this also removes Wnt-inhibitor steps of stem-cell cardiac differentiations).
   - Lift peaks to mm39 and compute, per experiment, RP_g = sum over peaks within 100 kb of exp(-(0.5 + 4d/100 kb)) (BETA-style decay, fixed in advance).
   - Convert RP_g to genome-wide percentiles, average them over admitted experiments to get s1_g, and set s_g = mean(s0_g, s1_g).
8. Panel weight: u_g = s_g * l_g * 1[pb_g(WT) >= 0.25 + margin], with pb_g and l_g taken from released WT E8.75 through the harmonised labels.
9. Use, by ruling only: a count-space factor exp(-lambda*u_g) on Mesp1-lineage cell types through the existing thinning primitive.
   - lambda is set from kappa = 0.25 of the entry's own delta SD and never tuned.
   - No cell is re-weighted, so the closed composition lever is untouched.
10. Analog: rerun steps 3-9 with the GATA family on WT E8.75, scored against Gata4 (lever L1's validation route).

**K. Validation logic.**
- No released data contain a Wnt loss, so no counting fold can score a beta-catenin-specific prior.
- The only check is the method family "propensity of the deleted factor's DNA-binding effector", run on the one released TF knockout (Gata4, GATA family). It would be a cross-knockout fold: WT E8.75 plus the external prior against both Gata4 replicates, with nothing from Gata4 tuned.
- That analog is the most favourable case: the deleted gene binds DNA itself, unlike the co-activator. Its other caveats are the Gata6 heterozygosity and the pre-registered check D9, which found no type-specific residual on Gata4 (residual share -0.04 and -0.15, replicate residual Spearman 0.73; `experiments/t3_decomp/results/gata4/gata4_decisions.json`).
- A pass is weak support; a fail is strong evidence against the bet.

## Candidate levers
**1. L1: gated TCF/LEF propensity factor in Mesp1-lineage cells (motif plus binding layer).**
- **Board and group:** T3:bcatenin. Response identity (de_score, 0.30) first, then response direction (de_direction, 0.25). Response identity is T3's weakest group in `experiments/harness/TRIES.md` (15.18 Overall points left across knockouts).
- **Mechanism:** u_g from recipe J, steps 1-8 (motif layer s0 plus the ChIP-Atlas layer s1), applied as exp(-lambda*u_g) to Mesp1-lineage types at kappa = 0.25.
- **Expected effect:** +0.2 to +0.8 T3 points if the binding-informed AUC is >= 0.7, S_w >= 0.25 and the entry carries little direction signal (bracket G); about 0 or negative otherwise. Weighted by the chance that the analog passes (A6: about 0.25) and that the scenario is favourable, EV is about +0.1 to +0.3 T3 (= Overall) points, SD about 1. This is not a top-3 lever.
- **Legality:**
  - It touches no closed lever: expression facet only, per-type count-space factors, no cell re-weighting, no T1 scRNA, no leaderboard input.
  - It departs from the frozen T3 rules and from `docs/PLAN.md` 4.4 item 2 ("released data only"), so it needs the captain's ruling and a `METHOD_SUMMARY.md` update.
  - The binding layer is the captain's judgement call (TCF/LEF binding; the held-out protein's binding is out by default).
  - The window union and the genotype/treatment filters are coded; each source gets its `EXTERNAL_DATA.md` row first; no list is hand-written.
  - The organisers' 2026-10-15 ruling allows generic motif and gene-set resources not built from perturbing the held-out genes.
- **Validation route:**
  - The analog of section K, registered before its run as one T3 try (the fourth), at kappa = 0.25, 10 seeds, on both replicates. It is kept only on a paired gain above 0.3 with at least 8 of 10 seeds ahead, sitting on de_score and de_direction.
  - The Gata4 -> Mab21l2 fold has no TF motif, so the captain must rule the analog the only counting fold.
  - TRIES rule 5 prunes rows expected below +0.3: the motif-only analog would likely be pruned, and the binding-informed one is borderline.
  - The beta-catenin application stays a blind transfer. Because the gate reads Gata4, it can enter only s1 (allowed branch), never s2.
- **Effort:** 3-4 days (scan 0.5; binding filter and RP 1.5; analog registration, run and review 1-1.5).
- **Main risk:** the severity r^2 cliff (E); double counting with the frozen depletion of second-heart-field and outflow-tract states, whose markers may be motif-rich (overshoot); a GATA-analog pass that is really a motif-by-processing coincidence (the D11b lesson); free-text ChIP metadata letting an in-window or mutant sample through.

**2. L0: provenance gate for any gene-set, regulon, motif or binding input (all boards; prerequisite to L1 and L2).**
- **Board and group:** all boards; no scored group (compliance).
- **Mechanism:** a coded admission rule and register template. It admits only sequence, in-vitro binding or unperturbed resources. It refuses anything built from Ctnnb1, Gata4, Gata6, TCF/LEF or Wnt-pathway perturbation in either direction, and anything that keeps protected-window samples (window union, Theiler ranges included). It records each licence.
  - Never admitted: hallmark, CGP or C6 sets, PROGENy-type footprints, and curated signed edges or regulons of the held-out regulators.
  - GO enters only from GAF files with the six codes and annotation extensions dropped.
  - No output is ever printed as a gene list.
- **Expected effect:** 0 points. It removes the disqualification tail that a naive build would carry: the common tools default to exactly the excluded resources (hallmark Wnt, CollecTRI, PROGENy).
- **Legality:** inside every rule; reads no target.
- **Validation route:** target-free audit, plus a unit test that the filter rejects each excluded class.
- **Effort:** 0.5 day.
- **Main risk:** over-strictness may drop a clean resource; negligible.

**3. L2: motif-only propensity (sequence layer s0 only) as an ungated small-kappa blind bet in T3 s2.**
- **Board and group:** T3; response identity and direction.
- **Mechanism:** L1 without the binding layer, placed in s2 (WT cells plus composition). By D, a prior pays first-order only on an entry without direction signal of its own, which s2 is when the factor bet fails. Being ungated, it keeps s2 free of Gata4-derived decisions.
- **Expected effect:** -0.1 to +0.2 in the bracket's motif-like cells without a factor, weighted by the chance that s2 is the counted upload: EV about 0 to +0.1.
- **Legality:** as L1 without the binding call; it changes a frozen s2 rule (`docs/S2_RULES.md`), so the captain's ruling is needed.
- **Validation route:** captain's blind bet only (an analog gate would make s2 depend on Gata4).
- **Effort:** 1.5-2 days.
- **Main risk:** near-zero signal with a small but certain noise cost; it spends the captain's attention on about 0.

**4. L3 (closed): a Wnt activity score from hallmark, CGP, PROGENy or CollecTRI/DoRothEA-curated regulons, for T3.**
- Closed by provenance (classes B and G): the hallmark is expression-refined from founders that include perturbation signatures [1]; CGP sets are perturbation signatures [2]; curated regulons of the held-out regulators cannot be cleaned edge by edge [6].
- It is also mis-specified: the classic readouts are not on the panel.

**5. L4 (closed): regulon or motif smoothing of T1 and T2-heart expression directions.**
- Closed: no counting check exists for a change to Feng's direction (TRIES, next-try item 5).
- Curated regulons carry stage-less developmental findings (W risk) [6].
- A motif-only regulon is too weak to denoise a direction whose sign agreement is already high (digest 4.3).

**Ranking by expected value to a top-3 place:** L1 (about +0.1 to +0.3), then L0 (0 points; do it first if L1 or L2 goes ahead), then L2 (about 0 to +0.1); L3 and L4 are closed. Together this resource class adds at most about +0.3 Overall in expectation and cannot supply the 15-25 points needed. Its main value is keeping a tainted Wnt score out of the entry.

## Deep Read Notes
- [6] CollecTRI (NAR 2023), read to about 95% of the article. It merges ExTRI, NTNU curation, ten curated databases, Pavlidis 2021 and DoRothEA A, harmonised to human symbols (Ensembl 109), with mouse and rat edges translated and treated as human. Each edge keeps its sources, per-source mode and pooled PMIDs; edges without a PMID are dropped; there is no species, stage, cell-type or method field. Signs come from the PMID balance, then the regulon majority, then a default of activation (19,541, 1,904 and 21,730 edges; 86% activating); co-regulators are retained. KnockTF perturbations serve only as a benchmark; the regulons are generic, and no resource licence is stated.
- [1] Hallmark collection (Cell Systems 2015), full text. Founders are 8,380 MSigDB v4.0 sets from C1-C6, clustered by Jaccard distance, then curator-assigned themes and the union of each cluster's genes. Refinement keeps up-regulated genes whose information coefficient against a phenotype is significant in GEO/ArrayExpress datasets (FDR 0.01), followed by ssGSEA validation on an independent dataset. The Wnt/beta-catenin hallmark has 49 founders and a one-line "canonical beta catenin pathway" description; its MSigDB page links the founders and refinement datasets (not opened).
- [8] JASPAR 2024 (NAR), full text. CORE matrices come from ChIP-seq, ChIP-exo and DAP-seq via GTRD, PBM and compiled matrices, admitted only with orthogonal literature evidence (otherwise UNVALIDATED). Access is by web, REST API, the JASPAR2024 Bioconductor package, pyJASPAR and flat files, with genome-wide TFBS tracks for human and mouse in UCSC (mm39 and mm10 per [9]). No perturbation data build the matrices, and no data licence is stated.
- [14] ChIP-Atlas README (with [13]). Mouse is on mm9/mm10 only, so mm39 needs a liftOver. Peaks come at q 1e-5, 1e-10 and 1e-20 per experiment (BED, bigBed, bigWig), with files merged per antigen and cell-type class and target-gene tables at TSS +-1/5/10 kb. The experiment list gives ids, assembly, antigen class and antigen, cell-type class, type and description, and author titles and free-text metadata; there are no structured stage, genotype or treatment fields. Licence: CC BY-SA 4.0, with share-alike derivatives.

## Gaps
- **Licences not confirmed:** JASPAR's data licence is not in the 2024 article or the UCSC notice, and the CC BY 4.0 in a secondary catalogue was not confirmed on the JASPAR site. The CIS-BP licence text and the Cistrome DB terms were not found.
- **Not verified this session:** GO's licence, the MGI GAF terms and GO's annotation-extension rules (the GO Consortium wiki returned HTTP 403; not worked around); Reactome's mouse orthology inference; MSigDB download registration.
- **DoRothEA evidence flags:** the per-edge evidence flags were not seen (two vignette copies returned 404). Confirm the curated/ChIP/TFBS/co-expression flags in the installed package or the OmniPath export before relying on them.
- **PROGENy:** not searched; its exclusion rests on its known design as a perturbation footprint.
- **Motif and binding AUC:** no source was searched (search budget) on how well promoter TCF/LEF motif scores or ChIP proximity predict responsive genes. A5 and the bracket's AUCs are my assumptions.
- **Process note:** one WebFetch saved a PDF of the hallmark paper into the harness's tool-results store. It was not opened; it is not a barred item, but the screener may want it listed as a residue.
- **Counter-claims against L1:**
  1. D9 on Gata4 found no type-specific residual for a programme to carry, so the within-lineage part L1 targets may be near zero.
  2. By D, a prior's gain is second order whenever the entry already carries direction signal, so even a good score adds about +0.3 at most.
  3. A GATA-analog pass could be a coincidence between motif-rich genes and the WT stage or processing offset, as most of D11b's carry turned out to be.
  4. The frozen composition already depletes second-heart-field and outflow-tract states, so L1 can double count and overshoot severity.
  5. A binding-informed AUC of 0.8 is optimistic; at 0.6-0.7 the bracket's mean is near 0.
- **Considered and dropped:**
  - Stripping a GATA-motif component from s1's Gata4 factor: the factor is measured in non-Mesp1 compartments where Gata4 is not deleted, so it holds little GATA-specific content, and s1 is frozen.
  - A cell-level Wnt axis built from the propensity-weighted panel: it reduces to L1 with more freedom, and resampling along it would touch the closed composition lever.

## Compliance log
none

## END
