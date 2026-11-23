# Upload sheet: every file that is built and waiting (state at 2026-11-23, 22:32 UTC)

For the captain. Every file is in `runs/manual_upload/`; each sha256 below was recomputed from the file on disk when
this sheet was written. No file here has a score. The order is the order on record in `docs/S2_RULES.md` and the
decision log; the captain may use any other order. At most two uploads a board a day until 2026-11-26 11:59 UTC; an
unused daily slot is lost; a lower score never replaces a board's best. This sheet holds no score: where a file waits
on a condition, the condition is given in words.

Three files are skipped by the captain's word of 2026-11-23: they hold no place in any order below and stand in a
section of their own at the end.

## Task 1

Data used, for every Task 1 file in the order below: released challenge data and the outside heart atlas on
recomputed groups (GEO GSE193346, register row X002).

| Order | File | sha256 | Condition |
|---|---|---|---|
| 1 | `T1_s2rcC.h5ad` (C) | `31477ef7df5d5374a74a8249c570217020c593916ab4d3780049b7b3f87b0620` | none |
| 2 | `T1_s2rcB2.h5ad` (the 2.0 file) | `52897fe375b8a81afffa788592047c30361824148b00a8176514385060518837` | none |
| 3 | `T1_s2rcMSB.h5ad` (combined momentum) | `78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2` | none |
| 4 | `T1_s2rcA15.h5ad` | `e2334137b54d98b04e9dece90213987e61fe9aa15a3fb8d85e6c262d8fb5ba75` | none (its condition is already met) |
| then | `T1_s2rcB25.h5ad` | `ccb7c457c0e8ac8ac668bcea17bc61ed22f414f58ab5a63e946cb2f08c4f4e77` | only if the 2.0 file scores strictly above the board's best file (B) |
| then | `T1_s2rcB30.h5ad` | `b639de29bd13d7005450ca53e7cd3dacb2e362a55beca8545937858e6582a024` | only if the 2.5 file scores strictly above the 2.0 file |
| then | `T1_s2rcA20.h5ad` | `372713686782be273a52427140eab652517e35b679f54fc2a969484dee730218` | only if A15 scores strictly above A |
| then | `T1_s2rcCm20.h5ad` | `acce5f5f7097269bffebbfd583833e4c9c1c5ae77f1d5294dde6d13e197f863e` | only if C scores strictly above A |
| then | `T1_s2rcCm25.h5ad` | `c5fb96a64075b6fbfb2f3066604d622ea8b89ac3a5e4f8a99d0fea7d93011dfd` | only if Cm20 scores strictly above C |

Never uploaded: `builds/T1_s2.h5ad`.

## Task 3

| Order | File | sha256 | Condition | Data used |
|---|---|---|---|---|
| 1 | `T3_c1_s2.h5ad` (C1) | `8718330866429dc089eaffd2f96f66936fe54134eba2d50f8126266e6a3e54b8` | none | released challenge data only: the wild-type E8.75 cells, the released Mab21l2 knockout and the released Gata4 validation knockout files the team unlocked; no outside data |
| 2 | `T3_c2_s2.h5ad` (C2) | `47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7` | none | as C1 |
| 3 | `T3_k23_s2.h5ad` (K23) | `77ff0b947f7f409c9ce171bdc59a44de3d3c2ef2f3e5a7332011a459ce9ee9c3` | none | released challenge data (the wild-type E8.75 cells and the released Mab21l2 knockout) and two public single-cell knockout studies of other genes, GEO GSE204906 (Dnmt3a and Dnmt3b nulls, E8.5) and GEO GSE277839 (a Usb1 knockout, E9.5), each as knockout-over-its-own-control changes on the panel genes. No data on a perturbation of beta-catenin or Gata4 from outside the release |
| 4 | `T3_ck_s2.h5ad` (CK) | `e470aaeeab28918e157fb57cdc771122af35e46f0916cd1c6d63650c31faae26` | only if C2 and K23 each score strictly above the board's current best | released challenge data (the wild-type E8.75 cells, the released Mab21l2 knockout and the released Gata4 validation knockout files the team unlocked) and the two studies named for K23, in the same form. No data on a perturbation of beta-catenin or Gata4 from outside the release |
| then | `T3_k23_s2_e4.h5ad` | `9934f1b104fafed0514cff7bb1dec3f5e362e6a5068ae32711cbfc67f62b90f1` | only if K23 scores strictly above the board's current best | as K23 |

## Task 2, embryo

Data used, for every embryo file: released challenge data only; no outside data.

| Order | File | sha256 |
|---|---|---|
| 1 | `T2embryo_s2x_cr.h5ad` (C) | `f20730da83084702c50c43d4d4aaa98be8ee801ea443472792172804fea4228e` |
| 2 | `T2embryo_s2x_crs.h5ad` (sized copy-rows) | `323f2acc32f731819b1a3cc0da32478a51d9e9a10f13dbd940773b113b5644fa` |
| 3 | `T2embryo_s2xs_ov.h5ad` (SO) | `77062082d97555961901c811d8a1aaf9b5b2403417c19a6988c816ff553f1165` |
| 4 | `T2embryo_s2.h5ad` | `41982d8afd784ecb07a4021077f701c3806c3bdb06eda92c2cb99fcc41b19b73` |
| 5 | `T2embryo_s2x_ovhi.h5ad` (O_hi) | `7e4d1eaafc40080d9d9c461d4b9396f0cee6353686538e79fb92ca31d1d2deeb` |

## Task 2, heart

Data used, for every heart file: released Task 2 heart stages and the unlocked E10.5 validation file only; no outside data.

| Order | File | sha256 | Condition |
|---|---|---|---|
| 1 | `T2heart_organ_half.h5ad` (D) | `1b4342c151c19dfc827d4b7a310c6a0fa8908f1e369836ebffb22c27c93b4f00` | none |
| 2 | `T2heart_core_ts05_half.h5ad` (TS05, the heart search's first file) | `b6b8155be9b5469a2dbc8943a59df567f52c7cbe013bbf0b55e4da706fda27b7` | none; after D |
| then | `T2heart_core_ts10_half.h5ad` (TS10) | `cf53b0c0afe7d18b8f7c84d7595729e4fc46e15e6077e645b6dc0a2fbfed39af` | only if TS05 scores strictly above the full trim (`T2heart_s2core_half.h5ad`), whatever D scores |
| then | `T2heart_organ_ts05_half.h5ad` (OTS05, the second heart search's first file) | `f918e76731aac18eda3b46c8a35c18f26ae00a1a6a1db36d33e7515704cc43ea` | only after D has been uploaded and its score reported; when a first-search file is due at the same slot, that file goes first |
| then | `T2heart_organ_ts10_half.h5ad` (OTS10) | `ee9e0e68d62e8d1448f1dbbbce1146dc5fb4877e94d77fc289fa82744acb06a9` | only if OTS05 scores strictly above D |

## Skipped by the captain's word of 2026-11-23: not uploaded unless the captain says otherwise

The captain's words: "skip those three and update the upload sheet" (`docs/PLAN.md` decision log 2026-11-23 (bb)). The
files stay in `runs/manual_upload/` as built; nothing was deleted or rebuilt.

| Board | File | sha256 | Why it is skipped, in the record's terms | Data used |
|---|---|---|---|---|
| Task 1 | `T1_s2rc_MS.h5ad` (plain momentum) | `6103d10535665634ee0cad3390590a5df132d1805682fb1d309bb19581fdaa5e` | the record's own range for its gain ends under the board's best headline | as the Task 1 files above |
| Task 1 | `T1_s2none.h5ad` (F) | `55eaee547ff18d7477da5b730c34e7de9ab74f85cc12fd01a9eff6450367611c` | on practice it never stands above the atlas file it was cut from, whose headline is under the board's best | released challenge data only; no outside data |
| Task 3 | `T3_refused_s1.h5ad` | `0394575a11fcd4f915e0daf4972e5a783b6b89efa065ed4f85fe2445cbbb1545` | its expected board score on record is under the board's best headline | released challenge data only (released wild-type stage files; its composition at the released Mab21l2 knockout's magnitude; no Gata4 validation file and no Gata4-derived number); no outside data |

## When a score comes back

Paste the headline per board to the main session. It records the headline only and applies each board's rule; the
conditions above then decide which file is next.
