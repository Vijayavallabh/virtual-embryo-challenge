# Upload sheet: every file that is built and waiting (state at 2026-11-24, 16:57 UTC)

For the captain. Every file is in `runs/manual_upload/`; each sha256 below was recomputed from the file on disk when
its row was written. No file here has a score. The allowance, on the captain's word of 2026-11-24: 20 uploads over 24
and 25 Nov, no cap per board (decision log 2026-11-24 (n), row R34). A lower score never replaces a board's best. The
order within a board is the order on record in `docs/S2_RULES.md` and the decision log; the captain may use any other
order. This sheet holds no score: where a file waits on a condition, the condition is given in words.

## The 20 uploads: what is due now and what waits on a score

Four of the 20 uploads were used on 2026-11-24, one per board (decision log 2026-11-24 (nn)); 16 are left. Nine files
are due now: their conditions are met or they have none, and with no cap per board none of them has to wait for
another. Six more wait on one reported score each. Three T1 files stand last by the captain's words of 2026-11-24 and
take what is left. That is 18 if every condition is met, more than the 16 left; conditions that fail free their
uploads. Between boards the captain decides; by the challenge's formula a point on T1 or T3 counts twice a point on
the heart or the embryo board.

**Due now (9), in the order of each board's table below:**

- T1: `T1_s2rcB30.h5ad`, `T1_s2rcCm20.h5ad`
- T3: `T3_ck_s2.h5ad`, `T3_k23_s2_e4.h5ad`, `T3_c1_s2_g2.h5ad`, `T3_c2_s2_g2.h5ad`
- Embryo: `T2embryo_s2.h5ad`, `T2embryo_s2x_ovhi.h5ad`
- Heart: `T2heart_organ_ts10_half.h5ad`

**Waiting on a score (6):**

- `T1_s2rcCm25.h5ad`: only if `T1_s2rcCm20.h5ad` scores strictly above C
- `T1_s2rcB2d15.h5ad`: only if `T1_s2rcCm20.h5ad` scores strictly above the 2.0 file; then `T1_s2rcB2d20.h5ad`: only
  if B2d15 scores strictly above the higher of the 2.0 file and Cm20
- `T1_s2rcB25d15.h5ad`: only if `T1_s2rcCm25.h5ad` scores strictly above both the 2.5 file and the 2.0 file; then
  `T1_s2rcB25d20.h5ad`: only if B25d15 scores strictly above Cm25
- `T1_s2rcCm30.h5ad`: only if `T1_s2rcCm25.h5ad` scores strictly above Cm20

**Last on T1, if uploads are left (3):** `T1_s2rcMSB.h5ad`, `T1_s2rcA15.h5ad`, then `T1_s2rcA20.h5ad` only if A15 scores
strictly above A.

**In work, not built:** none.

Eight files were uploaded and their scores reported on 2026-11-23, two per board (`docs/PLAN.md` decision log
2026-11-23 (cc) and (gg)), and four on 2026-11-24, one per board (decision log 2026-11-24 (nn)). The boards' best
headlines stand with `T1_s2rcB25.h5ad` (the 2.5 file), `T3_c2_s2.h5ad` (C2), `T2heart_organ_ts05_half.h5ad` (OTS05)
and `T2embryo_s2x_ov.h5ad` (O).

Five files are skipped by the captain's words of 2026-11-23 and 2026-11-24: they hold no place in any order below and
stand in a section of their own at the end.

## Task 1

Data used, for every Task 1 file in the order below: released challenge data and the outside heart atlas on recomputed
groups (GEO GSE193346, register row X002).

A search's next file is due once the score its condition needs has been reported; with no cap per board the due files
of different searches do not wait for each other. The first search's 2.5 file (`T1_s2rcB25.h5ad`) was uploaded on
2026-11-24 and holds the board's best; that search's next file, `T1_s2rcB30.h5ad`, is due. `T1_s2rcCm20.h5ad` (the 2.0
file's atlas step with g_mom's damping at 1.0 instead of 0.5) is due; its headline decides the third and fifth
searches' next files (`T1_s2rcCm25.h5ad`, `T1_s2rcB2d15.h5ad`), and Cm25's headline in turn decides the sixth and
seventh searches' (`T1_s2rcB25d15.h5ad`, `T1_s2rcCm30.h5ad`). The three files marked 'kept, last' were skipped and
then kept by the captain's words of 2026-11-24: by the main session's reading of those words they come after every
search file; the captain may place them anywhere.

**The fifth search stands by the captain's words of 2026-11-24** ("Yes: T1_s2rcCm20's headline may serve line 5; line
5 stands."; ruling row R35): `T1_s2rcCm20.h5ad`'s headline serves the fifth search as well as the third.

| Order | File | sha256 | Condition |
|---|---|---|---|
| first search | `T1_s2rcB30.h5ad` | `b639de29bd13d7005450ca53e7cd3dacb2e362a55beca8545937858e6582a024` | none (its condition is met: the 2.5 file scored strictly above the 2.0 file) |
| third search | `T1_s2rcCm20.h5ad` | `acce5f5f7097269bffebbfd583833e4c9c1c5ae77f1d5294dde6d13e197f863e` | none (its condition is met: C scored strictly above A) |
| third search | `T1_s2rcCm25.h5ad` | `c5fb96a64075b6fbfb2f3066604d622ea8b89ac3a5e4f8a99d0fea7d93011dfd` | only if Cm20 scores strictly above C |
| fifth search | `T1_s2rcB2d15.h5ad` | `fe4e7403d6257cc265a02d0a13d40cdb60a97e9fdc0647d6b5c7e318aedaed6c` | only if Cm20 scores strictly above the 2.0 file (`T1_s2rcB2.h5ad`) |
| fifth search | `T1_s2rcB2d20.h5ad` | `bd12c66312ba32ab938b4ee71bd5c79e58a154d6fc56434012e6dad907e668dd` | only if B2d15 scores strictly above the higher of the 2.0 file and Cm20 |
| sixth search | `T1_s2rcB25d15.h5ad` | `5bc2b32e4688439f29c54c1b1203325595d996c37df30ba2b9ff26f27a60f4c6` | only if `T1_s2rcCm25.h5ad` scores strictly above both the 2.5 file (`T1_s2rcB25.h5ad`) and the 2.0 file |
| sixth search | `T1_s2rcB25d20.h5ad` | `e6520dbee0fdc3f1336bc3f9de1e99e8a501de2d440b3f9085c88b5c7ecfb3f6` | only if B25d15 scores strictly above Cm25 |
| seventh search | `T1_s2rcCm30.h5ad` | `807fe097c88ee779430896002a34635aeffbc99b433c68ebdb7bd41d8af2cbbe` | only if `T1_s2rcCm25.h5ad` scores strictly above Cm20 |
| kept, last | `T1_s2rcMSB.h5ad` (combined momentum) | `78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2` | none; last on T1, after every search file |
| kept, last | `T1_s2rcA15.h5ad` (the second search's file) | `e2334137b54d98b04e9dece90213987e61fe9aa15a3fb8d85e6c262d8fb5ba75` | none (its condition is met); last on T1, after every search file |
| kept, last | `T1_s2rcA20.h5ad` | `372713686782be273a52427140eab652517e35b679f54fc2a969484dee730218` | only if A15 scores strictly above A |

Never uploaded: `builds/T1_s2.h5ad`.

## Task 3

| Order | File | sha256 | Condition | Data used |
|---|---|---|---|---|
| 1 | `T3_ck_s2.h5ad` (CK) | `e470aaeeab28918e157fb57cdc771122af35e46f0916cd1c6d63650c31faae26` | none (its condition is met: K23 scored strictly above the barred s2) | released challenge data (the wild-type E8.75 cells, the released Mab21l2 knockout and the released Gata4 validation knockout files the team unlocked) and two public single-cell knockout studies of other genes, GEO GSE204906 (Dnmt3a and Dnmt3b nulls, E8.5) and GEO GSE277839 (a Usb1 knockout, E9.5), each as knockout-over-its-own-control changes on the panel genes. No data on a perturbation of beta-catenin or Gata4 from outside the release |
| then | `T3_k23_s2_e4.h5ad` | `9934f1b104fafed0514cff7bb1dec3f5e362e6a5068ae32711cbfc67f62b90f1` | none (its condition is met: K23 scored strictly above the barred s2) | released challenge data (the wild-type E8.75 cells and the released Mab21l2 knockout) and two public single-cell knockout studies of other genes, GEO GSE204906 (Dnmt3a and Dnmt3b nulls, E8.5) and GEO GSE277839 (a Usb1 knockout, E9.5), each as knockout-over-its-own-control changes on the panel genes. No data on a perturbation of beta-catenin or Gata4 from outside the release |
| then | `T3_c1_s2_g2.h5ad` (the second search's file: C1 with the Gata4 table at double strength) | `e01e24facb6841a18a503655bcc5e9816304a993fe9a0dcfeca458c42ead1e16` | none; by its rule it comes after the files above (after CK and the E = 4 file when their condition is met), unless the captain places it earlier; with no cap per board it need not wait for a day | released challenge data only: the wild-type E8.75 cells, the released Mab21l2 knockout and the released Gata4 validation knockout files the team unlocked; no outside data |
| then | `T3_c2_s2_g2.h5ad` (the third search's one file: C2 with the Gata4 table at weight 2) | `014d2bf9d8ba7fdc56dcfba2a3b8aa77f0ea20bb806ebc2c1193e483d62712ac` | none; by the captain's yes of 2026-11-24 (row R41) it comes after the files above; the uploads are the captain's own acts | released challenge data only: wild-type E8.75 cells and the released knockout tables transfer_mab21l2_all and gata4_pooled, in C2's size-held product, gata4_pooled at weight 2 (applied powers 0.6159 and 1.2318; C2: c on both). No outside data, no gene set. |

## Task 2, embryo

Data used, for every embryo file: released challenge data only; no outside data.

By the embryo's order rule no score decides whether, when or in what order another embryo file goes; the captain may
use fewer, none, or another order.

| Order | File | sha256 |
|---|---|---|
| 1 | `T2embryo_s2.h5ad` | `41982d8afd784ecb07a4021077f701c3806c3bdb06eda92c2cb99fcc41b19b73` |
| 2 | `T2embryo_s2x_ovhi.h5ad` (O_hi) | `7e4d1eaafc40080d9d9c461d4b9396f0cee6353686538e79fb92ca31d1d2deeb` |

## Task 2, heart

Data used, for every heart file: released Task 2 heart stages and the unlocked E10.5 validation file only; no outside data.

`T2heart_organ_ts05_half.h5ad` (OTS05, built on D) was uploaded on 2026-11-24 and now holds the board's best; the
line's next file, OTS10, is due, and the line ends at OTS10 in any case. The full-trim file
`T2heart_core_ts10_half.h5ad` is skipped by the captain's word of 2026-11-24 and stands in the skipped section.

| Order | File | sha256 | Condition |
|---|---|---|---|
| 1 | `T2heart_organ_ts10_half.h5ad` (OTS10) | `ee9e0e68d62e8d1448f1dbbbce1146dc5fb4877e94d77fc289fa82744acb06a9` | none (its condition is met: OTS05 scored strictly above D); the line ends at OTS10 in any case |

## Skipped by the captain's words: not uploaded unless the captain says otherwise

The captain's words: "skip those three and update the upload sheet" (2026-11-23; `docs/PLAN.md` decision log
2026-11-23 (bb)) and "skip TS10, TS15 and TS20" (2026-11-24; decision log 2026-11-24 (p); TS15 and TS20 were never
built). The files stay in `runs/manual_upload/` as built; nothing was deleted or rebuilt.

| Board | File | sha256 | Why it is skipped, in the record's terms | Data used |
|---|---|---|---|---|
| Task 1 | `T1_s2rc_MS.h5ad` (plain momentum) | `6103d10535665634ee0cad3390590a5df132d1805682fb1d309bb19581fdaa5e` | the record's own range for its gain ends under the board's best headline | as the Task 1 files above |
| Task 1 | `T1_s2none.h5ad` (F) | `55eaee547ff18d7477da5b730c34e7de9ab74f85cc12fd01a9eff6450367611c` | on practice it never stands above the atlas file it was cut from, whose headline is under the board's best | released challenge data only; no outside data |
| Task 3 | `T3_refused_s1.h5ad` | `0394575a11fcd4f915e0daf4972e5a783b6b89efa065ed4f85fe2445cbbb1545` | its expected board score on record is under the board's best headline | released challenge data only (released wild-type stage files; its composition at the released Mab21l2 knockout's magnitude; no Gata4 validation file and no Gata4-derived number); no outside data |
| Task 2, heart | `T2heart_core_ts10_half.h5ad` (TS10, the full trim with the trend step at strength 1.0) | `cf53b0c0afe7d18b8f7c84d7595729e4fc46e15e6077e645b6dc0a2fbfed39af` | the captain's word of 2026-11-24, "skip TS10, TS15 and TS20" | released Task 2 heart stages and the unlocked E10.5 validation file only; no outside data |
| Task 1 | `T1_s2rc_pooled.h5ad` (R13's command on the pooled base) | `91341a5f8c71dcc103a47f9c19a5d7794b1358c1c642c5c9251624eb6feb8bdc` | the captain's word of 2026-11-24, "skip the pooled base file" | as the Task 1 files above |

## When a score comes back

Paste the headline per board to the main session. It records the headline only and applies each board's rule; the
conditions above then decide which file is next.
