# Session hand-off

State on 2026-11-26, the last upload day. The full previous hand-off is in git history:
`git show c10e079:session-handoff.md`.

## Where things stand

Uploads close on 2026-11-26 at 11:59 UTC. The captain marks one final submission per board on the account page
(until 2026-12-02, `docs/CHALLENGE.md`), and winners are announced on 2026-12-11.

Best headline per board, which is also the final pick (row R81 of `experiments/harness/TRIES.md`):

| Board | File | Headline |
|---|---|---|
| T1:test | `T1_s2rcMSBm10.h5ad` | 74.75 |
| T3:bcatenin | `T3_c2_s2_g15.h5ad` | 64.22 |
| T2:heart:test_extrap | `T2heart_organ_ts125_tk_half.h5ad` | 62.45 |
| T2:embryo:test_interp | `T2embryo_s2x_ovd.h5ad` | 70.85 |

Overall (T1 + T3 + the mean of the two T2 boards): 205.62.

- Every upload and its headline: `results/manual_uploads.csv`.
- Every built file and its sha256: `experiments/builds/results/`.
- Built files not yet uploaded: `docs/UPLOAD_SHEET.md`.

## What is open

- The final marks on the account page (the captain, by 2026-12-02).
- Row R80 (four embryo files) is ruled but unbuilt: its engineer was stopped and its unfinished work was not committed.

## How to resume

1. `./init.sh`
2. Read `AGENTS.md`, then the last entries of `docs/PLAN.md`.
