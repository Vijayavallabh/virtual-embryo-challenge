# virtual-embryo-challenge

Our entry to the [Virtual Embryo Challenge](https://virtualembryo.ai/challenge) (NeurIPS 2026), Track 1 (Human Team).

| Task | Predict | Test target |
|---|---|---|
| T1 temporal | a population of single cells, 32,285 genes | E12.5 heart-centred dissection |
| T2 spatial-temporal | cells with 500-gene expression and 3D coordinates | heart E12.5 (extrapolation), embryo E7.75 (interpolation) |
| T3 perturbation | the mutant embryo's cells, expression and 3D coordinates | beta-catenin KO at E8.75 |

Final submissions are due **25 November 2026, 23:59 AoE (26 Nov 11:59 UTC)**, moved from 2026-12-02 by the
organisers' notice of 2026-11-17; the team captain marks, on the account page, which submission counts on each board
by 25 Nov 12:00 UTC (`docs/P3_RUNBOOK.md`, section 5; the captain's ruling of [[DATE]]:
[[THE CAPTAIN'S WORDS AND DATE GO HERE BEFORE THIS IS MERGED]]). Start with `AGENTS.md` (it opens with the team
captain's guiding principle of 2026-11-04), then `docs/CHALLENGE.md` and `docs/PLAN.md`.

## Setup

```bash
uv sync                    # local: core + dev tools (no GPU stack)
uv sync --group gpu --group jax   # DGX: CUDA torch (cu130), scvi-tools, moscot
./init.sh                  # lock check, lint, tests
```

## Layout

```
src/vec/        package: io, boards, evaluate (local board emulator), baselines, cli, the shared core, the T1 and T2
                heads, the pseudo-task harness
experiments/    one folder per study (scripts, committed results); experiments/harness/TRIES.md is the tries ledger;
                experiments/builds/ holds the upload-file builders and the reference list of the ten upload files
boards/         public board specs from the challenge site (panels; test-board anchors are stripped at intake)
configs/        experiment configs
results/        registry.csv (one row per scored harness run), submissions.csv (the upload ledger), upload_gate.json
docs/           challenge reference, plan and decision log, data policy, runbook, compute, research notes
tests/          the test suite ./init.sh runs
data/           competition and external data (gitignored; on the DGX, and since 2026-11-12 on the local machine)
runs/           run outputs (gitignored); on the local machine also a checked copy of the ten reference upload files
                and of the T3 run folder (docs/COMPUTE.md)
```

Where things stand: `session-handoff.md` (what waits), `progress.md` (the log), `feature_list.json` (feature state),
`ASSUMPTIONS.md` (the bets), `METHOD_SUMMARY.md` (what each entry predicts and why).

Status on 2026-11-13: the entry for each board is chosen and its two upload files are built. After the captain's
rulings of that day (Task 1's second upload; the heart's expression step, which now follows each lineage's trend
over all released heart stages), all ten files (the first and second uploads and the Task 3 switch files) were
rebuilt three times on the DGX with identical bytes; their hashes are in
`experiments/builds/results/builds_dc7b6ae.SHA256SUMS`, and every one has passed the server's format check.
Nothing is uploaded yet. The first round may go
in from 21 Nov 00:00 UTC (moved from 2026-11-27 on the captain's ruling of [[DATE]]), each upload on the team
captain's go-ahead (`docs/P3_RUNBOOK.md`, section 5).

Competition data is unpublished and must not leave the team; see `docs/DATA_POLICY.md`.
