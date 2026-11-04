# virtual-embryo-challenge

Our entry to the [Virtual Embryo Challenge](https://virtualembryo.ai/challenge) (NeurIPS 2026), Track 1 (Human Team).

| Task | Predict | Test target |
|---|---|---|
| T1 temporal | a population of single cells, 32,285 genes | E12.5 heart-centred dissection |
| T2 spatial-temporal | cells with 500-gene expression and 3D coordinates | heart E12.5 (extrapolation), embryo E7.75 (interpolation) |
| T3 perturbation | the mutant embryo's cells, expression and 3D coordinates | beta-catenin KO at E8.75 |

Final submissions are due **through 2026-12-02 (UTC)**; by then the team captain also marks, on the account page,
which submission counts on each board (`docs/P3_RUNBOOK.md`, section 5). Start with `AGENTS.md` (it opens with the team
captain's guiding principle of 2026-11-04), then `docs/CHALLENGE.md` and `docs/PLAN.md`.

## Setup

```bash
uv sync                    # local: core + dev tools (no GPU stack)
uv sync --group gpu --group jax   # DGX: CUDA torch (cu130), scvi-tools, moscot
./init.sh                  # lock check, lint, tests
```

## Layout

```
src/vec/        package: io, boards, evaluate (local board emulator), baselines, cli
boards/         public board specs from the challenge site (panels, anchors)
configs/        experiment configs
results/        registry.csv, one row per scored run
docs/           challenge reference, plan, data policy, compute, research notes
data/           competition and external data (gitignored, DGX only)
runs/           run outputs (gitignored)
```

Competition data is unpublished and must not leave the team; see `docs/DATA_POLICY.md`.
