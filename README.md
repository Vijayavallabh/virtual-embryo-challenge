# virtual-embryo-challenge

Team Mithril's entry to the [Virtual Embryo Challenge](https://virtualembryo.ai/challenge) (NeurIPS 2026, Human Team
track) and `vec`, the toolkit it runs on. We predict single cells (T1, heart E12.5), cells with 3D coordinates (T2,
heart E12.5 and embryo E7.75) and a knockout embryo (T3, beta-catenin at E8.75).

## Quick start

You need [uv](https://docs.astral.sh/uv/) and a challenge account; the competition data is not in this repository.

```bash
git clone https://github.com/Vijayavallabh/virtual-embryo-challenge && cd virtual-embryo-challenge
uv sync
uv run vec boards                                   # every board and its limits
uv run vec data fetch --host "" --dest data/raw     # your token in ~/.config/vec/token (chmod 600)
uv run vec baseline copy_last --board T1:val --last data/raw/T1/E9.5_RNA.h5ad --n-cells 5000 --out pred.h5ad
uv run vec check pred.h5ad --board T1:val           # OK (5000 cells)
uv run vec eval --pred pred.h5ad --board T1:val \
    --target data/raw/T1/E10.5_RNA.h5ad --ref data/raw/T1/E9.5_RNA.h5ad --seeds 0
```

`vec baseline` writes any of the five reference baselines, `vec check` is a free format check, and `vec eval` is a
local board emulator on the organisers' metric code that reproduces all 66 published validation anchors exactly.
`./init.sh` runs the lint and the tests.

## Our best submission

Overall 205.62 (T1 + T3 + the mean of the two T2 boards):

| Board | File | Headline |
|---|---|---|
| T1:test | `T1_s2rcMSBm10.h5ad` | 74.75 |
| T3:bcatenin | `T3_c2_s2_g15.h5ad` | 64.22 |
| T2:heart:test_extrap | `T2heart_organ_ts125_tk_half.h5ad` | 62.45 |
| T2:embryo:test_interp | `T2embryo_s2x_ovd.h5ad` | 70.85 |

Every predicted cell is a real released cell, drawn to a target composition and edited only in count space.
`METHOD_SUMMARY.md` gives each method and the command that rebuilds each file.

## Where things are

- `src/vec/` holds the package, `experiments/` one folder per study (`experiments/builds/` builds the submitted files)
  and `tests/` the test suite.
- The records: `docs/PLAN.md` (decisions), `experiments/harness/TRIES.md` (tries), `EXTERNAL_DATA.md` (external
  sources), `docs/DATA_POLICY.md` (data rules).
- Agents start with `AGENTS.md`.
