# virtual-embryo-challenge

Team Mithril's entry to the [Virtual Embryo Challenge](https://virtualembryo.ai/challenge) (NeurIPS 2026, Human Team
track), and `vec`, the toolkit it runs on.

| Task | Predict | Test target |
|---|---|---|
| T1 temporal | single cells, 32,285 genes | E12.5 heart-centred dissection |
| T2 spatial-temporal | cells with 500-gene expression and 3D coordinates | heart E12.5 (extrapolation), embryo E7.75 (interpolation) |
| T3 perturbation | the mutant embryo's cells, expression and coordinates | beta-catenin knockout at E8.75 |

## Quick start for other teams

You need [uv](https://docs.astral.sh/uv/) and a challenge account. The competition data is not in this repository
(`data/` is gitignored): download it with your own token.

```bash
git clone https://github.com/Vijayavallabh/virtual-embryo-challenge
cd virtual-embryo-challenge
uv sync                                   # Python 3.12, CPU only

uv run vec boards                         # every board: gene panel, coordinates, cell limits, floor
# save your API token in ~/.config/vec/token (chmod 600), then:
uv run vec data manifest                  # the files your account can download
uv run vec data fetch --host "" --dest data/raw    # files land in data/raw/T1, T2 and T3

# a first valid prediction: a reference baseline, checked against its board
# (--n-cells keeps a random subset within the board's cell limits)
uv run vec baseline copy_last --board T1:val --last data/raw/T1/E9.5_RNA.h5ad --n-cells 5000 --out pred.h5ad
uv run vec check pred.h5ad --board T1:val          # prints OK (5000 cells)

# score it locally against the released validation answers (T1:val's target is E10.5)
uv run vec eval --pred pred.h5ad --board T1:val \
    --target data/raw/T1/E10.5_RNA.h5ad --ref data/raw/T1/E9.5_RNA.h5ad --seeds 0
```

## What `vec` does

| Command | What it does |
|---|---|
| `vec boards` | lists every board and its constraints |
| `vec data` | lists and downloads the released files your account can see |
| `vec baseline` | writes a reference baseline as a valid prediction: copy_last, pseudobulk_shift, wt_identity, shift_transfer or gene_ko |
| `vec check` | free format check of a prediction against a board |
| `vec eval` | local board emulator: the organisers' veckit metrics with the published skill map and weights, the T2 structure gate, the server's 10% working copies and several seeds. On released validation data it reproduces all 66 published validation anchors exactly. |
| `vec submit` | a free format-only check by the server (its official-upload path is tied to our team's own reference list) |

`./init.sh` checks the lock file, lints and runs the test suite.

## Our best submission

Overall 205.62 (T1 + T3 + the mean of the two T2 boards), one file per board:

| Board | File | Headline | Method in a line |
|---|---|---|---|
| T1:test | `T1_s2rcMSBm10.h5ad` | 74.75 | real E10.5 heart cells moved by released momentum, the Feng 2022 atlas direction and within-type transport |
| T3:bcatenin | `T3_c2_s2_g15.h5ad` | 64.22 | real wild-type cells scaled by one per-gene factor from the released Mab21l2 and Gata4 knockouts |
| T2:heart:test_extrap | `T2heart_organ_ts125_tk_half.h5ad` | 62.45 | real E10.5 heart cells moved along each lineage's trend across the released heart stages |
| T2:embryo:test_interp | `T2embryo_s2x_ovd.h5ad` | 70.85 | real E7.5 and E8.0 cells, the later stage turned into the earlier frame and thinned to the midway depth |

Every predicted cell is a real released cell: we draw real cells to a target composition, change their expression
only in count space, and place them in space. `METHOD_SUMMARY.md` gives each method in full and the command that
rebuilds each file.

## Layout

```
src/vec/       the package: boards, data, evaluate (local emulator), baselines, the shared core, task heads, CLI
experiments/   one folder per study; experiments/builds/ holds the builders of the submitted files
tests/         the test suite (./init.sh)
boards/        board specs from the challenge site (test-board anchors stripped)
configs/       experiment configs
docs/          challenge notes, the decision log (PLAN.md), data policy, runbook
results/       run registry and upload ledgers
```

Records: `docs/PLAN.md` dates every decision, `experiments/harness/TRIES.md` registers every try before it ran,
`EXTERNAL_DATA.md` lists every external source, and `docs/DATA_POLICY.md` holds the data-use rules. Agents working in
this repository start with `AGENTS.md`.
