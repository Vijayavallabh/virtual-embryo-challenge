# Compute

## Where things run

| Machine | Use |
|---|---|
| Local (shared, overloaded) | editing, git, `./init.sh`, quick checks on tiny data |
| `PrakashDGX_H2` (8x H100 80GB) | everything else: data processing, training, sampling, local board scoring on real data |

The DGX workspace is `~/v/virtual-embryo-challenge` (a git clone of this repo). GPU allocation (user, 2026-10-19):
until the end of October 2026, use only GPUs with free memory (other jobs run on the rest; check before launching and
never kill a job you did not start); from November 2026 all eight H100s are available to this project. `/home` is
shared and about 95% full: ask before writing more than 50 GB.

## Setup on the DGX

```bash
ssh -o BatchMode=yes PrakashDGX_H2
cd ~/v && git clone https://github.com/Vijayavallabh/virtual-embryo-challenge.git
cd virtual-embryo-challenge
export PATH=$HOME/.local/bin:$PATH    # uv is not on the non-login ssh PATH
uv sync --locked --group gpu --group jax   # torch 2.14 from the cu130 index (driver 580 = CUDA 13.0), scvi-tools, moscot
uv run python -c "import torch; print(torch.cuda.is_available(), torch.cuda.device_count())"
```

Competition data goes in `~/v/virtual-embryo-challenge/data/raw/` on the DGX only (gitignored). Code reaches the DGX by
`rsync` from the local clone (no GitHub credentials on the shared account):

```bash
rsync -az --exclude .venv --exclude __pycache__ --exclude data --exclude runs -e ssh ./ PrakashDGX_H2:v/virtual-embryo-challenge/
```

The DGX shell exports `VIRTUAL_ENV=~/v/.venv`; uv ignores it and uses the project `.venv` (the warning is harmless).

## Running a job

```bash
# pick GPUs with free memory
nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader
# avoid clashing session names
tmux ls
# launch detached, with a log
tmux new -d -s vec-<name> 'CUDA_VISIBLE_DEVICES=<ids> uv run vec <command> 2>&1 | tee runs/logs/<name>.log'
```

Monitor with `tail -n 50 runs/logs/<name>.log`. Bring back only scores and small figures; checkpoints stay on the DGX.

## Constraints

- `/home` on the DGX is about 96% full (roughly 600 GB free, shared). Check `df -h /home` before large writes and ask
  before writing more than 50 GB.
- Stay inside `~/v`; the `prachh` account is shared with other lab members.
- No system packages and no changes to the shared conda base. uv only.
