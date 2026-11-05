# Compute

## Where things run

| Machine | Use |
|---|---|
| Local (shared, overloaded) | editing, git, `./init.sh`, quick checks on tiny data |
| `PrakashDGX_H2` (8x H100 80GB) | everything else: data processing, training, sampling, local board scoring on real data |

The DGX workspace is `~/v/virtual-embryo-challenge` (a git clone of this repo). GPU allocation (user, 2026-10-19):
until the end of October 2026, use only GPUs with free memory (other jobs run on the rest; check before launching and
never kill a job you did not start); from November 2026 all eight H100s are available to this project. The final
phase opened on 2026-10-27, so the release-week runs (Oct 27-31) still use free GPUs only. `/home` is
shared and 89% full (1.6 TB free on 2026-10-28): ask before writing more than 50 GB.

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
rsync -az --exclude .venv --exclude __pycache__ --exclude data --exclude runs --exclude .claude -e ssh ./ PrakashDGX_H2:v/virtual-embryo-challenge/
```

`.claude/` holds Claude Code's workflow worktrees (gitignored); it stays off the DGX, whose tree must stay clean.

Never sync while a build that writes prediction files runs on the DGX: the sync rewrites `.git` (HEAD's reflog), so
the file's build stamp records that the code may have changed mid-build, and an official upload of it is refused
(rebuild after the sync). Commit before building official candidates; the stamp must say the tree was clean.

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

Jobs that use numpy or scikit-learn need BLAS thread caps on the DGX: a run without `OPENBLAS_NUM_THREADS` crashed on
2026-10-28 ("OpenBLAS warning: precompiled NUM_THREADS exceeded", then "malloc(): corrupted top size"); keep it at 64
or less. The repo's shell runners set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` themselves (at
most 16 per process); a command started directly with `uv run` needs them on its command line, for example
`OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 uv run python ...`.

Large numpy jobs also need `NUMPY_MADVISE_HUGEPAGE=0`. NumPy asks for transparent huge pages on large arrays, and with
the DGX's THP settings (`madvise`, defrag `madvise`) each such allocation compacts memory first. On 2026-11-01, 27
payoff processes (about 520 GB, filling NUMA node 0) spent 99% of their CPU time in the kernel for 7 hours. The
compaction isolated about 435,000 pages a second and almost every migration failed, so a heart seed took 4.8 hours.
Killing the jobs stopped it at once. To check, sample `grep -E 'pgmigrate_fail|compact_isolated' /proc/vmstat` twice a
few seconds apart, and sample a process's user and system ticks with `awk '{print $14, $15}' /proc/<pid>/stat`. A job
whose system ticks grow while its user ticks stay flat is stuck. Keep a batch's memory well under one NUMA node
(about 500 GB): the account is shared.

Detached jobs over ssh: put `cd` on its own line (or use absolute paths) before `setsid nohup ... &`. In `cd DIR && ...
&& job & <wait loop>` the `&` backgrounds the whole chain, so the wait loop runs in the login directory: a waiter of
2026-11-04 polled the wrong folder for a day before it was stopped.

## Constraints

- `/home` on the DGX is 89% full (1.6 TB free, shared; `df -h /home` on 2026-10-28). Check it before large writes
  and ask before writing more than 50 GB.
- Stay inside `~/v`; the `prachh` account is shared with other lab members.
- No system packages and no changes to the shared conda base. uv only.
