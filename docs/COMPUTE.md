# Compute

## Where things run

| Machine | Use |
|---|---|
| Local (shared, overloaded) | editing, git, `./init.sh`, quick checks on tiny data; since 2026-11-12 also the campaign's runs on otherwise-unused GPUs and cores (below) |
| `PrakashDGX_H2` (8x H100 80GB) | everything else: data processing, training, sampling, local board scoring on real data |

The DGX workspace is `~/v/virtual-embryo-challenge` (a git clone of this repo). GPU allocation (user, 2026-10-19):
until the end of October 2026, use only GPUs with free memory (other jobs run on the rest; check before launching and
never kill a job you did not start); from November 2026 all eight H100s are available to this project. The final
phase opened on 2026-10-27, so the release-week runs (Oct 27-31) still use free GPUs only. `/home` is
shared and 89% full (1.6 TB free on 2026-10-28): ask before writing more than 50 GB.

### The local machine during the campaign of 2026-11-12

By the captain's instruction of 2026-11-12 (`docs/PLAN.md` decision log 2026-11-12 (a)) the campaign of that day may
use this machine's otherwise-unused GPUs and cores; the DGX had no free GPU and saturated CPUs that day. What holds
here:

- The competition files were copied from the DGX to `data/` and checked against `data/raw/SHA256SUMS`. `data/` and
  `runs/` are owner-only (mode 700). The team API key stays local, as before, and no experiment uses it.
- The standard environment `.venv` is what `./init.sh` builds; the campaign's GPU work used a separate `.venv-gpu`
  (the locked `gpu` and `jax` groups), so `.venv` stays as it is.
- Run repository Python as `uv run --frozen --no-sync python -I -B ...` with `OMP_NUM_THREADS`,
  `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` at 16 or less and `NUMPY_MADVISE_HUGEPAGE=0`. Check the load first
  (`uptime`, `nvidia-smi`); other users' jobs run here and are never touched.
- A multi-task run starts each task as its own process, all together: provenance (commit and clean tree) is read at
  each task's start. Commit before starting, and do not edit tracked code while a run or a build is going.
- Test builds made here are checks only (`runs/builds_check_local/`): an official file is rebuilt at the freeze commit
  by `build_all.sh` at 16 BLAS threads. A rendered embryo file is byte-reproducible only at a fixed thread count
  (`experiments/builds/results/determinism_1f8956b.txt`).
- The DGX checkout was not synced during the campaign (it stayed at 6402b64). It was synced on 2026-11-13 to fd2683a,
  with no build running there, for the rebuild rehearsal of that day
  (`experiments/builds/results/rehearsal_fd2683a.txt`), and again that evening to 38616fe for the second rehearsal,
  after the captain's ruling on T1's second upload (`experiments/builds/results/rehearsal_38616fe.txt`). Sync it
  only when no build runs there.
- `.venv-gpu` exists on this machine only: the sync line below leaves it out.
- Since 2026-11-13 this machine holds a verified copy of what the uploads need from the DGX, so that a DGX outage in
  the upload week does not mean rebuilding under deadline: `runs/archive/builds_dc7b6ae/` (the ten reference files
  of the third rehearsal, with their stamps, sidecars and SHA256SUMS; check with `sha256sum -c SHA256SUMS` there;
  the copies of the first two rehearsals' sets, `runs/archive/builds_fd2683a/` and `builds_38616fe/`, are still on
  disk and are superseded: their two heart files are the ones the ruling of 2026-11-13 replaced) and
  `runs/t3_decomp/gata4_20261028T122425Z/` (the run folder every T3 build reads). Both are owner-only and ignored by
  git. This is the one exception to "bring back only scores and small figures" below. If a build changes, the copy
  of that file is stale until the new reference is copied and checked (`docs/P3_RUNBOOK.md`, "A copy off the DGX").
  The copy sets aside the "no lasting local copy" part of the choice of 2026-10-19 for these files; the captain
  confirmed it on 2026-11-13 ("Keep it (Recommended)", `docs/PLAN.md` decision log 2026-11-13 (k)).

### The campaign checkout on the DGX (2026-11-13)

Since the captain's renewal of the campaign on 2026-11-13 (`docs/PLAN.md` decision log 2026-11-13 (g)), campaign
runs on the DGX go in a second copy of the repository, `~/v/vec-campaign`, so that the reference checkout
(`~/v/virtual-embryo-challenge`), its `.venv` and `uv.lock` stay as the rehearsal left them.

- It has its own `.venv`, built once with `uv sync --frozen --group gpu --group jax` from the same `uv.lock`
  (7.3 GB), and its own `runs/`.
- `data/` there is a folder holding links to the reference checkout's `data/raw` and `data/external`. A bare link
  named `data` is not covered by `.gitignore`'s `data/`, and build stamps then say unclean
  (`experiments/builds/results/lib3_pair_ea17e39.txt`).
- Sync it with the line below and the target `PrakashDGX_H2:v/vec-campaign/`, only when nothing runs in it. A
  branch is checked out there detached (`git checkout --detach campaign/h7`): the copied `.git` remembers this
  machine's worktrees and refuses a plain checkout of a branch one of them holds.
- Start runs from a script file copied with `scp`, detached, and wait with a background poller for an exit file,
  as for the rehearsal.
- A file built there is a check, never an upload source: uploads are built in the reference checkout by
  `build_all.sh` at the freeze commit.
- There are three checkouts in all since the heart try of 2026-11-13: the reference one, `~/v/vec-campaign` and
  `~/v/vec-campaign-h7`, each with its own `.venv` from the same lock. Two campaign checkouts let runs at two
  commits go on at once (the generator tests: T1 at 069d6fc in `vec-campaign-h7`, the heart at 459916f in
  `vec-campaign`); nothing is ever checked out or synced in a checkout while a run is going in it.
- Two things stop a branch checkout after a sync, both met on 2026-11-13. The sync line has no `--delete`, so
  files that only an earlier branch tracked stay behind as untracked files and block the next branch that tracks
  them: compare each with the committed copy (`git show <commit>:<path> | sha256sum`) and remove it only if equal.
  And the sync copies this machine's uncommitted edits: commit first, or restore the copied file on the DGX
  (`git checkout -- <path>` there; this machine's file is not touched).
- A run that needs torch sees one GPU through `CUDA_VISIBLE_DEVICES`; GPUs 0 and 1 belong to other jobs. A count
  model fitted on an H100 reproduced bit for bit between two processes on two of the DGX's GPUs (the generator
  tests' `latent_sha256`), and differed from the fit on this machine's A100: a GPU fit is a per-hardware object,
  so a kept arm would need an archived fit, as T3's run folder is archived.

## What reproduces where (2026-11-13)

- **On one machine**, at 16 BLAS threads, every upload file rebuilds byte for byte: three builds of all ten on the DGX
  on 2026-11-13 (`experiments/builds/results/rehearsal_fd2683a.txt`), and repeated builds of each file on the local
  machine. T1's two files and the rendered embryo files also move with the thread count, so every build uses
  `build_all.sh`'s default of 16.
- **Across the two machines** seven of the ten files of the first rehearsal's set (fd2683a) are byte-identical; six
  of them are in today's set, and the seventh, T1's second file, was replaced by the ruling of 2026-11-13 (its
  successor has been built on the DGX only). The seven: T1's two, the embryo's two, the heart's two
  and T3's barred s2. The three T3 files that apply a per-gene factor (allowed s1, refused s1 with the wild-type-only
  factor, barred s1) hold the same cells in the same order with the same coordinates, but a few cells are drawn
  differently (3, 84 and 14 cells). The builder smooths each cell's profile over its nearest neighbours in a PCA, and
  the PCA's scores differ by up to 3e-4 between OpenBLAS's kernels: Haswell on this machine's EPYC 7742, SkylakeX on
  the DGX's Xeon Platinum 8468. A few cells then get other neighbours. T1's s1 as built until 2026-11-13 differed between the machines too.
- **So a hash checks a file only on the machine that built it.** The uploads are built on the DGX, and the reference
  list is the DGX's `runs/archive/builds_dc7b6ae/SHA256SUMS` (copy in the repository:
  `experiments/builds/results/builds_dc7b6ae.SHA256SUMS`). T1's second file as ruled on 2026-11-13 (the mid base)
  has been built on the DGX only; whether this machine gives the same bytes is not known. The two heart files as
  ruled that evening (the trend) came out with the same bytes on both machines.
- **Keep the DGX's environment as it is until the uploads are done.** The numerical libraries come from `uv.lock`; a
  change of the lock, of the Python environment there or of the processor can change those three files' bytes, and
  then the freeze-morning check names them. A change to `uv.lock` before the deadline (26 Nov 11:59 UTC;
  was 2 December, moved on the captain's ruling of [[DATE]] after the organisers' notice of 2026-11-17) needs a new
  rehearsal.

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
rsync -az --exclude .venv --exclude .venv-gpu --exclude __pycache__ --exclude data --exclude runs --exclude .claude -e ssh ./ PrakashDGX_H2:v/virtual-embryo-challenge/
```

`.claude/` holds Claude Code's workflow worktrees (gitignored); it stays off the DGX, whose tree must stay clean.
`.venv-gpu` is the campaign's local environment (2026-11-12; several GB): it stays off the DGX too. Before a sync,
check that no build runs there (`pgrep -af "[b]uild_all|[b]uild_test|[b]uild_submission"` prints nothing) and, after
it, that `git status --porcelain --untracked-files=no` is empty at the expected commit. To see first what a sync
would send, add `-n --itemize-changes`.

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
