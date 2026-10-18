# VEC 2026 (Human track) engineering report: environment, compute, experiment infrastructure

Compiled 2026-10-18. Scope: the uv environment on the PrakashDGX_H2 H100 box, library versions and conflicts, compute and memory estimates, and the reproducibility and evaluation tooling (run registry, CLI, board emulator, submission safety).

**How the claims were checked.** Package facts come from the PyPI JSON API, the download.pytorch.org simple index (including PEP 658 metadata sidecars) and the GitHub API, all fetched on 2026-10-18. Environment claims were validated by actually locking, syncing and importing on the DGX, in `~/v/vec2026-env-probe`, using uv 0.12.13, Ubuntu 22.04, glibc 2.35 and driver 580.178.04 (CUDA 13.0).

Timings were measured on the DGX CPU only. No GPU benchmarks were run, because all 8 GPUs were running your other jobs (details in §3.6). GPU contact was limited to a device query and one tiny matmul on GPU 7.

**Confidence labels:**
- **[High]**: verified by running it, or read directly from the source or metadata.
- **[Medium]**: derived from verified facts with a stated inference.
- **[Low]**: estimate or extrapolation.

---

## 0. Decisions and key findings

1. **Use Python 3.12 (`==3.12.*`)** [High]. The current numpy 2.5.4, scipy 1.18.1, anndata 0.13.4, scanpy 1.12.4, scvi-tools 1.5.1, cellrank 2.3.3, jax/jaxlib 0.11.2, flax 0.12.10, mudata 0.4.x and zarr 3.4.1 all require Python ≥3.12.
   - Under 3.11, uv silently falls back to numpy 2.4.6, scipy 1.17.1, anndata 0.12.19, scvi-tools 1.4.2 and jax 0.10.2.
2. **Use torch 2.14.1 from the cu130 index on the DGX and the cpu index locally** [High].
   - torch 2.14.1 is published only for the cpu, cu126, cu130 and cu132 variants. **cu128 stops at torch 2.11.0** and cu129 at 2.13.0.
   - Driver 580.178.04 is CUDA 13.0, so cu130 is the exact match.
   - On Linux, PyPI's default torch 2.14.1 is itself the CUDA 13.0 build (it requires `cuda-toolkit==13.0.3`).
   - Verified on an H100: `torch.cuda.is_available()` is true, compute capability 9.0, cuDNN 9.24, and matmul, POT-on-CUDA and torchcfm all work.
3. **The spateo/dynamo/Navigo stack cannot share an environment with the modern scverse stack** [High].
   - spateo's GitHub main pins `numpy<2` and `scipy<1.14` on Python ≥3.11.
   - dynamo 1.5.3 and Navigo pin `anndata<0.11` (dynamo also pins `pandas<3`).
   - cellrank 2.3.3 and mudata 0.4 need `anndata>=0.11`; anndata 0.13 and jax 0.11 need `numpy>=2.1`.
   - **Solution, validated:** one `pyproject.toml` and one `uv.lock`. The `spatial` dependency group is declared conflicting with the other groups, so uv resolves it in a separate fork (numpy 1.26.4, scipy 1.13.1, anndata 0.10.9) and the team syncs it into a second venv.
4. **Five pins or overrides were needed to make the stack resolve and import** [High, each reproduced on the DGX]:
   1. `statsmodels<0.15` in `spatial`. dynamo 1.5.3 fails to import with statsmodels 0.15.0: `multipletests` was removed from `statsmodels.sandbox.stats.multicomp`.
   2. ott-jax from git main. PyPI's ott-jax 0.6.0 (2025-11) crashes the online `PointCloud(batch_size=...)` under jax 0.11.2 (ott issue #700); main has the fix.
   3. `exclude-dependencies = ["torchdyn"]`. torchcfm's core API never imports torchdyn, and torchdyn pulls in torchvision, `torchcde<0.3` and `poethepoet<0.11`.
   4. `override-dependencies = ["jinja2>=3.1.6"]`. pygpcca 1.0.4, a cellrank dependency whose last PyPI release was in 2022, pins `jinja2==3.0.3`, but its code never imports jinja2.
   5. spateo from a GitHub **archive tarball** instead of `git+`. The repo stores two GRN CSVs in Git LFS, and the DGX has no `git-lfs`.
5. **Install veckit from git at `cc3859ab2a2305efc4544815d3f44f660efa77d4`, not from PyPI** [High]. PyPI's `veckit 0.1.2` (2026-08-25) differs from GitHub main:
   - It adds `severity_ratio`, `scale_ratio` and `size_fidelity`, which its comments call "what ranks" (open PR #6).
   - It lacks the issue-#7 right-handed-frame fix in `_canonicalise`, which main merged on 2026-10-02 and which affects `occupancy_dice` and sliced-Wasserstein.
   - The live leaderboard API still ranks `severity_slope` and `scale_log_ratio`.
   - *Note for the other plans:* `architecture.md` lists "veckit 0.1.2"; this should be the git pin.
6. **Board aggregation is known exactly** [High]. I fitted the public leaderboard API's per-metric skills (x100) to the headline scores:
   - T1 and T3 weights reproduce to ≤1e-4 points (203 and 204 rows).
   - The **T2 structure gate is `shape_group × min(1, 2·skill(neighborhood_mmd))`**, verified to ≤1e-4 on all 608 T2 rows. The metrics agent independently found the same formula.
   - The emulator in §4.4 implements this.
7. **A safety hazard in the challenge API** [High]. The site's own API guide says `format_only` "defaults to false": a bare upload is an **official** submission.
   - In P3 the server chooses the test split, and each board allows only **2 official submissions for the whole phase**.
   - Tooling must default to format-only and always send an `Idempotency-Key` on scored uploads. A client that does this is in §4.5 and was tested in dry-run mode.
8. **Scorer cost at server scale is small** [High, measured on DGX CPU, 16 threads, synthetic data at the exact board shapes].
   - A T1 board (truth 1,249 × ref 1,706 × prediction 2,000-5,118 cells × 32,285 genes) takes about **19 s** at 4.6-5.0 GB peak RSS: 6 s for the frozen probe plus 12.5 s for the panel.
   - Without the 10% subsample (full ~17K-cell stages), the same panel takes about 2 minutes and peaks at **48.6 GB RSS**.
   - The frozen probe, which only feeds the unranked composition_JSD, is the slowest step on T2/T3 (40-99 s). The emulator skips it by default.
9. **Raw counts are recoverable exactly from the released log1p data** [High on the public 150-cell samples; to re-check on the full files].
   - Each cell's size factor is `s = 1e4 / expm1(min nonzero)`, then `c = expm1(x)·s/1e4`.
   - Recovered counts are integral to within 1e-4 for 100% of nonzeros, with a median of 10.8K (E8.5) and 11.5K (E9.5) UMIs per cell.
   - This enables count-likelihood models (scVI NB/ZINB) on T1. The MERFISH files already ship `layers["counts"]`.
   - Samples store `.X` as **CSC** with about 8% **explicit stored zeros**. Convert to CSR and run `eliminate_zeros()` once at ingest.
10. **Compute is not the binding constraint** [Medium]. Every GPU model in scope (PCA, scVI, latent or full-gene flow matching, Sinkhorn up to 50K cells) needs seconds to tens of minutes on one H100.
    - The practical bottlenecks are CPU-side: exact minibatch OT is measured at 179 ms/step at batch 1024, versus 7 ms at batch 256, and CPU scoring takes about 20 s per board.
    - The bigger risk is GPU availability on the shared DGX.

---

## 1. Platform facts

| Item | Value | Source / confidence |
|---|---|---|
| DGX host | `p01r10srv05`, 8× NVIDIA H100 80GB HBM3 (SXM, compute capability 9.0), 192 CPU threads, 1 TB RAM (about 893 GB available), Ubuntu 22.04.5, glibc 2.35 | `nvidia-smi`, `/etc/os-release`, `free` [High] |
| Driver | 580.178.04, CUDA 13.0 | `nvidia-smi` [High] |
| uv on DGX | 0.12.13 at `~/.local/bin/uv` (**not on the non-login ssh PATH**; prefix `export PATH=$HOME/.local/bin:$PATH`) | [High] |
| uv locally | 0.11.7 (`/home/sports/.local/bin/uv`). Latest is 0.13.1 (2026-10-14). Local default `python3` is conda 3.12.11, hence `python-preference = "only-managed"` | [High] |
| Disk | `/home` 14 TB, 96% used, ~625 GB free | `df` [High] |
| H100 SXM peak | FP32 67 TFLOPS. TF32 989, BF16 1,979 and FP8 3,958 TFLOPS **with sparsity** (half that dense). 80 GB at 3.35 TB/s | NVIDIA H100 page [High] |

---

## 2. Environment (Q1)

### 2.1 Current versions (PyPI JSON, 2026-10-18)

| Package | Latest | Requires Python | Notes that matter here |
|---|---|---|---|
| numpy | 2.5.4 (2026-10-10) | ≥3.12 | Last 3.11-compatible line is 2.4.6. numba 0.68 caps numpy <2.6 |
| scipy | 1.18.1 | ≥3.12 | **The modern env resolves 1.16.3**, because pygam 0.12.0 (via cellrank) requires `scipy<1.17` |
| pandas | 3.0.6 (3.1.0rc0 out) | ≥3.11 | dynamo and spateo require `<3`. I pin `pandas>=2.2,<3` everywhere; resolved 2.3.3 |
| scikit-learn | 1.9.1 | ≥3.11 | |
| anndata | 0.13.4 | ≥3.12 | Needs numpy ≥2.1 and zarr ≥3.1. Extras `cu12`/`cu13` add cupy |
| scanpy | 1.12.4 | ≥3.12 | numpy ≥2. `[leiden]` adds igraph 1.0.0 and leidenalg 0.12.0 |
| h5py | 3.16.0 | ≥3.10 | |
| POT | 0.9.7.post1 | ≥3.7 | Imports every installed backend at import, **including JAX** (see §2.7) |
| torch | 2.14.1 (2026-09-30) | ≥3.10 | PyPI Linux wheel = CUDA 13.0 (`cuda-toolkit==13.0.3`, `nvidia-cudnn-cu13==9.24.0.43`, `nccl-cu13==2.30.7`, `triton~=3.8.0`) |
| torchcfm | 1.0.7 (2025-03) | none declared | Depends on torchdyn 1.0.6 → torchvision, `torchcde<0.3`, pytorch-lightning, ipywidgets, `poethepoet<0.11`. Only `torchcfm.utils` imports torchdyn |
| scvi-tools | 1.5.1 | ≥3.12 | No JAX dependency any more. Extras `cuda`/`cuda13` add cupy, torchvision and torchaudio (not needed) |
| moscot | 0.5.2 | ≥3.11 | Requires jax ≥0.6.1, ott-jax ≥0.6.0 (incl. `[neural]`), flax, diffrax and optax |
| ott-jax | 0.6.0 (2025-11-04) | ≥3.9 | Depends on **jaxopt ≥0.8** (unmaintained, import warning). **Online PointCloud broken on jax ≥0.9** (issue #700); fixed on main (#701, #696) but unreleased |
| jax / jaxlib | 0.11.2 | ≥3.12 | `jax[cuda13]` → `jax-cuda13-plugin[with-cuda]`; all its CUDA-13 lower bounds are met by torch-cu130's exact pins (verified in the lock) |
| cellrank | 2.3.3 | ≥3.12 | anndata ≥0.11, numpy ≥2. Brings pygam (scipy <1.17) and pygpcca 1.0.4 (`jinja2==3.0.3` pin) |
| dynamo-release | 1.5.3 (2026-03-17) | ≥3.9 | **`anndata<0.11`, `pandas<3`**. numpy-2 support since #744. Fails on statsmodels 0.15 |
| spateo-release | 1.1.1 on PyPI (2025-03); main `1bd8a35` (2026-02-25) | — | main pins **numpy `>=1.26,<2.0`, scipy `>=1.11,<1.14`, pandas `<3`, setuptools `<81`** for Python ≥3.11. Heavy deps: vtk, pyvista, opencv, geopandas, pysal, gpytorch, kornia |
| navigo (Navigo-release, "dynode_flow") | git main `ac21018` | — | `anndata>=0.8,<0.11`, POT ≥0.9, torch ≥2.0. Not on PyPI |
| veckit | 0.1.2 on PyPI ≠ git main `cc3859a` (0.1.1) | ≥3.9 | **Use git main** (§0.5). Installs top-level modules `common`, `T1`, `T2`, `T3`, `score_h5ad`, `veckit`, so never name your own top-level modules `metrics`, `common` or `T1`-`T3` |
| typer | 0.27.3 | ≥3.10 | |
| pydantic | 2.14.1 | ≥3.10 | |
| uv | 0.13.1 (DGX 0.12.13, local 0.11.7) | — | |

### 2.2 Incompatibility map and how each is resolved

| Conflict | Evidence | Resolution in the validated config | Confidence |
|---|---|---|---|
| numpy 1.x vs 2.x: spateo main needs `<2`; anndata 0.13, scanpy 1.12, cellrank, jax 0.11, zarr 3, geopandas 1.2 and opencv 5 need ≥2 | spateo `requirements.txt` (#360); PyPI metadata | `spatial` group in its own fork (numpy 1.26.4; opencv backtracks to 4.11.0.86, geopandas to 1.1.4, pysal to 25.7) | High |
| anndata `<0.11` (dynamo, Navigo) vs `>=0.11` (cellrank, mudata 0.4) | requirements files / PyPI | Forks: spatial gets anndata 0.10.9, modern gets 0.13.4. uv chose **mudata 0.3.10 in both** (works with scvi-tools 1.5.1: import OK) | High |
| pandas 3 vs `<3` | dynamo/spateo pins | `pandas>=2.2,<3` in base deps (2.3.3 everywhere, so h5ad interoperates) | High |
| scipy ≥1.18 vs pygam `<1.17` | pygam 0.12.0 metadata | Accept scipy 1.16.3 in the modern env | High |
| statsmodels 0.15 breaks dynamo | `ImportError: cannot import name 'multipletests'` reproduced on DGX | `statsmodels>=0.14.5,<0.15` in `spatial` (0.14.6) | High |
| ott-jax 0.6.0 online mode vs jax 0.11 | `AttributeError: module 'jax.interpreters.batching' has no attribute 'is_vmappable'` reproduced; fixed by ott #701 (merged 2026-06-08) | `ott-jax` from git `36338c9…` (0.6.1.dev30); online Sinkhorn verified | High |
| torchcfm → torchdyn baggage | torchdyn 1.0.6 metadata; torchcfm source imports | `exclude-dependencies = ["torchdyn"]`; `import torchcfm` and the OT-CFM sampler verified | High |
| pygpcca `jinja2==3.0.3` (old, has CVEs; conflicts with Flask/mlflow) | PyPI metadata; GPCCA/pyGPCCA main no longer pins it | `override-dependencies = ["jinja2>=3.1.6"]`; `import cellrank` OK | High |
| spateo git source needs git-lfs | `git-lfs filter-process: git-lfs: not found` on DGX | Archive URL of the same commit. **Caveat:** `spateo/tools/database/{human,mouse}_GRN.csv` arrive as LFS pointer stubs (only spateo's GRN tools use them) | High |
| torch CUDA variant: cu128 has no 2.14.x | download.pytorch.org index | cu130 (or cu126 if JAX cuda12 were ever needed) | High |
| torch cu130 + JAX in one env | lock shows nvidia-cudnn-cu13 9.24.0.43, nvidia-cublas 13.1.1.3 and nccl-cu13 2.30.7 shared by both | Works. `jax-cuda13-plugin` pulls `nvidia-cuda-nvcc` 13.4.92 (newer than the driver's 13.0); no ptxas/driver warning appeared in the tests | High (works) / Medium (heavier XLA kernels untested) |
| setuptools ≥81 drops `pkg_resources` (spateo uses it) | spateo pin `<81` | Resolved 80.10.2 | High |
| uv version skew (local 0.11.7, DGX 0.12.13) | `required-version` guard fired locally: "Required uv version `>=0.12` does not match … 0.11.7" | Run `uv self update` locally. Don't upgrade the DGX binary unilaterally: the `prachh` account is shared | High |

### 2.3 Why one project with a conflicting fork (and the fallback)

- **One lock, two environments.**
  - The main env (`.venv`) holds the `core`, `gpu`, `dev` and optionally `jax`/`traj` groups.
  - The spatial env (`.venv-spatial`, selected via `UV_PROJECT_ENVIRONMENT`) holds `spatial` and `dev`.
  - Both install the team package `vec`, so CLI and data-I/O code is shared.
  - Exchange data between the two envs as h5ad with plain dtypes (categorical obs, float32 X), or as `.npz`/parquet. Never pickle across them.
- **Validated:** the lock resolves 310 packages in 10.7 s on the first cold resolve, 31.9 s when a new git source has to be fetched, and about 1 s warm. The three syncs ran in 5-12 s from a warm DGX uv cache. Venv sizes: main (cu130 + jax + traj) about 6.5-7.4 GB, spatial about 7.3 GB on its own, CPU about 1.8 GB. These are mostly hardlinks into the uv cache.
- **Fallback:** if the spatial fork ever blocks a lock update, move `spatial` into a separate project at `envs/spatial/pyproject.toml` with its own lock. The cost is two locks to maintain.

### 2.4 The validated `pyproject.toml`

File: `SCRATCH/probe/pyproject.toml`. Lock: `SCRATCH/probe/uv.lock`, sha256 `8a8f1261dca4e86c481f0bbe8350ed6019fa08f3e46c3cbd1f1dc3dac1e1b032`. Both are also on the DGX at `~/v/vec2026-env-probe/`. `SCRATCH` is `SCRATCH`.

```toml
[build-system]
requires = ["uv_build>=0.12,<0.14"]
build-backend = "uv_build"

[project]
name = "vec"
version = "0.1.0"
description = "VEC 2026 Human-track team code: models, board emulator, CLI"
requires-python = "==3.12.*"
# Shared by EVERY environment (modern main env and the numpy<2 spatial env): keep these loose.
dependencies = [
  "numpy>=1.26",
  "scipy>=1.11",
  "pandas>=2.2,<3",            # one pandas major everywhere: h5ad round-trips between envs, no pandas-3 string dtype surprises
  "scikit-learn>=1.4",
  "anndata>=0.10.8",
  "h5py>=3.10",
  "typer>=0.12",
  "pydantic>=2.6",
  "pyyaml>=6.0",
  "httpx>=0.27",               # challenge API client (format checks, submissions, leaderboard snapshots)
  "filelock>=3.13",            # append-safe run registry / submission ledger
  "veckit",                    # pinned to the git commit the live board matches (see sources)
]

[project.optional-dependencies]
# Pick exactly one: --extra cpu (laptop / shared box) or --extra cu130 (DGX, driver 580 = CUDA 13.0).
cpu = ["torch==2.14.*"]
cu130 = ["torch==2.14.*"]

[project.scripts]
vec = "vec.cli:app"

[dependency-groups]
core = ["numpy>=2.1", "anndata>=0.12", "scanpy[leiden]>=1.11"]
gpu = ["pot>=0.9.6", "torchcfm>=1.0.7", "scvi-tools>=1.4", "torchdiffeq>=0.2.5", "geomloss>=0.3.1"]
jax = ["jax[cuda13]>=0.11", "ott-jax>=0.6", "moscot>=0.5.2"]
traj = ["cellrank>=2.3"]
spatial = [
  "spateo-release",
  "dynamo-release>=1.5.3",
  "navigo",
  # dynamo 1.5.3 imports statsmodels.sandbox.stats.multicomp.multipletests, removed in statsmodels 0.15
  "statsmodels>=0.14.5,<0.15",
]
dev = ["pytest>=8", "ruff>=0.6", "ipykernel>=6.29"]

[tool.uv]
required-version = ">=0.12"
python-preference = "only-managed"
environments = ["sys_platform == 'linux' and platform_machine == 'x86_64'"]
default-groups = ["core", "gpu", "dev"]
conflicts = [
  [{ extra = "cpu" }, { extra = "cu130" }],
  [{ group = "spatial" }, { group = "core" }],
  [{ group = "spatial" }, { group = "gpu" }],
  [{ group = "spatial" }, { group = "jax" }],
  [{ group = "spatial" }, { group = "traj" }],
]
# torchcfm's core API (ConditionalFlowMatcher, OTPlanSampler) never imports torchdyn; torchdyn drags in
# torchvision, torchcde<0.3, pytorch-lightning, ipywidgets and poethepoet<0.11.
exclude-dependencies = ["torchdyn"]
# pygpcca 1.0.4 (cellrank dep, last PyPI release 2022) pins jinja2==3.0.3; its code never imports jinja2.
override-dependencies = ["jinja2>=3.1.6"]

[tool.uv.sources]
torch = [
  { index = "pytorch-cpu", extra = "cpu" },
  { index = "pytorch-cu130", extra = "cu130" },
]
# PyPI veckit 0.1.2 is NOT the live scorer (it ranks severity_ratio/size_fidelity and lacks the issue-#7 frame fix).
veckit = { git = "https://github.com/aristoteleo/veckit", rev = "cc3859ab2a2305efc4544815d3f44f660efa77d4" }
# git source needs git-lfs (two GRN csvs are LFS-tracked); the archive tarball of the same commit does not.
spateo-release = { url = "https://github.com/aristoteleo/spateo-release/archive/1bd8a35e6e7b0dfc712d76b1c32987dc06652774.tar.gz" }
navigo = { git = "https://github.com/aristoteleo/Navigo-release", rev = "ac210180e565900b7b396db2a9deab8d8457e910" }
# ott-jax 0.6.0 (PyPI, 2025-11) breaks online PointCloud on jax>=0.9 (ott issue #700); fixed on main by #701/#696.
ott-jax = { git = "https://github.com/ott-jax/ott", rev = "36338c94707002772f2392ffa0a3acf11ab23857" }

[[tool.uv.index]]
name = "pytorch-cpu"
url = "https://download.pytorch.org/whl/cpu"
explicit = true

[[tool.uv.index]]
name = "pytorch-cu130"
url = "https://download.pytorch.org/whl/cu130"
explicit = true
```

(The file on disk also carries short `[tool.ruff]` and `[tool.pytest.ini_options]` blocks. `.python-version` contains `3.12`; uv installs a managed CPython 3.12.12.)

Resolved versions in the lock [High]:

| Package | Modern fork | Spatial fork |
|---|---|---|
| numpy | 2.5.4 | 1.26.4 |
| scipy | 1.16.3 | 1.13.1 |
| anndata | 0.13.4 | 0.10.9 |
| pandas | 2.3.3 | 2.3.3 |
| statsmodels | 0.15.0 | 0.14.6 |
| torch | 2.14.1+cu130 / 2.14.1+cpu, or 2.14.1 from PyPI if no extra is given | same |
| Same in both forks | scanpy 1.12.4, scikit-learn 1.9.1, POT 0.9.7.post1, torchcfm 1.0.7, scvi-tools 1.5.1, lightning 2.6.6, jax 0.11.2 (+cuda13 plugin), ott-jax 0.6.1.dev30, moscot 0.5.2, cellrank 2.3.3, spateo 1.1.1 (main archive), dynamo 1.5.3, navigo 0.1.0, veckit 0.1.1@cc3859a | |

### 2.5 Sync and run commands

```bash
# ---- local workstation (CPU only; editing, unit tests, small checks) ----
uv self update                                   # local uv is 0.11.7; the project requires >=0.12
uv sync --locked --extra cpu                     # core + gpu + dev groups, CPU torch (~1.8 GB venv)
uv run --locked --extra cpu vec info

# ---- DGX main environment ----
export PATH=$HOME/.local/bin:$PATH               # uv is not on the non-login PATH there
cd ~/v/vec2026                                   # team repo (copy pyproject.toml + uv.lock from the probe)
uv sync --locked --extra cu130 --group jax --group traj
# ---- DGX spatial environment (spateo / dynamo / Navigo; numpy<2 fork of the same lock) ----
UV_PROJECT_ENVIRONMENT=.venv-spatial uv sync --locked --no-default-groups --group spatial --group dev --extra cu130

# ---- running jobs (tmux, pinned GPU) ----
export UV_NO_SYNC=1                              # never let a bare `uv run` re-sync (see below)
CUDA_VISIBLE_DEVICES=3 tmux new -d -s vec-t1 '.venv/bin/vec train configs/t1/fm.yaml 2>&1 | tee logs/vec-t1.log'

# ---- freezing the environment for the method report ----
uv export --locked --extra cu130 --group jax --group traj --no-hashes > reports/requirements-main.txt
```

**Operational gotchas, all reproduced with `--dry-run` on the DGX** [High]:
- `uv sync` is exact. Re-running it **without** `--group jax --group traj` uninstalls 36 packages, including jax and cellrank.
- Forgetting `--extra cpu` swaps torch+cpu for the PyPI CUDA-13 torch plus about 20 NVIDIA wheels.
- `uv run` syncs inexactly, but still installs whatever the flags imply.
- Put the per-machine flags in one `scripts/sync.sh` that picks `cu130` if `nvidia-smi` works and `cpu` otherwise. Use `UV_NO_SYNC=1` (verified to leave the env untouched) or `.venv/bin/...` in job scripts.

### 2.6 Smoke test that passed on the DGX [High]

Script: `SCRATCH/probe/scripts/check_env.py`.
- **Main env:** numpy, anndata, scanpy, h5py, sklearn, POT, torch, torchcfm, scvi, lightning, veckit, cellrank, jax, ott, moscot, geomloss, torchdiffeq, httpx, filelock and vec all import. Also verified:
  - `torch.cuda` on an H100 with a 1024² matmul
  - `ot.sinkhorn` with CUDA tensors
  - the torchcfm `ExactOptimalTransportConditionalFlowMatcher` sample
  - `jax.devices()` returning `CudaDevice(id=0)`
  - ott Sinkhorn converging on GPU, in both dense and online `batch_size=4096` mode
- **Spatial env:** dynamo (9.8-10.9 s import), spateo, `spateo.align`, `spateo.tdr`, `spateo.tools`, navigo and veckit all import. torch-cu130 sees the GPU there too.

### 2.7 Runtime settings to put in `env.sh` (and at the top of the CLI)

```bash
export POT_BACKEND_DISABLE_JAX=1             # `import ot` otherwise initialises JAX on the GPU (verified: 5.5 s import + device_put)
export XLA_PYTHON_CLIENT_PREALLOCATE=false   # JAX default grabs 75% of the GPU: fatal on a shared H100
export CUBLAS_WORKSPACE_CONFIG=:4096:8       # needed for torch deterministic algorithms
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16   # don't spawn 192 BLAS threads
export UV_NO_SYNC=1
export VE_API_BASE=https://api.virtualembryo.ai   # VE_CHALLENGE_API_KEY comes from a secret file, never git
```

The POT behaviour is confirmed in POT's source: `ot/backend.py` reads `POT_BACKEND_DISABLE_{PYTORCH,JAX,CUPY,TENSORFLOW}` [High].

If XLA ever prints "the NVIDIA driver's CUDA version is 13.0 which is older than the PTX compiler version", add `constraint-dependencies = ["nvidia-cuda-nvcc<13.1"]` [Low; not observed].

### 2.8 Data-format facts that shape the pipeline

From the public veckit sample h5ads; re-check on the full release.

- **T1 samples** (150 × 32,285):
  - `.X` is **CSC** float32, uncompressed, with **~8% explicit stored zeros** (49,511 of 639,243 entries in the E8.5 sample).
  - `obsm['X_umap.harmony.rna']`; no `layers`.
  - `expm1(X).sum(1)` = 10,000.00 for every cell.
  - Genes per cell after `eliminate_zeros()`: 3,862 (E8.5) and 4,016 (E9.5).
  - [High for the samples; Medium that the full files match]
- **Count recovery works exactly:** see §0.9. Script: `SCRATCH/probe/scripts/inspect_h5ad.py` [High on samples].
- **MERFISH samples** (500 genes):
  - `layers['counts']` and `layers['log1p']` (plus `normalized_counts` in the KO sample).
  - `obsm['spatial_2D']` and `obsm['spatial_3D']`.
  - Extra obs such as `cm_celltype`, `slide`, `section`, `z_coord` and `time`.
  - About 35-43 detected genes and 69-102 counts per cell.
  - The Mab21l2 sample stores `.X` dense, the others CSC.
- **Panel identity:** index.json's 16-hex `genes_sha256` = `sha256("\n".join(genes))[:16]` (no trailing newline). It matched all 5 boards [High], so validators can check `var_names` in one line.
- **What to do at ingest:**
  - `X = X.tocsr(); X.eliminate_zeros()`
  - add `layers['counts']` (recovered) for T1
  - cache as CSR h5ad or zarr
  - keep raw files read-only with a sha256 manifest
  - CSC row-slicing (cell subsampling) is slow, and backed-mode row access on CSC is very slow.

---

## 3. Compute and memory (Q2)

### 3.1 Data sizes

File sizes come from the Data page. Cell counts for unreleased or other stages are inferred from index.json: the working copy is about 10% of the stage, and truth is half of the target's working copy.

| File | Size | Cells | Stored nonzeros (est.) | Sparse in RAM | Dense float32 |
|---|---|---|---|---|---|
| T1 E8.5_RNA | 571 MB | 16,787 | ≈71M (≈4.2K/cell incl. explicit zeros) | ≈0.57 GB | 2.17 GB |
| T1 E9.5_RNA | 590 MB | 17,057 | ≈74M | ≈0.59 GB | 2.20 GB |
| T1 E10.5 (P3) | 647 MB | ≈24,980 (truth 1,249 → 10% copy ≈ 2,498) [Medium] | ≈81M | ≈0.65 GB | 3.23 GB |
| T2 heart E10.5 (P3) | [redacted 2026-10-19: A] | [redacted 2026-10-19: A] [Medium] | 500 genes | — | [redacted 2026-10-19: A] |
| T2 heart E9.5 / E8.25 / E8.75 / E8.5 | 164 / 225 / 84 / [redacted 2026-10-19: A] MB | ≈53.7K / ≈58.7K / ≈24.8K / [redacted 2026-10-19: A] [Medium] | | | ≤0.12 GB each |
| T2 embryo E7.5 / E7.25 | not published | [redacted 2026-10-19: A] / ≈13.3K [Medium] | | | tiny |
| T3 Gata4 KO (2 reps), Mab21l2, WT E8.75/E9.5 | 374+409, 458, 84/164 MB | ≈46K per scored replicate (truth 2,320) [Low-Medium] | | | <0.1 GB |

- **Whole-download budget:** under 6 GB including the P3 releases. All of it should stay on the DGX.
- **Memory:** all T1 training stages fit densely in under 8 GB, and therefore entirely in one H100's memory.
- **Measured on a synthetic E9.5-sized stage** (CSR, 9.5% density, 52M nnz) [High]:
  - h5ad file: 420 MB. Extrapolating, the real 590 MB file implies about 13% stored density.
  - full read: 0.2 s (warm page cache)
  - densify: 0.9 s, giving 2.2 GB
  - `backed='r'` open: 0.03 s; taking 10% of rows: 0.11 s

### 3.2 PCA

Measured on DGX CPU, 16 threads, 17K cells [High]:

| Operation | Time |
|---|---|
| `scanpy.pp.pca` (50 PCs) on sparse, all 32,285 genes | 36 s |
| HVG selection (3,000, seurat flavour) | 2.7 s |
| PCA (50 PCs) on 3K HVGs | 3.1 s |
| sklearn randomized PCA on dense 32K genes | 9.0 s (5.5 GB peak) |

On GPU [Low-Medium, estimate]: `torch.pca_lowrank(q=50, niter=4)` on 34K × 32K dense float32 (4.4 GB) makes about 10 memory-bound passes. At 3.35 TB/s that is tens of milliseconds, so under 1 s including overheads. PCA is never the bottleneck.

### 3.3 scVI on an H100

- T1 training data: E8.5 + E9.5 = 33,844 cells, rising to about 58.8K with E10.5 at P3.
- scvi-tools' default `max_epochs = min(round(20000/n·400), 400)` gives **236 epochs** at 33.8K cells and about 136 at 58.8K. At batch 128 that is about 265 minibatches per epoch, or about 62K steps.
- The default model (1×128 hidden, 10 latent) is tiny. Steps are dataloader-bound at roughly 3-10 ms, giving **about 3-10 minutes per model** with under 5 GB of GPU memory.
- Using all 32K genes instead of 2-5K HVGs widens only the decoder output layer, which is negligible on an H100.
- [Low-Medium: estimate. I found no published H100 timing for this size; measure it once with the default settings.]
- **Prerequisite:** raw counts, recoverable exactly per §0.9.

### 3.4 OT / Sinkhorn on 10-50K cells

Dense cost or kernel matrix memory is n² × 4 B: 10K → 0.4 GB, 20K → 1.6 GB, 50K → 10 GB, 100K → 40 GB.

| Backend | How it scales | Estimated time, 1,000 iterations, d = 50 | Notes |
|---|---|---|---|
| POT `ot.sinkhorn` (torch backend, GPU, dense M) | 2 reads of n×m per iteration, memory-bound | 20K: ≈1-2 s; 50K: ≈6-10 s [Low-Medium] | Verified working on CUDA [High]. Small ε needs `method="sinkhorn_log"` (≈3-5× slower: ≈20-50 s at 50K) and more than the default `numItermax=1000`. POT warned "did not converge" at ε = 0.05 with defaults |
| ott-jax dense `PointCloud` | same as POT, plus jit | similar | Verified on GPU [High] |
| ott-jax **online** `PointCloud(batch_size=…)` | recomputes cost tiles: ≈4·n²·d FLOPs per iteration, O(n·d) memory | 50K: ≈10-20 s; scales to 100K+ [Low-Medium] | **Requires ott-jax git main with jax 0.11** (verified fixed) |
| GeomLoss `SamplesLoss("sinkhorn")` | tensorized backend = dense; `online` needs KeOps | ε-scaling converges in tens of iterations | pykeops 2.3 is sdist-only with JIT compilation, untested on this pip-only CUDA 13 stack; not installed [Medium] |
| moscot (wraps ott) | temporal and spatiotemporal problems | minutes at 34K-60K cells with `batch_size` [Low] | |
| Exact EMD (`ot.emd`, CPU) | network simplex, about O(n³) | **measured** (single CPU thread, d = 50): batch 128: 1.4 ms; 256: 7.0 ms; 512: 34 ms; 1024: 179 ms; 2048: 668 ms [High] | CPU Sinkhorn (200 iterations): 4.2 / 8.2 / 40 / 146 / 707 ms at the same sizes |

To replace the estimates with numbers, run `SCRATCH/probe/scripts/bench_gpu.py` on a free GPU (it was smoke-tested on CPU):

```bash
CUDA_VISIBLE_DEVICES=<free> XLA_PYTHON_CLIENT_PREALLOCATE=false POT_BACKEND_DISABLE_JAX=1 .venv/bin/python -I scripts/bench_gpu.py --sizes 10000,20000,50000
```

### 3.5 Flow-matching training

- **Latent model** (30-50 PCs or a scVI latent):
  - A 3×512 MLP has about 0.6M parameters; one step at batch 1024 is well under 1 ms of GPU compute.
  - 20-50K steps take about 1-3 minutes **if the coupling is not the bottleneck**.
  - With torchcfm's `ExactOptimalTransportConditionalFlowMatcher`, `ot.emd` runs on the CPU every step: 179 ms at batch 1024, so 20K steps would take about an hour.
  - Use batch ≤256-512 (7-34 ms), a GPU Sinkhorn coupling, or a precomputed global coupling sampled per batch. The organisers' `neural_ode` used a global entropic-OT coupling. [High for the OT timings; Medium for the GPU step estimate]
- **Full-gene model** (32,285 → 1K → 1K → 32,285):
  - About 67M parameters; ≈6 × 67M × 512 ≈ 2e11 FLOPs per step, so about 1-3 ms in bf16 on an H100.
  - 50K steps take about 1-3 minutes. GPU memory is a few GB (Adam states ≈0.8 GB). [Low-Medium]
- **Sampling** 5K cells with 100 Euler steps: under 1-2 s.
- **Throughput:** one full cycle (cached prep, train, sample, 5-seed emulator scoring) is about 10-15 minutes, so ≥50 configurations per GPU-day. CPU scoring parallelises across the 192 threads.

### 3.6 Scorer (veckit) runtime per board

Measured [High]: DGX CPU, 16 threads, veckit `cc3859a`, numpy 2.5.4, scipy 1.16.3, sklearn 1.9.1. Synthetic log1p data at the exact board shapes; the timings are what matter, not the scores. Raw JSON: `SCRATCH/probe/bench/*.json`.

| Board shape (truth / ref / pred cells × genes) | Frozen probe | Full panel `score_task*_v2` | Main components | Peak RSS |
|---|---|---|---|---|
| T1 server scale: 1,249 / 1,706 / 2,000 × 32,285 | 6.3 s | 12.4 s | DE (Mann-Whitney + BH) 5.5 s, MMD 3.9 s, variogram 1.2 s, energy distance (3 BLAS pairwise calls on 1,500 × 1,500 × 32,285) 0.55 s | 4.6 GB |
| T1 server scale, pred 5,118 (the max) | 6.2 s | 12.6 s | same | 5.0 GB |
| T1 **without** the 10% subsample: 17,057 / 16,787 / 16,787 | 25 s | 96 s | DE 72 s, MMD 10.9 s | **48.6 GB** |
| T2 heart extrap: [redacted 2026-10-19: A] / 5,374 / 8,000 × 500 + coords | **98.7 s** | 3.3 s | neighborhood_mmd 1.0 s, d2/dice/SW < 0.2 s | 0.8 GB |
| T2 embryo interp: [redacted 2026-10-19: A] / 1,330 / 1,500 × 498 | 5.1 s | 5.2 s | neighborhood_mmd 2.1 s | 0.6 GB |
| T3 Gata4: 2,320 / 2,483 / 2,500 × 500 | 40.4 s | 5.5 s | | 0.7 GB |

What this means:
- Emulate at server scale (10% subsample). Full-stage scoring needs a large-memory node and about 2 minutes per board.
- The frozen probe, `LogisticRegression(max_iter=300)` on the truth, feeds only the **unranked** composition_JSD. Skip it, or train it once per board and cache it. Its time on real data may differ from synthetic.
- Truth-only work can be cached per board and seed: DE sets via `de_genes(true, ref)`, the PCA fit on the truth, and the median-heuristic bandwidth when truth ≤2,000 cells. With that, per-candidate T1 scoring drops to about 2-5 s.
- Exact parity needs care. veckit draws the truth subsample (and the variogram gene pairs) *after* the prediction's draws from the same RNG, so those draws depend on the prediction's cell count. Either call veckit as-is, or replay the RNG consumption.

**DGX GPU contention at 2026-10-18 18:08 IST** [High, observed]:
- GPUs 0, 5 and 6 each had about 65 GB in use by Qwen3.5-9B SFT jobs.
- GPUs 1, 2 and 4 had about 25 GB in use by Qwen3.5-2B SFT jobs plus MD jobs.
- GPUs 3 and 7 each ran one alchemy_ddg MD job (about 0.6 GB, 92-94% utilisation).
- VEC jobs need under 10 GB, so they can co-locate on GPUs 3 or 7 memory-wise, but they will share SMs (roughly 2× slower, and they slow the MD jobs). **Agree 1-2 dedicated GPUs with the user** for 2026-10-20 to 2026-12-02.

---

## 4. Reproducibility and experiment management (Q3)

### 4.1 Repository layout (private repo; data never in git)

```
vec2026/
  pyproject.toml  uv.lock  .python-version  env.sh  scripts/sync.sh
  src/vec/        cli.py  emulator.py  validate.py  api.py  registry.py  config.py  data.py  models/...
  configs/        t1/*.yaml  t2/*.yaml  t3/*.yaml  boards/*.yaml   (board specs: task, ref stage, target stage, panel, frac, seeds)
  data/panels/    index.json + *.genes.txt (re-pull at P3; verify genes_sha256)
  data/EXTERNAL_SOURCES.yaml   (disclosure: dataset, URL/DOI, licence, stages kept/removed, filter code path)
  runs/registry.jsonl          (git) ; runs/<run_id>/... (DGX only: configs, env.json, metrics, preds, ckpts)
  submissions/ledger.jsonl     (git) ; submissions/files/<sha256>.h5ad (DGX + one off-box copy)
  leaderboard/<UTC>_<board>.json  (daily public snapshots)
  reports/        method-card drafts, uv export snapshots
```

A starter implementation is in `SCRATCH/probe/src/vec/`: `emulator.py`, `validate.py`, `api.py`, `registry.py` and `cli.py`. All of them were exercised on the DGX (tests in §4.4-4.6).

### 4.2 Configs, seeds, determinism

- **Configs:** pydantic v2 models loaded from YAML, with `--set key=value` overrides. The *resolved* config is written into the run directory.
  - Hydra is unnecessary for a small team; typer plus pydantic is explicit and type-checked.
- **Seeds:** record a separate seed for each of data subsampling (emulator), model init, batch order and sampling.
  - `seed_everything()` (in `registry.py`) seeds `random`, numpy's legacy global state, `torch.manual_seed` and `torch.cuda.manual_seed_all`.
  - It sets `cudnn.benchmark=False` and `torch.use_deterministic_algorithms(True, warn_only=True)`.
  - JAX: pass explicit `PRNGKey`s. Optionally `XLA_FLAGS=--xla_gpu_deterministic_ops=true` for final runs.
  - Use numpy `Generator`s in your own code.
- **What reproducible means here:** same config and seeds give the same scores within noise. Bitwise GPU determinism is not guaranteed.
  - The **submitted file itself** is the artefact of record. Keep it, plus its sha256, which the server also returns in the submission detail.
  - Before each official submission, re-run the final config once and compare emulator scores.

### 4.3 Run registry

- **Run directory:** `runs/<YYYYMMDD-HHMMSS>_<task>_<slug>_<gitsha7>/` containing:
  - `config.yaml` (resolved)
  - `env.json`: git SHA, dirty flag, `uv.lock` sha256, Python and package versions, host, `CUDA_VISIBLE_DEVICES`. Written by `registry.environment_fingerprint()`.
  - `metrics/<board>.json`: raw metrics, skills and score per seed, plus mean ± sd
  - `pred.h5ad` (or its path), `ckpt/`, `log.txt`
- **Registry index:** `runs/registry.jsonl` gets one line per scored (run, board) with the summary and pointers.
  - Appends go through a `filelock`, so concurrent tmux jobs are safe. The file is git-tracked so the team shares one view.
  - Query it with pandas, or DuckDB if you like.
- **Avoid MLflow and W&B as the source of truth.**
  - The MLflow server (Flask) conflicts with cellrank's jinja2 pin unless overridden.
  - W&B is cloud-hosted: never upload h5ads or derived data there. The data is unpublished, and the API guide says the download integration "does not authorize redistribution or hosting copies of the data". Scalar metrics only, if at all.

### 4.4 Local board emulator (`src/vec/emulator.py`, tested)

**Pipeline, mirroring the documented server process:**
1. Load the full target and reference stages. The reference is the preceding stage (T1/T2) or the matched WT (T3).
2. Subsample each to **10%** with `round(0.1·n)` cells. This matches index.json: E9.5's 17,057 cells give `ref_cells` 1,706; `int()` would give 1,705 [Medium].
3. Split the target with `veckit.common.core_metrics.split_half(n, seed)`: half A is the truth, half B is the ceiling prediction.
4. The floor is the reference working copy verbatim, including coordinates (copy_last), or the WT (wt_identity).
5. Compute raw metrics for prediction, floor and ceiling with veckit's own `score_task{1,2,3}_v2` at `seed=0`.
6. Compute skill per metric with `veckit.skill`, on |m| for `severity_slope` and `scale_log_ratio`. NaN maps to 0.
7. Aggregate with the exact weights, including the T2 gate:

```python
WEIGHTS = {"T1": [(.25, {"de_score": 1}), (.25, {"de_direction": 1}), (.30, {"mmd_u": 1}), (.20, {"variogram": 1})],
           "T2": [(.25, {"de_score": .5, "de_direction": .5}), (.25, {"mmd_u": .6, "variogram": .4}),
                  (.25, "SHAPE_GATED"), (.25, {"neighborhood_mmd": 1})],
           "T3": [(.30, {"de_score": 1}), (.25, {"de_direction": 1}), (.25, {"severity_slope": 1}),
                  (.20, {"mmd_u": .6, "variogram": .4})]}
# SHAPE_GATED = mean(skill d2_shape, occupancy_dice, |scale_log_ratio|) * min(1, 2 * skill(neighborhood_mmd))
```

**Verification:**
- Fitting the public API rows (200+ per board) reproduces every headline to ≤1e-4. The emulator's `aggregate()` gives 67.428 for a public T2 heart-extrap row whose headline is 67.4261.
- End-to-end on the public 150-cell samples (with `frac=1.0`, since they are too small to subsample): T1 copy_last = **50.00**, T3 wt_identity = **50.00**, and a crude shift "model" = 71.4 ± 1.2 [High].
- The T2 sample floor scored 54.17. Its half-vs-half "ceiling" for occupancy_dice (0.350) was *worse* than the floor (0.362), and `veckit.skill` then returns 1.0. **The emulator must flag any metric whose ceiling does not beat its floor**, since this will happen on small pseudo-boards [High].

**Usage policy:**
- **Seeds:** run 5 subsample/split seeds and report mean ± sd. The server uses one unknown fixed draw. Treat gaps under about 2 sd as ties.
- **Pseudo-boards now (P2):**
  - T1: E8.5 → E9.5 (weak: only one input stage)
  - T2 heart, extrapolation: {E8.25, E8.75} → E9.5 (reference E8.75)
  - T2 heart, interpolation: {E8.25, E9.5} → E8.75 (reference E8.25)
  - T2 embryo, interpolation: {E6.75, E8.0} → E7.25 (reference E6.75)
  - T3: none yet, since only one KO is released.
- **From P3 (2026-10-20), exact validation boards:**
  - T1 {E8.5, E9.5} → E10.5
  - T2 heart E9.5 → E10.5 and heart interpolation E8.5
  - T2 embryo E7.5
  - T3 Mab21l2 → Gata4, scoring **both replicates**, since the server publishes both and the scorer fixes which one ranks
  - plus leave-one-transition-out variants
- **Calibration at P3:**
  - Try plausible subsampling implementations and seeds, and check whether the published `index.json` floor/ceiling anchors reproduce to 3-4 digits. If they do, the emulator matches the server's working copy exactly.
  - Independently, re-implement the published `pseudobulk_shift` (T1 = 48.5751, with per-metric skills 46.22 / 53.90 / 52.55 / 38.91), `shift_transfer` (T3 = 42.4843) and `gene_ko` (42.5033). Confirm the emulator lands within noise of them.
- **Rules boundary:** emulate only from released data, including validation answers once they are released. **Never** back out target properties from returned server scores (Rules s.10). Scores may select between predictions, never compute them.
- **Performance:** at server scale the full T1 panel is about 12.5 s and the probe is skipped by default (`_NoProbe`). Add per-board caching of the truth-side work (§3.6) and run seeds in parallel with joblib (for example 4 seeds × 4 BLAS threads).

### 4.5 Submission safety and the challenge API

All [High]: from `https://api.virtualembryo.ai/openapi.json` and the site's API guide (account → API keys page; text read from its JS bundle).
- **Auth:** `Authorization: Bearer $VE_CHALLENGE_API_KEY`.
  - Only the captain creates keys, at most 10 active per team.
  - Keys expire within 90 days and no later than the deadline; scopes cover submissions and quota.
- **Data:**
  - `GET /challenge/data/manifest` lists the released files and their keys.
  - `GET /challenge/data/link?key=…` returns a storage URL valid for **≤15 minutes**. Fetch that URL *without* the API key. Downloads do not consume quota.
  - Do this on the DGX.
- **Upload:** `POST /challenge/submissions` (multipart) with `task`, `model`, `setting`, `mode`, `target`, `architecture`, `format_only` and `file`.
  - The site maps board keys to fields as follows: `T1:*` gives `task=T1`; `T2:heart:val_extrap` gives `setting=heart, mode=extrap`; `T3:gata4` gives `target=gata4`.
  - `architecture` is the name to use consistently across tasks for the Generality Award.
- **`format_only` defaults to false.** For a format check, send `format_only=true` and **omit** the Idempotency-Key.
  - Format checks are unlimited and spend no attempts, but every upload counts against a daily *upload* (storage) limit.
- **Idempotency:** 1-128 characters from `[A-Za-z0-9._:-]`.
  - An identical retry returns the original response without spending another attempt, even after the quota is exhausted.
  - The same key with a different payload returns 409.
  - Back off 2/4/8 s on 5xx; do not retry 4xx responses unchanged.
- **Responses:** a rejected upload can still return HTTP 200, so read `status` and `problems`.
  - Status values: queued, scoring, scored, rejected (`reject_reason`), cancelled, eligibility_review.
  - Poll `GET /challenge/submissions/{id}` every 15-30 s. Responses include the server's `sha256` of the uploaded bytes.
- **Quota:** `GET /challenge/integration/quota` reports daily used/limit and, in P3, `official_per_board`.
  - In P3 "the server chooses the phase's validation or test split", so **every scored P3 upload is a test attempt**.
  - Heart interpolation is validation-only and is omitted in the final phase.
  - Final boards: T1 test (E12.5), T2 heart extrapolation (E12.5), T2 embryo interpolation (E7.75), T3 β-catenin.
- **Public endpoints (no key):** `/challenge/phase`, `/challenge/boards`, `/challenge/leaderboard?task=&board=&track=` (per-metric skills × 100 for 200 rows plus reference rows) and `/challenge/overall`.

**Client behaviour** (`src/vec/api.py` with `vec submit`, tested in dry-run mode):
- Default is `format_only=true`.
- An official upload requires `--scored --confirm-board <board>` plus a local format check.
- The ledger guard refuses a third official P3 upload on any board.
- The idempotency key is `<board>.<sha256[:40]>` of the exact file bytes, so a retry can never spend a second attempt.
- Every upload is appended to `submissions/ledger.jsonl` with the file sha256 and the run id.

**Two-shot P3 protocol:**
1. Freeze candidates on emulated validation and pseudo-boards (mean ± sd).
2. Run the free API format check on the exact final bytes.
3. Make the first official submission early, by about 2026-11-20, as an anchor.
4. Make the second by 2026-11-30, never in the last hour (Rules s.12: outage remedies are discretionary).
5. Tag git as `sub/<board>/<n>` and keep the file, checkpoint, config and `uv export` snapshot.

### 4.6 typer CLI

Implemented in `src/vec/cli.py`. `vec --help` responds in 0.27 s because heavy imports happen inside commands.
- `vec info`: environment fingerprint.
- `vec validate PRED --board T1:val`: offline mirror of the server's format check.
  - Checks: panel identity and order plus `genes_sha256`, cell range, finite and non-negative values, `spatial_3D` (n, ≥3), ≤1,200 MiB, and median `expm1(X).sum(1)` ≈ 1e4 (the server never renormalises).
  - Tested: it rejected the 150-cell samples for the cell range while accepting their panels.
- `vec phase` and `vec lb-snapshot`: public API. Tested; the 2026-10-18 tops were T1 76.0, heart interpolation 75.4, heart extrapolation 67.3 (human), embryo 77.7 and T3 80.3 (human).
- `vec submit … [--scored --confirm-board B] [--dry-run]` and `vec status ID`.
- Still to add: `vec data pull` (manifest → links → sha256-verified raw), `vec data prep` (CSC → CSR, drop zeros, recover counts, HVG/PCA caches), `vec board build`, `vec train CONFIG`, `vec predict RUN --board`, `vec score PRED --board --seeds 0-4` (emulator → registry), `vec registry ls/compare` and `vec report RUN`.

### 4.7 Kit for the 14-day method report (start it now)

Generate the following automatically from the ledger and registry for every official submission:
1. Data provenance: released files with sha256, plus every external source with its stage window and the removal code (`data/EXTERNAL_SOURCES.yaml`).
   - Watch for contamination: **any pretrained weights (for example Navigo checkpoints or scvi-hub models) trained on atlases spanning E9.5-E13.5, E7.25-E8.0, or the Gata4/β-catenin KOs are disallowed.** Retrain from allowed data instead.
2. Preprocessing: count recovery, HVGs, PCA and any spatial alignment.
3. Model and training config, seeds, and a "what was trained on what" stage table.
4. Compute used (GPU hours) and wall-clock time.
5. The `uv.lock` hash and `uv export`.
6. Emulator scores (mean ± sd) next to the server's returned score.
7. Exact commands to reproduce the run.

The organisers may ask to repeat a run, so keep checkpoints on the DGX plus one off-box copy (code and models only, never the data).

---

## 5. Risks and open questions

1. **Scorer drift.** veckit PR #6 (`severity_ratio`, `size_fidelity`) is open and already shipped in PyPI 0.1.2.
   - Re-check `index.json` metric names and the leaderboard's metric keys at P3 start and weekly.
   - Keep both forms computable in the emulator (0.1.2's code shows the formulas) [Medium probability of change, unknown timing].
2. **Unknown server subsampling seed and procedure, and which half is truth.** Calibrate at P3 (§4.4). Until then, emulator numbers are expectations with a spread [Medium].
3. **Full release files may differ from the samples** (CSC, explicit zeros, count integrality). Verify at download, which is cheap [Medium].
4. **GPU availability on the shared DGX:** all 8 GPUs were busy at survey time.
5. **JAX on CUDA 13 with nvcc 13.4 versus driver 13.0:** not stress-tested beyond small ops [Low risk].
6. **ott-jax depends on unmaintained jaxopt,** and moscot's API may break with future jax. Keep jax pinned by the lock; upgrade deliberately.
7. **Spateo's GRN databases are LFS stubs in the archive install.** If T3 needs `mouse_GRN.csv`, fetch it from GitHub LFS separately (needs git-lfs or the GitHub media URL) [Medium].
8. **The frozen-probe runtime on real T2/T3 data is unmeasured** (synthetic: 40-99 s). It is unranked, so skip or cache it.

---

## 6. Artifacts produced

`SCRATCH` = `SCRATCH`.

| Path | What it is |
|---|---|
| `SCRATCH/probe/pyproject.toml`, `SCRATCH/probe/uv.lock` | Validated environment (lock sha256 `8a8f1261…b032`) |
| `SCRATCH/probe/src/vec/{cli,api,validate,registry,emulator}.py` | Tested starter tooling |
| `SCRATCH/probe/scripts/check_env.py` | Environment smoke test |
| `SCRATCH/probe/scripts/bench_veckit.py`, `SCRATCH/probe/bench/*.json` | Scorer benchmark and its results |
| `SCRATCH/probe/scripts/bench_gpu.py` | GPU benchmark to run when a GPU is free |
| `SCRATCH/probe/scripts/inspect_h5ad.py` | Layout inspection and count-recovery check |
| `SCRATCH/probe/scripts/lock_summary.py` | Lock version summary |
| `SCRATCH/probe/scripts/test_emulator.py` | Emulator smoke test |
| `SCRATCH/probe/data/panels/` | index.json and the panel files |
| DGX `~/v/vec2026-env-probe/` | Same files plus synced `.venv` (cu130 + jax + traj), `.venv-spatial` and `.venv-cpu`. Delete it after copying into the team repo; the venvs are mostly hardlinks into the uv cache |

---

## 7. References

**Challenge**
- Challenge site and pages (Data, Evaluation, Baselines, FAQ, Rules): https://virtualembryo.ai/challenge (local copies in `SCRATCH/site_text/`)
- Challenge API OpenAPI spec: https://api.virtualembryo.ai/openapi.json
- Phase: https://api.virtualembryo.ai/challenge/phase
- Boards: https://api.virtualembryo.ai/challenge/boards
- Leaderboard: https://api.virtualembryo.ai/challenge/leaderboard?task=T2&board=T2:heart:val_extrap
- Board index and panels: https://virtualembryo.ai/challenge/panels/index.json

**veckit (official scorer)**
- Repository, pinned commit `cc3859ab2a2305efc4544815d3f44f660efa77d4`: https://github.com/aristoteleo/veckit
- Issues and PRs: #2 (memory-safe scorer: before ≈14.5 GiB / 2:24, after ≈6.5 GiB / 3:37), #6 (PSS/TSR change), #7 and #9 (right-handed frame), #8 (MMD memory), #12 (no renormalisation): https://github.com/aristoteleo/veckit/issues
- PyPI veckit 0.1.2: https://pypi.org/project/veckit/0.1.2/

**Aristoteleo stack**
- spateo-release, main `1bd8a35` and requirements PR #360: https://github.com/aristoteleo/spateo-release/pull/360
- dynamo-release: requirements, plus #727 (temporary anndata pin), #744 (NumPy 2): https://github.com/aristoteleo/dynamo-release
- Navigo-release, main `ac21018`: https://github.com/aristoteleo/Navigo-release

**OT and flow matching**
- ott-jax issue #700, PR #701 (`is_vmappable`), PR #696 (JAX 0.10 compatibility): https://github.com/ott-jax/ott/issues/700 , https://github.com/ott-jax/ott/pull/701 , https://github.com/ott-jax/ott/pull/696
- torchcfm (conditional-flow-matching): https://github.com/atong01/conditional-flow-matching
- pyGPCCA (no jinja2 pin on main): https://github.com/GPCCA/pyGPCCA
- GeomLoss / Sinkhorn divergences: Feydy et al., AISTATS 2019, https://arxiv.org/abs/1810.08278
- FlashSinkhorn (possible fast Sinkhorn option, not evaluated): https://arxiv.org/abs/2602.03067 ; code https://github.com/ot-triton-lab/flash-sinkhorn

**Packaging and PyTorch**
- PyPI JSON API, `https://pypi.org/pypi/<pkg>/json`, for numpy, scipy, pandas, scikit-learn, anndata, scanpy, h5py, POT, torch, torchcfm, torchdyn, scvi-tools, moscot, ott-jax, jaxopt, jax, jaxlib, jax-cuda12/13-plugin, cellrank, pygam, pygpcca, dynamo-release, spateo-release, mudata, statsmodels, numba, zarr, typer, pydantic, uv, cuda-toolkit 13.0.3 (all fetched 2026-10-18)
- PyTorch wheel indexes: https://download.pytorch.org/whl/cu130/torch/ , https://download.pytorch.org/whl/cu128/torch/ , https://download.pytorch.org/whl/cu126/torch/ (PEP 658 `.metadata` sidecars read for the 2.14.1 pins)
- uv + PyTorch guide: https://docs.astral.sh/uv/guides/integration/pytorch/
- uv settings (`conflicts`, `exclude-dependencies`, `override-dependencies`, `required-version`, `environments`): https://docs.astral.sh/uv/reference/settings/
- uv syncing semantics (exact `sync`, inexact `run`): https://docs.astral.sh/uv/concepts/projects/sync/

**Hardware and scVI**
- NVIDIA H100 specifications: https://www.nvidia.com/en-us/data-center/h100/
- scvi-tools max-epochs heuristic: https://docs.scvi-tools.org/en/1.4.3/api/reference/scvi.model.get_max_epochs_heuristic.html ; discussion https://discourse.scverse.org/t/motivation-for-max-epochs-np-min-round-20000-n-cells-400-400/718
