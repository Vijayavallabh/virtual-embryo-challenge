"""Competition data download.

The site signs in through GitHub and keeps a short-lived bearer token in the browser (localStorage key
`ve-challenge-session`). The API lists released files at /challenge/data/manifest and hands out a signed URL per file at
/challenge/data/link?key=..., valid for 15 minutes.

The token is a team API key that can also make official submissions, so this module only ever issues GET requests for
the manifest and download links; uploads belong in a dedicated client that defaults to format-only (feat-016). The token
stays on the machine that runs this module (default ~/.config/vec/token, mode 600) and is never printed.
`fetch` asks for one signed link at a time, just before that file is downloaded, and passes it to curl on its standard
input so it never appears in a process list or a log.
"""

from __future__ import annotations

import json
import os
import re
import stat
import subprocess
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

API_BASE = "https://api.virtualembryo.ai"
DEFAULT_TOKEN_FILE = Path.home() / ".config" / "vec" / "token"
SAFE_KEY = re.compile(r"^[A-Za-z0-9._/+-]+$")


class AuthError(RuntimeError):
    pass


def load_token(path: Path = DEFAULT_TOKEN_FILE) -> str:
    env = os.environ.get("VEC_TOKEN")
    if env:
        return env.strip()
    path = Path(path).expanduser()
    if not path.exists():
        raise AuthError(f"no token at {path}; see docs/DATA_ACCESS.md")
    if path.stat().st_mode & (stat.S_IRWXG | stat.S_IRWXO):
        raise AuthError(f"{path} is readable by others; run: chmod 600 {path}")
    token = path.read_text().strip()
    if token.startswith("{"):  # the whole localStorage value was pasted
        token = json.loads(token)["token"]
    return token


def api_get(path: str, token: str) -> dict:
    req = urllib.request.Request(f"{API_BASE}{path}", headers={"Authorization": f"Bearer {token}"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        detail = e.read().decode(errors="replace")[:200]
        if e.code in (401, 403):
            raise AuthError(f"HTTP {e.code} from {path}: {detail} (token expired or registration incomplete?)") from e
        raise RuntimeError(f"HTTP {e.code} from {path}: {detail}") from e


@dataclass
class DataFile:
    key: str
    task: str
    role: str
    stage: str
    modality: str
    bytes: int

    @property
    def relpath(self) -> str:
        if not SAFE_KEY.match(self.key) or ".." in self.key.split("/"):
            raise ValueError(f"refusing unsafe file key {self.key!r}")
        return self.key


def manifest(token: str) -> list[DataFile]:
    m = api_get("/challenge/data/manifest", token)
    return [DataFile(key=f["key"], task=str(f.get("task", "")), role=str(f.get("role", "")),
                     stage=str(f.get("stage", "")), modality=str(f.get("modality", "")),
                     bytes=int(f.get("bytes", 0))) for f in m["files"]]


def signed_link(key: str, token: str) -> str:
    return api_get(f"/challenge/data/link?key={urllib.parse.quote(key)}", token)["url"]


def inventory(path: Path, max_cells_for_checks: int = 20000, seed: int = 0) -> dict:
    """Summary of one h5ad: shape, storage, layers, metadata, normalisation check, cell types, coordinate scale."""
    import anndata as ad
    import numpy as np
    from scipy import sparse

    a = ad.read_h5ad(path, backed="r")
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(a.n_obs, min(a.n_obs, max_cells_for_checks), replace=False))
    sub = a[idx].to_memory()
    X = sub.X
    nnz_frac = (X.nnz / (X.shape[0] * X.shape[1])) if sparse.issparse(X) else float((np.asarray(X) != 0).mean())
    Xd = X.toarray() if sparse.issparse(X) else np.asarray(X)
    lib = np.expm1(np.clip(Xd, 0, 50)).sum(1)
    out = {
        "file": str(path),
        "n_obs": int(a.n_obs),
        "n_vars": int(a.n_vars),
        "x_storage": type(a.X).__name__ if not hasattr(a.X, "format") else f"{a.X.format}",
        "x_dtype": str(X.dtype),
        "nonzero_fraction": round(float(nnz_frac), 4),
        "x_min": float(Xd.min()),
        "x_max": float(Xd.max()),
        "expm1_rowsum_median": float(np.median(lib)),
        "expm1_rowsum_iqr": [float(np.percentile(lib, 25)), float(np.percentile(lib, 75))],
        "layers": [str(k) for k in a.layers.keys()],
        "obs_columns": list(a.obs.columns),
        "obsm_keys": [str(k) for k in a.obsm.keys()],
        "obs_names_unique": bool(a.obs_names.is_unique),
    }
    for col in ("celltype", "cm_celltype", "condition", "genotype", "sample_id", "slide"):
        if col in a.obs:
            vc = a.obs[col].astype(str).value_counts()
            out[f"{col}_n_levels"] = int(vc.size)
            out[f"{col}_top"] = {k: int(v) for k, v in vc.head(25).items()}
    if "spatial_3D" in a.obsm:
        C = np.asarray(a.obsm["spatial_3D"])[:, :3].astype(float)
        Xc = C - C.mean(0)
        z = np.unique(C[:, 2])
        out["spatial_3D"] = {
            "extent": [float(v) for v in (C.max(0) - C.min(0))],
            "rms_radius": float(np.sqrt((Xc ** 2).sum(1).mean())),
            "n_unique_z": int(z.size),
            "median_z_step": float(np.median(np.diff(z))) if z.size > 1 else None,
            "duplicate_coordinates": int(C.shape[0] - np.unique(C, axis=0).shape[0]),
        }
    return out


def _remote_size(host: str, remote_path: str) -> int:
    out = subprocess.run(["ssh", "-o", "BatchMode=yes", host, f"stat -c %s '{remote_path}' 2>/dev/null || echo -1"],
                         capture_output=True, text=True, check=True).stdout.strip()
    return int(out or -1)


def resolve_dest(dest: str, remote: bool) -> str:
    """Make a destination path safe to single-quote in a shell command.

    A quoted '~' is never expanded: on 2026-10-18 the first download landed in a directory literally named '~'. ssh
    commands start in the remote $HOME, so a home-relative path is written without the tilde; locally it is expanded.
    """
    if remote:
        return dest[2:] if dest.startswith("~/") else dest
    return os.path.expanduser(dest)


def fetch_one(f: DataFile, token: str, dest: str, host: str | None) -> str:
    """Download one file to `dest` (a local dir, or a dir on `host` over ssh). Skips files already complete."""
    target = f"{resolve_dest(dest, remote=bool(host)).rstrip('/')}/{f.relpath}"
    if host:
        if _remote_size(host, target) == f.bytes:
            return "present"
    elif Path(target).exists() and Path(target).stat().st_size == f.bytes:
        return "present"
    url = signed_link(f.key, token)
    curl_cfg = f'url = "{url}"\n'
    part = f"{target}.part"
    script = (f"mkdir -p \"$(dirname '{target}')\" && curl -sS -L --fail --retry 3 -C - -K - -o '{part}' "
              f"&& mv '{part}' '{target}'")
    cmd = ["ssh", "-o", "BatchMode=yes", host, script] if host else ["bash", "-c", script]
    subprocess.run(cmd, input=curl_cfg, text=True, check=True)
    size = _remote_size(host, target) if host else Path(target).stat().st_size
    if f.bytes and size != f.bytes:
        raise RuntimeError(f"{f.key}: got {size} bytes, manifest says {f.bytes}")
    return "downloaded"
