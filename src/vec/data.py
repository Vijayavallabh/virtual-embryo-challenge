"""Competition data download.

The site signs in through GitHub and keeps a short-lived bearer token in the browser (localStorage key
`ve-challenge-session`). The API lists released files at /challenge/data/manifest and hands out a signed URL per file at
/challenge/data/link?key=..., valid for 15 minutes. In the final phase the P2 validation answers have their own state,
list and links under /challenge/data/validation/ (state, manifest, link?key=...), served once the team has unlocked
them. The unlock is a POST that permanently ends P2 submissions for the whole team; a person makes it on the site's
Data page (done on 2026-10-28), and this module never sends it.

The token is a team API key that can also make official submissions, so this module only ever issues GET requests for
the lists and download links; uploads belong in a dedicated client that defaults to format-only (feat-016). The token
stays on the machine that runs this module (default ~/.config/vec/token, mode 600) and is never printed: load_token
refuses one that is not a single run of printable ASCII (http.client would quote such a header value in an error), and
the server's error text has it redacted before it goes into an exception.
`fetch` asks for one signed link at a time, just before that file is downloaded, and passes it to curl on its standard
input so it never appears in a process list or a log.
"""

from __future__ import annotations

import http.client
import json
import os
import re
import stat
import subprocess
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

API_BASE = "https://api.virtualembryo.ai"
VALIDATION = "validation"  # DataFile.source of the P2 validation answers (final phase, after the team's unlock)
DEFAULT_TOKEN_FILE = Path.home() / ".config" / "vec" / "token"
SAFE_KEY = re.compile(r"^[A-Za-z0-9._/+-]+$")
TOKEN_CHARS = re.compile(r"[\x21-\x7e]+")  # printable ASCII without spaces


class AuthError(RuntimeError):
    pass


def load_token(path: Path = DEFAULT_TOKEN_FILE) -> str:
    env = os.environ.get("VEC_TOKEN")
    if env:
        return _valid_token(env.strip(), "VEC_TOKEN")
    path = Path(path).expanduser()
    if not path.exists():
        raise AuthError(f"no token at {path}; see docs/DATA_ACCESS.md")
    if path.stat().st_mode & (stat.S_IRWXG | stat.S_IRWXO):
        raise AuthError(f"{path} is readable by others; run: chmod 600 {path}")
    token = path.read_text().strip()
    if token.startswith("{"):  # the whole localStorage value was pasted
        try:
            token = json.loads(token)["token"]
        except (ValueError, KeyError, TypeError):
            raise AuthError(f'{path} holds JSON without a "token" entry') from None
    return _valid_token(token, str(path))


def _valid_token(token, source: str) -> str:
    """The token, if it is one run of printable ASCII; the error never repeats it."""
    if not isinstance(token, str) or not TOKEN_CHARS.fullmatch(token):
        raise AuthError(f"the token in {source} is empty or has whitespace, control or non-ASCII characters inside it; "
                        "it must be a single line of printable characters")
    return token


def redact(value, token: str | None):
    """`value` with the token replaced by <token>: in a string, or in the keys and values of parsed JSON at any depth.
    The escaped form a repr shows is replaced too, since http.client quotes a bad header value in its ValueError."""
    if not token:
        return value
    if isinstance(value, str):
        for form in sorted({token, repr(token)[1:-1]}, key=len, reverse=True):
            value = value.replace(form, "<token>")
        return value
    if isinstance(value, dict):
        return {redact(k, token): redact(v, token) for k, v in value.items()}
    if isinstance(value, list):
        return [redact(v, token) for v in value]
    return value


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    """Never follow a redirect: urllib would resend the Authorization header to the new host."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise urllib.error.HTTPError(req.full_url, code, f"redirect to {newurl} refused (token not forwarded)",
                                     headers, fp)


def api_get(path: str, token: str) -> dict:
    """GET an API path with the token. Errors carry the server's text with the token redacted (before it is cut), and
    never chain the original exception, whose message may hold it."""
    req = urllib.request.Request(f"{API_BASE}{path}", headers={"Authorization": f"Bearer {token}"})
    try:
        with urllib.request.build_opener(_NoRedirect).open(req, timeout=60) as r:
            return json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        detail = redact(e.read().decode(errors="replace") or str(e.msg), token)[:200]
        if e.code in (401, 403):
            raise AuthError(f"HTTP {e.code} from {path}: {detail} (token expired or registration "
                            "incomplete?)") from None
        raise RuntimeError(f"HTTP {e.code} from {path}: {detail}") from None
    except (ValueError, http.client.InvalidURL) as e:  # a bad header value or URL, or an answer that is not JSON
        raise RuntimeError(f"request to {path} failed: {redact(f'{type(e).__name__}: {e}', token)}") from None


@dataclass
class DataFile:
    key: str
    task: str
    role: str
    stage: str
    modality: str
    bytes: int  # 0 when the list gives no size
    source: str = "train"  # or VALIDATION: which list the file came from, and so which link endpoint serves it
    extra: dict = field(default_factory=dict)  # any other fields the list gives (a checksum, say), as served

    @property
    def relpath(self) -> str:
        if not SAFE_KEY.match(self.key) or ".." in self.key.split("/"):
            raise ValueError(f"refusing unsafe file key {self.key!r}")
        return self.key


KNOWN_FIELDS = {"key", "task", "role", "stage", "modality", "bytes"}


def _files(m: dict, source: str) -> list[DataFile]:
    return [DataFile(key=f["key"], task=str(f.get("task", "")), role=str(f.get("role", "")),
                     stage=str(f.get("stage", "")), modality=str(f.get("modality", "")),
                     bytes=int(f.get("bytes") or 0), source=source,
                     extra={k: v for k, v in f.items() if k not in KNOWN_FIELDS}) for f in m["files"]]


def manifest(token: str) -> list[DataFile]:
    return _files(api_get("/challenge/data/manifest", token), "train")


def validation_state(token: str) -> dict:
    """The team's state for the P2 validation release, as served (the site reads `unlocked`, `available` and
    `lock_scope`). Asking changes nothing."""
    return api_get("/challenge/data/validation/state", token)


def validation_manifest(token: str) -> list[DataFile]:
    return _files(api_get("/challenge/data/validation/manifest", token), VALIDATION)


def released_files(token: str) -> tuple[list[DataFile], dict]:
    """The training files, then the validation answers if the team has unlocked them, and the validation state. The
    site lists the validation files only when the state says unlocked and available, and so does this."""
    files, state = manifest(token), validation_state(token)
    if state.get("unlocked") and state.get("available"):
        val = validation_manifest(token)
        clash = {f.key for f in files} & {f.key for f in val}
        if clash:  # one path on disk for two files: never let one overwrite the other
            raise RuntimeError(f"keys listed both as training and validation files: {sorted(clash)}")
        files += val
    return files, state


def signed_link(key: str, token: str, source: str = "train") -> str:
    base = "/challenge/data/validation/link" if source == VALIDATION else "/challenge/data/link"
    return api_get(f"{base}?key={urllib.parse.quote(key)}", token)["url"]


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
    """Download one file to `dest` (a local dir, or a dir on `host` over ssh). Skips files already complete: of the
    listed size, or, when the list gives no size, present at all (curl writes to `.part` and the file gets its name
    only after curl finished without error)."""
    target = f"{resolve_dest(dest, remote=bool(host)).rstrip('/')}/{f.relpath}"
    have = _remote_size(host, target) if host else (Path(target).stat().st_size if Path(target).exists() else -1)
    if have >= 0 and (have == f.bytes or not f.bytes):
        return "present"
    url = signed_link(f.key, token, f.source)
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
