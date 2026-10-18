"""Competition data download.

The site signs in through GitHub and keeps a short-lived bearer token in the browser (localStorage key
`ve-challenge-session`). The API lists released files at /challenge/data/manifest and hands out a signed URL per file at
/challenge/data/link?key=..., valid for 15 minutes.

The token stays on the machine that runs this module (default ~/.config/vec/token, mode 600) and is never printed.
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


def _remote_size(host: str, remote_path: str) -> int:
    out = subprocess.run(["ssh", "-o", "BatchMode=yes", host, f"stat -c %s '{remote_path}' 2>/dev/null || echo -1"],
                         capture_output=True, text=True, check=True).stdout.strip()
    return int(out or -1)


def fetch_one(f: DataFile, token: str, dest: str, host: str | None) -> str:
    """Download one file to `dest` (a local dir, or a dir on `host` over ssh). Skips files already complete."""
    target = f"{dest.rstrip('/')}/{f.relpath}"
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
