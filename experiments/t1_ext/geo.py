"""GEO access for the Task 1 outside atlases: the family SOFT parser, the sample-folder file check and the one
download function. Written for this folder (the T1 code never imports the T3 outside-knockout code; a test of that
code checks it), in the same form as that fetcher. Requests go only to ftp.ncbi.nlm.nih.gov over https with a generic
user agent and nothing that names a team, a person or a challenge.
"""

from __future__ import annotations

import hashlib
import re
import urllib.parse
import urllib.request
from pathlib import Path

import t1x_common as cm

Refused = cm.Refused
HOSTS = ("ftp.ncbi.nlm.nih.gov",)
USER_AGENT = "python-urllib (public metadata and data download)"
# never fetched, whatever an atlas's patterns say
DENY = re.compile(r"fastq|\.bam$|\.bai$|fragments|atac|peaks?[._]|\.bw$|bigwig|\.bed(\.gz)?$|cloupe|\.sra$|"
                  r"_raw\.tar$|\.tar(\.gz)?$|\.zip$|\.rds$|\.rdata$|\.h5seurat$|\.h5ad$|\.loom$", re.IGNORECASE)
OTHER = "_other_lines"
LINE_SEPARATORS = "\r\x0b\x0c\x1c\x1d\x1e\x85  "
_KEY = re.compile(r"[A-Za-z0-9_]+")
_FILE_NAME = re.compile(r"[A-Za-z0-9._+-]+")


def soft_url(series: str) -> str:
    if not re.fullmatch(r"GSE\d+", series):
        raise Refused(f"{series}: only GEO series are planned")
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series[:-3]}nnn/{series}/soft/{series}_family.soft.gz"


def https_url(url: str) -> str:
    u = url.strip()
    if u.startswith("ftp://"):
        u = "https://" + u[len("ftp://"):]
    p = urllib.parse.urlparse(u)
    if p.scheme != "https" or p.hostname not in HOSTS:
        raise Refused(f"{url}: not an allowed public repository URL")
    return u


def decode(raw: bytes, series: str) -> str:
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError as e:
        raise Refused(f"{series}: undecodable bytes in the family SOFT (byte {e.start}): refused") from None


def parse_soft(text: str) -> dict:
    """{'series': {key: [values]}, 'samples': {GSM: {key: [values]}}} from a family SOFT text, values verbatim
    (stripped). Lines are split on '\\n' only; a line holding another line separator, a block accession that is not a
    plain GSE or GSM, or a field name that is not a plain key refuses the series."""
    series: dict[str, list[str]] = {}
    samples: dict[str, dict[str, list[str]]] = {}
    cur, kind = None, None
    for i, line in enumerate(text.split("\n"), 1):
        if line.endswith("\r"):
            line = line[:-1]
        if any(c in LINE_SEPARATORS for c in line):
            raise Refused(f"a line separator inside SOFT line {i}: refused")
        if line.startswith("^"):
            k, _, ident = line[1:].partition("=")
            kind, cur = k.strip().upper(), ident.strip()
            if (kind == "SAMPLE" and not re.fullmatch(r"GSM\d+", cur)) or (kind == "SERIES" and
                                                                             not re.fullmatch(r"GSE\d+", cur)):
                raise Refused(f"a {kind} block whose accession is not plain (line {i}): refused")
            if kind == "SAMPLE":
                samples.setdefault(cur, {})
            continue
        if not line.strip() or kind not in ("SAMPLE", "SERIES"):
            continue
        rec = samples[cur] if kind == "SAMPLE" else series
        prefix = "Sample_" if kind == "SAMPLE" else "Series_"
        key, eq, val = line[1:].partition("=")
        key, val = key.strip(), val.strip()
        if line.startswith("!") and eq and key.startswith(prefix):
            if not _KEY.fullmatch(key[len(prefix):]) and key != "Series_contact_zip/postal_code":
                raise Refused(f"{cur}: a field name that is not a plain key (line {i}): refused")
            rec.setdefault(key[len(prefix):], []).append(val)
        else:
            rec.setdefault(OTHER, []).append(line.strip())
    return {"series": series, "samples": samples}


def own_file(gsm: str, url: str) -> tuple[str | None, str]:
    """(https url, '') when the link is a file in the sample's OWN GEO folder, named after the sample; else (None,
    reason)."""
    if not re.fullmatch(r"GSM\d+", gsm):
        return None, "not a GEO sample accession"
    try:
        u = https_url(url)
    except Refused:
        return None, "not an https or ftp link on the repository host"
    p = urllib.parse.urlparse(u)
    if p.netloc != HOSTS[0] or p.params or p.query or p.fragment:
        return None, "not a plain link on the repository host"
    name = p.path.rsplit("/", 1)[-1]
    if p.path != f"/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{name}":
        return None, "not under the sample's own GEO folder"
    if not _FILE_NAME.fullmatch(name) or ".." in name or not name.startswith(f"{gsm}_"):
        return None, "file name not plain or not named after the sample"
    return u, ""


def get(url: str, dst: Path, max_bytes: int) -> tuple[int, str]:
    """The one network call: stream an allowed URL to dst (outside the repository), refusing past max_bytes."""
    url = https_url(url)
    dst = cm.outside_repo(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    part = dst.with_name(dst.name + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    n, h = 0, hashlib.sha256()
    with urllib.request.urlopen(req, timeout=120) as r, open(part, "wb") as f:
        while block := r.read(1 << 20):
            n += len(block)
            if n > max_bytes:
                f.close()
                part.unlink(missing_ok=True)
                raise Refused(f"{url}: over the byte limit ({max_bytes}); nothing kept")
            h.update(block)
            f.write(block)
    part.replace(dst)
    return n, h.hexdigest()
