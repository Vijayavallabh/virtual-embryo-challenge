"""GEO access for the Task 1 outside atlases: the family SOFT parser, the sample-folder file check and the one
download function. Written for this folder (the T1 code never imports the T3 outside-knockout code; a test of that
code checks it), in the same form as that fetcher. Requests go only to ftp.ncbi.nlm.nih.gov over https with a generic
user agent and nothing that names a team, a person or a challenge; no redirect is followed (one refuses), and a 403,
429, 5xx or timeout is retried a bounded number of times with a pause (TRIES, PAUSE_S).
"""

from __future__ import annotations

import hashlib
import re
import time
import urllib.error
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
# The repository's standard field names that hold a slash, accepted after the block's prefix (Series_ or Sample_):
# a literal list, so every other odd key still refuses the series. The first real plan run (290da3f) met
# "Sample_contact_zip/postal_code" and refused the whole series on it.
SLASH_KEYS = frozenset({"contact_zip/postal_code"})
_FILE_NAME = re.compile(r"[A-Za-z0-9._+-]+")
# The fetch's bounded retry (the first real plan run met one HTTP 403 that was a 200 a moment later): at most TRIES
# requests per file, PAUSE_S seconds apart, retried only on HTTP 403, 429 and 5xx and on a timeout; any other error,
# and every redirect, refuses at once.
TRIES = 4
PAUSE_S = 20.0
RETRY_CODES = frozenset({403, 429})


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
            if not _KEY.fullmatch(key[len(prefix):]) and key[len(prefix):] not in SLASH_KEYS:
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


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    """Every redirect refuses: the fetch stays on the one host over https and follows no redirect at all (the
    reading that refuses more of 'follows no redirect off it'; a same-host redirect refuses too)."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise Refused(f"{req.full_url}: the repository answered with a redirect ({code} to {newurl!r}): not followed")


def opener() -> urllib.request.OpenerDirector:
    return urllib.request.build_opener(_NoRedirect())


def _retryable(e: BaseException) -> bool:
    if isinstance(e, urllib.error.HTTPError):
        return e.code in RETRY_CODES or 500 <= e.code <= 599
    if isinstance(e, TimeoutError):
        return True
    return isinstance(e, urllib.error.URLError) and isinstance(e.reason, TimeoutError)


def get(url: str, dst: Path, max_bytes: int, open_url=None, sleep=time.sleep, tries: int = TRIES,
        pause_s: float = PAUSE_S) -> tuple[int, str]:
    """The one network call: stream an allowed URL to dst (outside the repository), refusing past max_bytes. Up to
    `tries` requests, `pause_s` seconds apart, on HTTP 403, 429 and 5xx and on a timeout (each try starts the file
    again); any other error, a redirect (refused, never followed) and the byte limit end it at once. `open_url(request,
    timeout)` is the opener's open (a fake in the tests)."""
    url = https_url(url)
    dst = cm.outside_repo(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    part = dst.with_name(dst.name + ".part")
    open_url = open_url or opener().open
    errors = []
    for attempt in range(1, tries + 1):
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        n, h = 0, hashlib.sha256()
        try:
            with open_url(req, timeout=120) as r, open(part, "wb") as f:
                final = getattr(r, "geturl", lambda: url)()
                if https_url(final) != url:
                    raise Refused(f"{url}: the response came from {final!r}, not the requested address")
                while block := r.read(1 << 20):
                    n += len(block)
                    if n > max_bytes:
                        f.close()
                        part.unlink(missing_ok=True)
                        raise Refused(f"{url}: over the byte limit ({max_bytes}); nothing kept")
                    h.update(block)
                    f.write(block)
        except Refused:
            part.unlink(missing_ok=True)
            raise
        except (urllib.error.URLError, TimeoutError) as e:
            part.unlink(missing_ok=True)
            if not _retryable(e):
                raise Refused(f"{url}: {e!r}: not retried") from None
            errors.append(f"try {attempt}: {e!r}")
            if attempt < tries:
                sleep(pause_s)
            continue
        part.replace(dst)
        return n, h.hexdigest()
    raise Refused(f"{url}: no answer after {tries} tries ({'; '.join(errors)})")
