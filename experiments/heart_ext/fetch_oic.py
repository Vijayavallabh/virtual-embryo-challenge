#!/usr/bin/env python
"""Fetch the outside files of heart file OIC125 (row R70) by a pinned list, with resume. CPU and network only.

Fixed here, before any fetch (HEART_MOON1_challenge.md FIX 4 and FIX 6; register rows X021 and X022 of
EXTERNAL_DATA.md; X002 is refused for OIC125 and is not on this list, so OIC125 has two studies):

- **MOSTA bin size (X022).** The bin1 GEM files give one row per gene per DNB, at the chip's DNB pitch of 0.5 um
  (bin1 = one DNB; the authors' bin50 is 25 um). OIC125 bins them at **bin24: squares of 24 x 24 DNB, 12 um a side**
  (MOSTA_BIN_DNB, DNB_PITCH_UM), the Xenium grid's size, so both platforms' units have one physical size. Only the
  raw count column (MIDCounts or the file's UMI-count column) and the DNB coordinates are read; never a per-section
  MOSTA.h5ad, the all-stage object, a cell-bin file or any author label.
- **Xenium grid (X021).** The decoded-transcript table only, binned on a **12 um grid** (XENIUM_GRID_UM) of its
  x_location and y_location (microns). Never the vendor's cells, cell or nucleus boundaries, segmentation, cell-feature
  matrix, clusters or embeddings: of the table, only gene name, x, y, the decoding quality and the gene flag are read.

**The transcript table comes only inside an archive that holds refused files.** The dataset page offers it only inside
the output bundles: "Xenium Output Bundle (full)" (`..._outs.zip`, 83,505,883,137 bytes, md5 8b2a2da6...) and the
"Xenium Output Bundle (Xenium Explorer subset)" (`..._xe_outs.zip`, 30,910,928,431 bytes); both also hold the cells,
the boundaries, the cell-feature matrix and the analysis (clusters, embeddings). The smallest route, pinned here: read
the full bundle's zip directory with HTTP range requests (its end-of-directory records and the central directory, a
few MB), find the ONE member whose name ends in `transcripts.parquet`, and request that member's byte range alone. No
other member's bytes are requested. A stored member is written as it comes; a deflated one is kept as downloaded
(`<name>.deflate`) and inflated beside it; its CRC-32 and size must equal the directory's.

The pinned list (PINNED) holds the exact URLs found on the sources' own pages on 2026-11-25:
- the MOSTA download page https://db.cngb.org/stomics/mosta/download/ and its directory listing
  https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/Bin1_matrix/ (23 files: E15.5_E1S1-E1S4, E15.5_E2S1,
  E16.5_E1S1-E1S5, E16.5_E2S1-E2S13), each size the server's Content-Length on that day (48,844,857,975 bytes in all);
- the 10x dataset page "Whole Mouse Pup (FFPE) with 5K Mouse Pan Tissue and Pathways Panel and Snap25a/Snap25b
  Add-on" (xenium-prime-ffpe-neonatal-mouse), its bundle link on s3-us-west-2.amazonaws.com.

Behaviour:
- Any URL not on the list is refused (exit 2) before any request; so is a server size that differs from the pin.
- Each download goes to `<name>.part` in the output folder, resumed from its length with a Range request, and is
  renamed when complete. A finished file of the pinned size is not fetched again.
- `<out>/MANIFEST.json` records, per file: name, URL, size, sha256 (and, for the transcript member, the bundle's URL
  and size, the member's name, method, offsets, sizes and CRC-32), plus the constants above and this script's sha256.
- Nothing is ever deleted: no file is removed or truncated; a .part longer than its pin is refused.

  python3 -I -B experiments/heart_ext/fetch_oic.py --out ~/v/vec-ext/heart_oic            # every pinned file
  python3 -I -B experiments/heart_ext/fetch_oic.py --out D --only E16.5_E1S3_GEM_bin1.tsv.gz
  python3 -I -B experiments/heart_ext/fetch_oic.py --out D --url <a pinned URL>
  python3 -I -B experiments/heart_ext/fetch_oic.py --list                                 # no network

Requests carry a generic user agent and nothing that names a team, a person or a challenge. Exit 0 when every
selected file is complete and recorded, 2 on a refusal, 1 when a download fails after its retries.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import sys
import threading
import time
import urllib.error
import urllib.request
import zipfile
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# ---- constants fixed before any fetch -----------------------------------------------------------------------------
DNB_PITCH_UM = 0.5          # Stereo-seq DNB centre-to-centre distance (bin1 = one DNB)
MOSTA_BIN_DNB = 24          # OIC125's MOSTA unit: 24 x 24 DNB
MOSTA_BIN_UM = MOSTA_BIN_DNB * DNB_PITCH_UM   # = 12.0 um
XENIUM_GRID_UM = 12.0       # OIC125's Xenium unit: a 12 um square of the transcript table's x, y
assert MOSTA_BIN_UM == XENIUM_GRID_UM

USER_AGENT = "python-urllib (public data download)"
CHUNK = 8 << 20
RETRIES = 30
TIMEOUT = 120

MOSTA_BASE = "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/Bin1_matrix/"
MOSTA_SECTIONS = (
    ("E15.5_E1S1", 2224024877),
    ("E15.5_E1S2", 2001279746),
    ("E15.5_E1S3", 2884177092),
    ("E15.5_E1S4", 1785626749),
    ("E15.5_E2S1", 1199695752),
    ("E16.5_E1S1", 2805156097),
    ("E16.5_E1S2", 1845184454),
    ("E16.5_E1S3", 4053007636),
    ("E16.5_E1S4", 2369817047),
    ("E16.5_E1S5", 3715090384),
    ("E16.5_E2S1", 1151553811),
    ("E16.5_E2S2", 1552138129),
    ("E16.5_E2S3", 1537519288),
    ("E16.5_E2S4", 1503472533),
    ("E16.5_E2S5", 2170522623),
    ("E16.5_E2S6", 1743904353),
    ("E16.5_E2S7", 2464045197),
    ("E16.5_E2S8", 2448765042),
    ("E16.5_E2S9", 2873449795),
    ("E16.5_E2S10", 2100729384),
    ("E16.5_E2S11", 1725927936),
    ("E16.5_E2S12", 1332186011),
    ("E16.5_E2S13", 1357584039),
)
XENIUM_BUNDLE_URL = ("https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/3.0.0/Xenium_Prime_Mouse_Pup_FFPE/"
                     "Xenium_Prime_Mouse_Pup_FFPE_outs.zip")
XENIUM_BUNDLE_BYTES = 83505883137
XENIUM_BUNDLE_MD5 = "8b2a2da686a7ce369b434f9291d15d2e"   # the page's record of the whole bundle; not fetched whole
XENIUM_MEMBER = "transcripts.parquet"
XENIUM_NAME = "X021_Xenium_Prime_Mouse_Pup_FFPE_transcripts.parquet"


def _entries() -> list[dict]:
    out = []
    for section, size in MOSTA_SECTIONS:
        name = f"{section}_GEM_bin1.tsv.gz"
        out.append({"name": name, "url": MOSTA_BASE + name, "bytes": size, "kind": "file", "study": "X022",
                    "section": section, "stage": section.split("_")[0]})
    out.append({"name": XENIUM_NAME, "url": XENIUM_BUNDLE_URL, "bytes": None, "kind": "zip_member",
                "archive_bytes": XENIUM_BUNDLE_BYTES, "member": XENIUM_MEMBER, "study": "X021",
                "section": "pup_P1", "stage": "P1"})
    return out


PINNED = tuple(_entries())
MOSTA_TOTAL_BYTES = sum(s for _, s in MOSTA_SECTIONS)


class Refused(Exception):
    """A refusal: exit 2, nothing requested beyond what was needed to see it."""


class Failed(Exception):
    """A download that failed after its retries: exit 1; the .part is kept for the next call."""


# ---- HTTP ---------------------------------------------------------------------------------------------------------
class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: D401
        raise Refused(f"refused: {req.full_url} redirects to {newurl}; only the pinned URL is requested")


_OPENER = urllib.request.build_opener(_NoRedirect)


def _request(url: str, start: int | None = None, end: int | None = None):
    headers = {"User-Agent": USER_AGENT, "Accept-Encoding": "identity"}
    if start is not None:
        headers["Range"] = f"bytes={start}-" if end is None else f"bytes={start}-{end}"
    return _OPENER.open(urllib.request.Request(url, headers=headers), timeout=TIMEOUT)


def _total_from(resp, start: int | None) -> int:
    """The resource's total size from a response, checking that a range request got the range asked for."""
    if start is not None:
        if resp.status != 206:
            raise Failed(f"{resp.url}: asked for a range, got HTTP {resp.status}")
        cr = resp.headers.get("Content-Range", "")
        if not cr.startswith(f"bytes {start}-") or "/" not in cr:
            raise Failed(f"{resp.url}: unexpected Content-Range {cr!r}")
        return int(cr.rsplit("/", 1)[1])
    return int(resp.headers["Content-Length"])


def _range_bytes(url: str, start: int, end: int, total: int) -> bytes:
    """Bytes [start, end] inclusive of a resource of pinned total size (a directory read: small)."""
    last = None
    for attempt in range(RETRIES):
        try:
            with _request(url, start, end) as resp:
                got_total = _total_from(resp, start)
                if got_total != total:
                    raise Refused(f"refused: {url} is {got_total} bytes on the server, pinned {total}")
                data = resp.read()
            if len(data) != end - start + 1:
                raise Failed(f"{url}: short range read")
            return data
        except Refused:
            raise
        except (urllib.error.URLError, OSError, Failed, ValueError) as exc:
            last = exc
            time.sleep(min(60, 2 ** min(attempt, 6)))
    raise Failed(f"{url}: range {start}-{end} failed after {RETRIES} attempts: {last}")


def _download_range(url: str, part: Path, start: int, length: int, total: int, log) -> None:
    """Append bytes [start + have, start + length) of the resource to `part`, resuming from its length."""
    last = None
    for attempt in range(RETRIES):
        have = part.stat().st_size if part.exists() else 0
        if have > length:
            raise Refused(f"refused: {part} holds {have} bytes, more than the {length} pinned; nothing is deleted")
        if have == length:
            return
        try:
            with _request(url, start + have, start + length - 1) as resp:
                got_total = _total_from(resp, start + have)
                if got_total != total:
                    raise Refused(f"refused: {url} is {got_total} bytes on the server, pinned {total}")
                with open(part, "ab") as fh:
                    while True:
                        buf = resp.read(CHUNK)
                        if not buf:
                            break
                        fh.write(buf)
                        if fh.tell() > length:
                            raise Refused(f"refused: {part} grew past its pinned {length} bytes")
            last = None
        except Refused:
            raise
        except (urllib.error.URLError, OSError, Failed, ValueError) as exc:
            last = exc
            log(f"  retry {attempt + 1} for {part.name}: {exc}")
            time.sleep(min(60, 2 ** min(attempt, 6)))
    have = part.stat().st_size if part.exists() else 0
    if have != length:
        raise Failed(f"{url}: {have} of {length} bytes after {RETRIES} attempts: {last}")


def _sha256_crc(path: Path) -> tuple[str, int, int]:
    h = hashlib.sha256()
    crc = 0
    n = 0
    with open(path, "rb") as fh:
        while True:
            buf = fh.read(CHUNK)
            if not buf:
                break
            h.update(buf)
            crc = zlib.crc32(buf, crc)
            n += len(buf)
    return h.hexdigest(), crc & 0xFFFFFFFF, n


# ---- one whole file -----------------------------------------------------------------------------------------------
def fetch_file(entry: dict, out: Path, log=print) -> dict:
    final = out / entry["name"]
    part = out / (entry["name"] + ".part")
    size = entry["bytes"]
    if final.exists():
        if final.stat().st_size != size:
            raise Refused(f"refused: {final} exists with {final.stat().st_size} bytes, pinned {size}; "
                          "nothing is deleted")
    else:
        _download_range(entry["url"], part, 0, size, size, log)
        os.replace(part, final)
    sha, _crc, n = _sha256_crc(final)
    if n != size:
        raise Failed(f"{final}: {n} bytes, pinned {size}")
    return {"name": entry["name"], "url": entry["url"], "size": n, "sha256": sha, "study": entry["study"],
            "section": entry["section"], "stage": entry["stage"]}


# ---- one member of a remote zip -----------------------------------------------------------------------------------
class RangeFile(io.RawIOBase):
    """A read-only, seekable view of a remote file of pinned size, read with HTTP range requests (for zipfile)."""

    def __init__(self, url: str, total: int):
        super().__init__()
        self.url, self.total, self.pos = url, total, 0
        self.requests = 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        base = {0: 0, 1: self.pos, 2: self.total}[whence]
        self.pos = max(0, base + offset)
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.total - self.pos
        n = min(n, self.total - self.pos)
        if n <= 0:
            return b""
        self.requests += 1
        data = _range_bytes(self.url, self.pos, self.pos + n - 1, self.total)
        self.pos += len(data)
        return data

    def readinto(self, b):
        data = self.read(len(b))
        b[: len(data)] = data
        return len(data)


def member_info(url: str, total: int, member: str) -> dict:
    """The one member whose name is `member` or ends in '/' + member, from the zip directory alone."""
    rf = RangeFile(url, total)
    with zipfile.ZipFile(rf) as zf:
        infos = zf.infolist()
    hits = [i for i in infos if i.filename == member or i.filename.endswith("/" + member)]
    if len(hits) != 1:
        like = sorted(i.filename for i in infos if "transcript" in i.filename.lower())
        raise Refused(f"refused: {len(hits)} members named {member} in {url} (names with 'transcript': {like})")
    zi = hits[0]
    if zi.flag_bits & 0x1:
        raise Refused(f"refused: member {zi.filename} is encrypted")
    if zi.compress_type not in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED):
        raise Refused(f"refused: member {zi.filename} uses compression method {zi.compress_type}")
    head = _range_bytes(url, zi.header_offset, zi.header_offset + 29, total)
    if head[:4] != b"PK\x03\x04":
        raise Refused(f"refused: no local header at {zi.header_offset} for {zi.filename}")
    n_name = int.from_bytes(head[26:28], "little")
    n_extra = int.from_bytes(head[28:30], "little")
    name = _range_bytes(url, zi.header_offset + 30, zi.header_offset + 30 + n_name - 1, total)
    if name.decode("utf-8", "replace") != zi.filename and name.decode("cp437", "replace") != zi.filename:
        raise Refused(f"refused: local header name {name!r} is not {zi.filename!r}")
    data_offset = zi.header_offset + 30 + n_name + n_extra
    if data_offset + zi.compress_size > total:
        raise Refused(f"refused: member {zi.filename} runs past the end of the archive")
    return {"member": zi.filename, "method": int(zi.compress_type), "crc32": int(zi.CRC),
            "compressed_size": int(zi.compress_size), "size": int(zi.file_size), "header_offset": int(zi.header_offset),
            "data_offset": int(data_offset), "directory_requests": rf.requests, "n_members": len(infos)}


def _inflate(src: Path, dst_part: Path, size: int) -> None:
    d = zlib.decompressobj(-15)
    with open(src, "rb") as fi, open(dst_part, "wb") as fo:
        while True:
            buf = fi.read(CHUNK)
            if not buf:
                break
            fo.write(d.decompress(buf))
        fo.write(d.flush())
        if fo.tell() != size:
            raise Failed(f"{src}: inflated to {fo.tell()} bytes, the directory says {size}")


def fetch_member(entry: dict, out: Path, log=print) -> dict:
    final = out / entry["name"]
    url, total = entry["url"], entry["archive_bytes"]
    info = member_info(url, total, entry["member"])
    if final.exists():
        if final.stat().st_size != info["size"]:
            raise Refused(f"refused: {final} exists with {final.stat().st_size} bytes, the member has "
                          f"{info['size']}; nothing is deleted")
    elif info["method"] == zipfile.ZIP_STORED:
        part = out / (entry["name"] + ".part")
        _download_range(url, part, info["data_offset"], info["compressed_size"], total, log)
        os.replace(part, final)
    else:
        raw = out / (entry["name"] + ".deflate")
        raw_part = out / (entry["name"] + ".deflate.part")
        if not raw.exists():
            _download_range(url, raw_part, info["data_offset"], info["compressed_size"], total, log)
            os.replace(raw_part, raw)
        elif raw.stat().st_size != info["compressed_size"]:
            raise Refused(f"refused: {raw} has {raw.stat().st_size} bytes, the member {info['compressed_size']}")
        part = out / (entry["name"] + ".part")
        _inflate(raw, part, info["size"])
        os.replace(part, final)
    sha, crc, n = _sha256_crc(final)
    if n != info["size"] or crc != info["crc32"]:
        raise Failed(f"{final}: size {n} and CRC-32 {crc:08x}, the zip directory says {info['size']} and "
                     f"{info['crc32']:08x} (the file is kept; nothing is deleted)")
    rec = {"name": entry["name"], "url": url, "size": n, "sha256": sha, "study": entry["study"],
           "section": entry["section"], "stage": entry["stage"], "archive_bytes": total,
           "archive_md5_on_page": XENIUM_BUNDLE_MD5 if url == XENIUM_BUNDLE_URL else None}
    rec.update({k: info[k] for k in ("member", "method", "crc32", "compressed_size", "header_offset", "data_offset",
                                     "n_members")})
    return rec


# ---- the run ------------------------------------------------------------------------------------------------------
def select(pinned, only=(), urls=()) -> list[dict]:
    by_name = {e["name"]: e for e in pinned}
    by_url = {}
    for e in pinned:
        by_url.setdefault(e["url"], []).append(e)
    for u in urls:
        if u not in by_url:
            raise Refused(f"refused: {u} is not on the pinned list")
    for n in only:
        if n not in by_name:
            raise Refused(f"refused: {n} is not a pinned file name")
    if not only and not urls:
        return list(pinned)
    chosen = [by_name[n] for n in only] + [e for u in urls for e in by_url[u]]
    seen, out = set(), []
    for e in chosen:
        if e["name"] not in seen:
            seen.add(e["name"])
            out.append(e)
    return out


def _script_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def write_manifest(out: Path, records: dict) -> Path:
    path = out / "MANIFEST.json"
    doc = {"schema": 1, "script": "experiments/heart_ext/fetch_oic.py", "script_sha256": _script_sha256(),
           "constants": {"dnb_pitch_um": DNB_PITCH_UM, "mosta_bin_dnb": MOSTA_BIN_DNB, "mosta_bin_um": MOSTA_BIN_UM,
                         "xenium_grid_um": XENIUM_GRID_UM},
           "files": [records[k] for k in sorted(records)]}
    tmp = out / "MANIFEST.json.writing"
    tmp.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    os.replace(tmp, path)
    return path


def run(entries: list[dict], out: Path, jobs: int = 4, log=print) -> int:
    out.mkdir(parents=True, exist_ok=True)
    if out.is_symlink():
        raise Refused(f"refused: {out} is a symlink")
    man = out / "MANIFEST.json"
    records = {}
    if man.exists():
        for r in json.loads(man.read_text()).get("files", []):
            records[r["name"]] = r
    pinned_urls = {e["url"] for e in PINNED} | {e["url"] for e in entries}
    for r in records.values():
        if r.get("url") not in pinned_urls:
            raise Refused(f"refused: {man} records {r.get('url')}, not a pinned URL")
    lock = threading.Lock()
    errors = []

    def one(e):
        try:
            log(f"start {e['name']}")
            rec = fetch_member(e, out, log) if e["kind"] == "zip_member" else fetch_file(e, out, log)
            with lock:
                records[rec["name"]] = rec
                write_manifest(out, records)
            log(f"done  {e['name']} {rec['size']} bytes sha256 {rec['sha256']}")
        except Refused as exc:
            with lock:
                errors.append((2, str(exc)))
        except Failed as exc:
            with lock:
                errors.append((1, str(exc)))

    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        list(pool.map(one, entries))
    with lock:
        write_manifest(out, records)
    for _code, msg in errors:
        print(msg, file=sys.stderr)
    if any(c == 2 for c, _ in errors):
        return 2
    return 1 if errors else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, help="output folder (outside the repository)")
    ap.add_argument("--only", action="append", default=[], help="a pinned file name (repeatable)")
    ap.add_argument("--url", action="append", default=[], help="a pinned URL (repeatable)")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--list", action="store_true", help="print the pinned list and exit (no network)")
    a = ap.parse_args(argv)
    try:
        entries = select(PINNED, a.only, a.url)
    except Refused as exc:
        print(exc, file=sys.stderr)
        return 2
    if a.list:
        for e in entries:
            size = e["bytes"] if e["bytes"] is not None else f"member {e['member']} of {e['archive_bytes']}"
            print(f"{e['study']}\t{e['name']}\t{size}\t{e['url']}")
        print(f"MOSTA total {MOSTA_TOTAL_BYTES} bytes; the Xenium member's size is read from the bundle's directory")
        print(f"MOSTA unit bin{MOSTA_BIN_DNB} = {MOSTA_BIN_UM} um; Xenium grid {XENIUM_GRID_UM} um")
        return 0
    if a.out is None:
        print("refused: --out is required", file=sys.stderr)
        return 2
    out = a.out.expanduser().resolve()
    repo = Path(__file__).resolve().parents[2]
    in_checkout = (repo / "pyproject.toml").is_file() and (repo / "experiments" / "heart_ext").is_dir()
    if in_checkout and (out == repo or repo in out.parents):
        print(f"refused: {out} is inside the repository", file=sys.stderr)
        return 2
    try:
        return run(entries, out, a.jobs)
    except Refused as exc:
        print(exc, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
