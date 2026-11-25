"""Synthetic tests of experiments/heart_ext/fetch_oic.py: a local HTTP server with range support, no outside request."""

from __future__ import annotations

import hashlib
import http.server
import importlib.util
import io
import json
import re
import subprocess
import sys
import threading
import zipfile
import zlib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "heart_ext" / "fetch_oic.py"


def _load():
    spec = importlib.util.spec_from_file_location("fetch_oic_t", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.RETRIES = 3
    return mod


class _Handler(http.server.BaseHTTPRequestHandler):
    files: dict = {}
    log: list = []
    cut_after: dict = {}

    def log_message(self, *a):  # quiet
        pass

    def _send(self, body_ok: bool):
        data = self.files.get(self.path)
        if data is None:
            self.send_response(404)
            self.end_headers()
            return
        rng = self.headers.get("Range")
        self.log.append((self.command, self.path, rng))
        if rng:
            m = re.fullmatch(r"bytes=(\d+)-(\d*)", rng)
            a = int(m.group(1))
            b = int(m.group(2)) if m.group(2) else len(data) - 1
            b = min(b, len(data) - 1)
            chunk = data[a: b + 1]
            cut = self.cut_after.pop(self.path, None)
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {a}-{b}/{len(data)}")
            self.send_header("Content-Length", str(len(chunk)))
            self.end_headers()
            if body_ok:
                if cut is not None:
                    self.wfile.write(chunk[:cut])
                    self.wfile.flush()
                    self.connection.shutdown(2)
                    return
                self.wfile.write(chunk)
            return
        self.send_response(200)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        if body_ok:
            self.wfile.write(data)

    def do_GET(self):
        self._send(True)

    def do_HEAD(self):
        self._send(False)


@pytest.fixture()
def server():
    _Handler.files = {}
    _Handler.log = []
    _Handler.cut_after = {}
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    yield httpd, f"http://127.0.0.1:{httpd.server_address[1]}"
    httpd.shutdown()


def _bundle(method, force64=False):
    """A bundle like the vendor's: refused members around the one transcript table."""
    table = bytes(range(256)) * 997 + b"tail"
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr(zipfile.ZipInfo("outs/cells.parquet"), b"C" * 5000, compress_type=zipfile.ZIP_DEFLATED)
        zf.writestr(zipfile.ZipInfo("outs/transcripts.zarr.zip"), b"Z" * 3000, compress_type=zipfile.ZIP_DEFLATED)
        zi = zipfile.ZipInfo("outs/transcripts.parquet")
        zi.compress_type = method
        with zf.open(zi, "w", force_zip64=force64) as fh:
            fh.write(table)
        zf.writestr(zipfile.ZipInfo("outs/cell_boundaries.parquet"), b"B" * 7000, compress_type=zipfile.ZIP_DEFLATED)
    return buf.getvalue(), table


def _entries(mod, base, files):
    out = []
    for name, data in files.items():
        out.append({"name": name, "url": base + "/" + name, "bytes": len(data), "kind": "file", "study": "X022",
                    "section": name.split("_GEM")[0], "stage": name.split("_")[0]})
    return out


def test_constants_named_before_fetch():
    mod = _load()
    assert mod.MOSTA_BIN_DNB == 24 and mod.DNB_PITCH_UM == 0.5 and mod.XENIUM_GRID_UM == 12.0
    doc = mod.__doc__
    assert "bin24" in doc and "12 um" in doc and doc.index("bin24") < doc.index("Behaviour")
    names = [e["name"] for e in mod.PINNED]
    assert len([n for n in names if n.endswith("_GEM_bin1.tsv.gz")]) == 23
    assert sum(1 for n in names if n.startswith("E15.5_")) == 5
    assert sum(1 for n in names if n.startswith("E16.5_")) == 18
    assert all(e["url"].startswith("https://") for e in mod.PINNED)
    assert mod.MOSTA_TOTAL_BYTES == 48844857975
    assert not any("h5ad" in e["url"] or "cell_bin" in e["url"].lower() for e in mod.PINNED)
    assert [e["member"] for e in mod.PINNED if e["kind"] == "zip_member"] == ["transcripts.parquet"]


def test_source_deletes_nothing():
    src = SCRIPT.read_text()
    for word in ("unlink(", "rmtree(", "os.remove(", "rmdir(", ".truncate(", "shutil"):
        assert word not in src


def test_whole_file_resume_and_manifest(server, tmp_path):
    mod = _load()
    _, base = server
    files = {"E15.5_E1S1_GEM_bin1.tsv.gz": bytes(range(256)) * 1000, "E16.5_E2S9_GEM_bin1.tsv.gz": b"x" * 77777}
    _Handler.files = {"/" + k: v for k, v in files.items()}
    ents = _entries(mod, base, files)
    out = tmp_path / "out"
    out.mkdir()
    (out / "E15.5_E1S1_GEM_bin1.tsv.gz.part").write_bytes(files["E15.5_E1S1_GEM_bin1.tsv.gz"][:1000])
    _Handler.cut_after["/E16.5_E2S9_GEM_bin1.tsv.gz"] = 5000      # a dropped connection, then resumed
    mod.PINNED = tuple(ents)
    assert mod.run(ents, out, jobs=2, log=lambda *_: None) == 0
    for k, v in files.items():
        assert (out / k).read_bytes() == v
        assert not (out / (k + ".part")).exists()
    ranges = [r for c, p, r in _Handler.log if p == "/E15.5_E1S1_GEM_bin1.tsv.gz"]
    assert ranges[0].startswith("bytes=1000-")
    man = json.loads((out / "MANIFEST.json").read_text())
    recs = {r["name"]: r for r in man["files"]}
    for k, v in files.items():
        assert recs[k]["sha256"] == hashlib.sha256(v).hexdigest() and recs[k]["size"] == len(v)
        assert recs[k]["url"] == base + "/" + k
    assert man["constants"]["mosta_bin_dnb"] == 24 and man["constants"]["xenium_grid_um"] == 12.0
    # a second call fetches nothing again and keeps the records
    n_before = len(_Handler.log)
    assert mod.run(ents, out, jobs=1, log=lambda *_: None) == 0
    assert len(_Handler.log) == n_before


def test_size_mismatch_refused_and_nothing_deleted(server, tmp_path):
    mod = _load()
    _, base = server
    files = {"E15.5_E1S2_GEM_bin1.tsv.gz": b"y" * 4000}
    _Handler.files = {"/" + k: v for k, v in files.items()}
    ents = _entries(mod, base, files)
    ents[0]["bytes"] = 3999                       # the pin differs from the server
    out = tmp_path / "o"
    assert mod.run(ents, out, log=lambda *_: None) == 2
    assert not (out / "E15.5_E1S2_GEM_bin1.tsv.gz").exists()
    # a .part longer than the pin is refused and kept
    ents[0]["bytes"] = 4000
    (out / "E15.5_E1S2_GEM_bin1.tsv.gz.part").write_bytes(b"y" * 4500)
    assert mod.run(ents, out, log=lambda *_: None) == 2
    assert (out / "E15.5_E1S2_GEM_bin1.tsv.gz.part").stat().st_size == 4500


@pytest.mark.parametrize("method,force64", [(zipfile.ZIP_STORED, False), (zipfile.ZIP_DEFLATED, False),
                                            (zipfile.ZIP_DEFLATED, True)])
def test_zip_member_only_its_bytes(server, tmp_path, method, force64):
    mod = _load()
    _, base = server
    bundle, table = _bundle(method, force64)
    _Handler.files = {"/outs.zip": bundle}
    ent = {"name": "X021_t.parquet", "url": base + "/outs.zip", "bytes": None, "kind": "zip_member",
           "archive_bytes": len(bundle), "member": "transcripts.parquet", "study": "X021", "section": "pup_P1",
           "stage": "P1"}
    mod.PINNED = (ent,)
    out = tmp_path / "x"
    assert mod.run([ent], out, log=lambda *_: None) == 0
    assert (out / "X021_t.parquet").read_bytes() == table
    rec = json.loads((out / "MANIFEST.json").read_text())["files"][0]
    assert rec["member"] == "outs/transcripts.parquet" and rec["crc32"] == zlib.crc32(table)
    assert rec["sha256"] == hashlib.sha256(table).hexdigest()
    # no byte of another member was requested: every range lies in the directory or in this member's span
    with zipfile.ZipFile(io.BytesIO(bundle)) as zf:
        spans = []
        for zi in zf.infolist():
            if zi.filename != "outs/transcripts.parquet":
                spans.append((zi.header_offset + 30 + len(zi.filename), zi.header_offset + 30 + len(zi.filename)
                              + zi.compress_size))
    for _c, _p, r in _Handler.log:
        a, b = r[6:].split("-")
        a, b = int(a), int(b) if b else len(bundle) - 1
        for s, e in spans:
            assert b < s or a >= e, (r, s, e)
    if method == zipfile.ZIP_DEFLATED:
        assert (out / "X021_t.parquet.deflate").exists()      # kept, nothing deleted


def test_zip_member_missing_refused(server, tmp_path):
    mod = _load()
    _, base = server
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("outs/cells.parquet", b"c" * 100)
    _Handler.files = {"/b.zip": buf.getvalue()}
    ent = {"name": "t.parquet", "url": base + "/b.zip", "bytes": None, "kind": "zip_member",
           "archive_bytes": len(buf.getvalue()), "member": "transcripts.parquet", "study": "X021",
           "section": "pup_P1", "stage": "P1"}
    mod.PINNED = (ent,)
    assert mod.run([ent], tmp_path / "y", log=lambda *_: None) == 2


def test_cli_refuses_unpinned_url_and_names(tmp_path):
    for extra in (["--url", "https://example.org/E15.5_E1S1_GEM_bin1.tsv.gz"],
                  ["--url", "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/E15.5_E1S1.MOSTA.h5ad"],
                  ["--only", "E14.5_E1S1_GEM_bin1.tsv.gz"]):
        p = subprocess.run([sys.executable, "-I", "-B", str(SCRIPT), "--out", str(tmp_path / "z"), *extra],
                           capture_output=True, text=True, timeout=60)
        assert p.returncode == 2 and "refused" in p.stderr
    assert not (tmp_path / "z").exists()
    p = subprocess.run([sys.executable, "-I", "-B", str(SCRIPT), "--out", str(ROOT / "data" / "oic")],
                       capture_output=True, text=True, timeout=60)
    assert p.returncode == 2 and "inside the repository" in p.stderr


def test_cli_list_no_network():
    p = subprocess.run([sys.executable, "-I", "-B", str(SCRIPT), "--list"], capture_output=True, text=True, timeout=60)
    assert p.returncode == 0
    lines = p.stdout.splitlines()
    assert sum(1 for x in lines if x.startswith("X022\t") and x.endswith("_GEM_bin1.tsv.gz")) == 23
    assert sum(1 for x in lines if x.startswith("X021\t") and "transcripts.parquet" in x) == 1
    assert "bin24 = 12.0 um" in p.stdout


def test_check_folder_against_the_manifest_and_the_chords(tmp_path):
    import os

    mod = _load()
    e = next(x for x in mod.PINNED if x["kind"] == "file")
    f = tmp_path / e["name"]
    with open(f, "wb"):
        pass
    os.truncate(f, e["bytes"])  # sparse: the pinned size, no disk used
    rec = {"name": e["name"], "url": e["url"], "size": e["bytes"], "sha256": "a" * 64}
    (tmp_path / "MANIFEST.json").write_text(json.dumps({"files": [rec]}))
    assert mod.check_folder(tmp_path) == []
    ch = tmp_path / "chords.json"
    ch.write_text(json.dumps({"sources": {e["name"]: "a" * 64}}))
    why = mod.check_folder(tmp_path, [ch])  # one file of X022's 23: the source is incomplete, so dropped
    assert len(why) == 1 and "whose source is dropped" in why[0]
    ch.write_text(json.dumps({"sources": {e["name"]: "b" * 64}}))
    assert mod.check_folder(tmp_path, [ch])
    ch.write_text(json.dumps({"sources": {}}))
    assert mod.check_folder(tmp_path, [ch])
    for bad in ({**rec, "url": "https://example.org/x"}, {**rec, "size": 5}, {**rec, "sha256": None}):
        (tmp_path / "MANIFEST.json").write_text(json.dumps({"files": [bad]}))
        assert mod.check_folder(tmp_path), bad
    (tmp_path / "MANIFEST.json").write_text(json.dumps({"files": [rec]}))
    os.truncate(f, 10)
    assert mod.check_folder(tmp_path)
    p = subprocess.run([sys.executable, "-I", "-B", str(SCRIPT), "--check", str(tmp_path)], capture_output=True,
                       text=True, timeout=60)
    assert p.returncode == 2 and "refused" in p.stderr


def test_records_carry_their_completion_time(server, tmp_path):
    import os
    import time

    mod = _load()
    _, base = server
    files = {"E15.5_E1S3_GEM_bin1.tsv.gz": b"z" * 999}
    _Handler.files = {"/" + k: v for k, v in files.items()}
    ents = _entries(mod, base, files)
    mod.PINNED = tuple(ents)
    out = tmp_path / "t"
    assert mod.run(ents, out, log=lambda *_: None) == 0
    rec = json.loads((out / "MANIFEST.json").read_text())["files"][0]
    mtime = os.stat(out / "E15.5_E1S3_GEM_bin1.tsv.gz").st_mtime
    assert rec["completed_utc"] == time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(mtime))


def test_source_status_drops_a_late_or_incomplete_source(tmp_path):
    import os

    mod = _load()
    before, cut = mod.SOURCE_CUTOFF_EPOCH - 3600, mod.SOURCE_CUTOFF_EPOCH
    files = []
    for e in (x for x in mod.PINNED if x["study"] == "X022"):
        f = tmp_path / e["name"]
        with open(f, "wb"):
            pass
        os.truncate(f, e["bytes"])
        os.utime(f, (before, before))
        files.append({"name": e["name"], "url": e["url"], "size": e["bytes"], "sha256": "a" * 64,
                      "completed_utc": mod.completed_utc(f)})
    (tmp_path / "MANIFEST.json").write_text(json.dumps({"files": files}))
    st = mod.source_status(tmp_path)
    assert st["X022"]["kept"] and not st["X021"]["kept"]  # X021 not recorded: dropped
    last = tmp_path / files[-1]["name"]
    os.utime(last, (cut, cut))  # exactly at the cut-off: still complete by it, but it changed after its record
    st = mod.source_status(tmp_path)
    assert not st["X022"]["kept"] and "changed after it was recorded" in st["X022"]["why"][0]
    files[-1]["completed_utc"] = mod.completed_utc(last)
    (tmp_path / "MANIFEST.json").write_text(json.dumps({"files": files}))
    assert mod.source_status(tmp_path)["X022"]["kept"]
    os.utime(last, (cut + 1, cut + 1))
    files[-1]["completed_utc"] = mod.completed_utc(last)
    (tmp_path / "MANIFEST.json").write_text(json.dumps({"files": files}))
    st = mod.source_status(tmp_path)
    assert not st["X022"]["kept"] and "completed after" in st["X022"]["why"][0]
    p = subprocess.run([sys.executable, "-I", "-B", str(SCRIPT), "--sources", str(tmp_path)], capture_output=True,
                       text=True, timeout=60)
    assert p.returncode == 0 and json.loads(p.stdout)["X022"]["kept"] is False
