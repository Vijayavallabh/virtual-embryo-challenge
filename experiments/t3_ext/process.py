#!/usr/bin/env python
"""One fetched outside-knockout dataset to its factor tables and one QC file (experiments/t3_ext/README.md).

  reference --wt WT_E8.75.h5ad --out D/reference/refstats.npz
      the WT E8.75 reference's panel genes (its gene order), lineage groups (obs['group'] of vec.prep.load_stage,
      'Unknown' left out), their shares and their mean within-cell proportions. Read by run and by combine.
  run ACC --dest D --wt WT_E8.75.h5ad --refstats D/reference/refstats.npz --out D/processed
      reads D/meta/ACC/plan.json and D/raw/ACC/MANIFEST.json (every file's sha256 checked). It touches ONLY the
      samples in datasets.yaml's pinned block: it refuses a dataset with no pins, a manifest fetched under other pins,
      a pinned series whose stored text no longer hashes to its pin, a pinned accession absent from its own series'
      plan table, a pinned sample whose record's series_id lists a series outside screen_series (or a held one, or
      none), a pinned sample whose stored metadata no longer hashes to its pin, a file whose manifest record was
      written under another metadata_sha256, url or pins, a file in a pinned sample's folder that is not one of its
      pinned files with a record, and screens each pinned sample's record again with t3x_screen BEFORE any count is
      read (a sample that does not pass refuses the run). Arms, genes and stages come from the pins only; no file of an
      unpinned sample is opened. It reads each pinned sample's counts restricted to the panel genes, labels every
      cell with a lineage group by transfer
      from the WT reference, and computes one log-ratio table per knockout gene with the fixed guards, minimum-cells
      rule and shrinkage below. It writes per gene <out>/<ACC>/<ACC>__<gene>.unit.json (log space, for combine) and
      its finished tables (pooled and per group), and <out>/<ACC>/QC.json. It prints no gene-level value.
      Stage 5: it first deletes <out>/<ACC> (so a run that refuses leaves no output combine could read); it refuses a
      symlinked raw, dataset or sample folder, manifest or file; it takes each pinned file's url only from the list
      verify_pins derives from the hashed fields; QC.json records plan_sha256 (the sha256 of the plan it ran under,
      checked by combine), the manifest's sha256 and each count file's sha256.
      Stage 6: a sample pinned with files [] takes its cells from its pinned series file (read_series_unit: an archive
      read member by member in memory after archive_members has checked every member, or a loose triple), split by
      the pin's rule; a unit with a control from the named control series takes d = t3x_common.SISTER_D, and QC.json
      says so (series_file_units, sister_series_controls, per unit d and control_source).
      Stage 7: a pinned sample whose file a raw_matrices declaration names is an unfiltered matrix: read sparsely in
      two streamed passes (read_raw_sample), the ONE fixed cell call of t3x_common.raw_cell_call keeps its cells, and
      the ordinary per-cell rule and barcode cap then apply to them; QC.json holds raw_cell_call (the rule and each
      raw sample's record) and each unit's QC raw_called and raw_called_samples. Nothing else downstream differs.
      Stage 8: a unit whose samples are not all raw-called or all filtered is refused (a second guard behind
      fetch.check_unit_forms); a raw-called unit's QC holds raw_call_arms, the per-arm medians of the cutoff, the cells
      called and the UMI share in called cells and their knockout-to-control ratios (a record only); a matrix line
      that is not three fields, and more than 1,000 comment lines before the size line, are refused.

The rules (constants in t3x_common.py, fixed 2026-11-20, none from a score):
  * Cells: stage capped at E9.5 (the T1 window applied to every cell); fewer than 200 detected genes is not a cell;
    a per-sample matrix of more than 60,000 barcodes is refused as unfiltered.
  * Groups: a logistic-regression transfer (vec.harmonise.fit_transfer, the frozen settings: 50 PCs, C 1.0, clip 10,
    5 depth draws, seed 0) trained on the reference's lineage GROUPS directly, so the outside call and the reference
    cell's own obs['group'] mean the same thing (the challenger's change 2). Confident: label_conf >= 0.6.
  * A stage drops if fewer than half of either arm's cells are confident. A (group, stage) drops if: fewer than 100
    confident cells in any arm; the group's share of the arm's cells differs across arms by more than 2x (a knockout
    that moves cells between calls is not read as expression change); the control arm's group profile does not have
    the reference group as its best Spearman match, ahead of the next by 0.05; or, where author labels and a committed
    dictionary exist, fewer than half the control cells' dictionary groups agree.
  * Each sample's confident and unconfident shares go into QC.json; a knockout gene's table is refused (status
    refused_low_confidence, nothing written) if any of its samples has a confident share under datasets.yaml's
    qc_floors.min_confident_share_per_sample.
  * Log ratio per gene on summed UMIs with a pseudocount of 10: log((x_ko+10)/N_ko) - log((x_ctl+10)/N_ctl); NaN (not
    measured) for a gene with fewer than 20 UMIs in any arm; clipped to log [0.25, 4]. Chimeras only with a
    wild-type-chimera second difference: [mut vs host] in knockout chimeras minus [mut vs host] in wild-type
    chimeras.
  * Shrinkage L = k * log ratio, k = min(n_h/(n_h+200), m/(m+1)) * m/(m+1) * d (t3x_common.shrink_k).
  * Per group, its passing stages are averaged with equal weight. No group passing = the dataset is left out
    (status 'left_out'), never given a pooled fallback.
"""

from __future__ import annotations

import gzip
import io
import json
import re
import shutil
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import fetch as fe  # noqa: E402  (experiments/t3_ext/fetch.py: the plan loader and the pin check)
import t3x_common as cm  # noqa: E402
import t3x_screen as sc  # noqa: E402
import t3x_tables as tb  # noqa: E402

Refused = cm.Refused
CONTROL_ARMS = ("control", "ko_host", "wt_host")
PLAIN_ARMS = ("ko", "control")
CHIMERA_ARMS = ("ko_mut", "ko_host", "wt_mut", "wt_host")


# ---------------------------------------------------------------------------------------------------------------------
# Reading one sample's counts (no extraction into any folder)


MEMBER_MAX_BYTES = 8 * 1024 ** 3   # one archive member read into memory: at most 8 GiB (fixed 2026-11-20, no score)
MAC_APPLEDOUBLE_MAX = 4096   # a macOS tar companion member ('._name'): never read; a larger one is refused (2026-11-20)
MAC_DS_STORE_MAX = 65536     # a macOS Finder '.DS_Store' member: never read; a larger one is refused (2026-11-20)
MAC_SKIPPED_TOTAL_MAX = 1024 ** 2   # all skipped macOS members of one archive together: at most 1 MiB, else refused


def _gunzip(raw: bytes, name: str) -> bytes:
    """A gzip stream decompressed in chunks, refused once it passes MEMBER_MAX_BYTES (stage 7: the first reviewer of
    853b214, lesser point 2; the 8 GiB check on an archive member covers its compressed size only)."""
    out, n = io.BytesIO(), 0
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as g:
        while block := g.read(1 << 24):
            n += len(block)
            if n > MEMBER_MAX_BYTES:
                raise Refused(f"{name}: decompresses to more than {MEMBER_MAX_BYTES} bytes: refused")
            out.write(block)
    return out.getvalue()


def _read_text(raw: bytes, name: str) -> list[str]:
    if name.endswith(".gz"):
        raw = _gunzip(raw, name)
    return raw.decode("utf-8", errors="replace").splitlines()


_NOT_COMMENT = rb"(?m)^(?![ \t]*%)[ \t]*\S+"
# a line that is not a comment with four or more fields, or with one or two
_FOUR_FIELDS = re.compile(_NOT_COMMENT + rb"(?:[ \t]+\S+){3,}")
_FEW_FIELDS = re.compile(_NOT_COMMENT + rb"(?:[ \t]+\S+)?[ \t]*\r?$")


def _three_fields(data: bytes, mtx_name: str) -> None:
    """Stage 9 round 2 (the verification of 9edee21, lesser point 1): on the filtered route too, every line after the
    size line that is not a comment or empty is 'feature barcode value' (three fields), checked before scipy's mmread,
    which takes the first three fields of a longer line and drops the rest."""
    fh = io.BytesIO(data)
    _mtx_header(fh, mtx_name)
    body = data[fh.tell():]
    if _FOUR_FIELDS.search(body) or _FEW_FIELDS.search(body):
        raise Refused(f"{mtx_name}: a matrix line is not 'feature barcode value': refused")


def _from_parts(mtx: bytes, mtx_name: str, feats: list[str], bcs: list[str]):
    data = _gunzip(mtx, mtx_name) if mtx_name.endswith(".gz") else mtx
    _three_fields(data, mtx_name)
    M = sparse.csr_matrix(spio.mmread(io.BytesIO(data), spmatrix=False).T)  # cells x features
    rows = [r.split("\t") for r in feats if r.strip()]
    symbols = [r[1] if len(r) > 1 else r[0] for r in rows]
    types = [r[2] if len(r) > 2 else "Gene Expression" for r in rows]
    barcodes = [b.split("\t")[0] for b in bcs if b.strip()]
    if M.shape != (len(barcodes), len(symbols)):
        raise Refused(f"matrix shape {M.shape} does not match {len(barcodes)} barcodes x {len(symbols)} features")
    return M, symbols, types, barcodes


def _pick(names: list[str], key: str) -> str:
    hits = [n for n in names if Path(n).name.lower().endswith((key, key + ".gz"))]
    if len(hits) != 1:
        raise Refused(f"expected exactly one '{key}' file, found {len(hits)}")
    return hits[0]


def read_sample(paths: list[Path]):
    """(cells x features CSR, symbols, feature types, barcodes) from a sample's mtx triple, its 10x h5, or a tar.gz
    holding a mtx triple (read member by member in memory)."""
    paths = [Path(p) for p in paths]
    h5 = [p for p in paths if p.name.endswith(".h5")]
    tars = [p for p in paths if p.name.endswith((".tar.gz", ".tgz", ".tar"))]
    if h5:
        if len(h5) != 1:
            raise Refused(f"{len(h5)} h5 files for one sample")
        import h5py

        with h5py.File(h5[0], "r") as f:
            g = f["matrix"]
            shape = tuple(g["shape"][()])
            M = sparse.csc_matrix((g["data"][()], g["indices"][()], g["indptr"][()]), shape=shape).T.tocsr()

            def dec(a):
                return [x.decode() if isinstance(x, bytes) else str(x) for x in a]

            symbols = dec(g["features"]["name"][()])
            types = dec(g["features"]["feature_type"][()])
            barcodes = dec(g["barcodes"][()])
        return M, symbols, types, barcodes
    if tars:
        if len(tars) != 1 or len(paths) != 1:
            raise Refused(f"{len(tars)} archives (and {len(paths) - len(tars)} other files) for one sample: one "
                          "archive alone")
        # stage 7 (both reviewers of 853b214): the per-sample archive route checks every member with the same rules
        # as a series archive, and holds exactly one 10x triple in one folder, nothing else
        return _read_sample_archive(tars[0])
    names = [str(p) for p in paths]
    mtx = _pick(names, "matrix.mtx")
    feat = [n for n in names if "features.tsv" in Path(n).name.lower() or "genes.tsv" in Path(n).name.lower()]
    if len(feat) != 1:
        raise Refused(f"expected one features/genes file, found {len(feat)}")
    bc = _pick(names, "barcodes.tsv")
    return _from_parts(Path(mtx).read_bytes(), mtx, _read_text(Path(feat[0]).read_bytes(), feat[0]),
                       _read_text(Path(bc).read_bytes(), bc))


# ---------------------------------------------------------------------------------------------------------------------
# Stage 6: a pinned series file (one archive or a loose 10x triple), read only inside itself and split by its pin's rule

_MEMBER_PART = re.compile(r"[A-Za-z0-9._+-]+")


def archive_members(t: tarfile.TarFile, label: str, raw_ok: bool = False,
                    skipped: list[str] | None = None) -> list[tuple[list[str], str, tarfile.TarInfo]]:
    """The archive's file members as (path parts, triple kind, member), every member checked first: refused for an
    absolute path, a '..' component, an empty component, a backslash, a name part outside letters, digits and
    . _ + -, a symlink, a hard link, a device or FIFO, a duplicate name, a name the held-out screen hits (held out or
    review), a member larger than MEMBER_MAX_BYTES, and a file member that is not part of a 10x triple (a member it has
    no rule for). Directory members are allowed and read nothing. Nothing is extracted to disk. raw_ok (stage 7: a
    declared raw file only, read_raw_sample) lifts the unfiltered-matrix deny term alone for the member names.

    macOS companion members (2026-11-20, after a real archive held them): a file member named '._X' (the AppleDouble
    companion macOS tar writes beside a member X) and a file member named '.DS_Store' are NEVER READ and never classed
    as part of a triple; they pass every check above first, and are skipped only when a '._X' member is at most
    MAC_APPLEDOUBLE_MAX bytes and the member X it names (a file or a folder) is in the archive beside it, and a
    '.DS_Store' member is at most MAC_DS_STORE_MAX bytes; together they hold at most MAC_SKIPPED_TOTAL_MAX bytes;
    otherwise the archive is refused. Their names are appended to `skipped` when a list is given."""
    out, seen, companions, mac, mac_bytes = [], set(), [], [], 0
    for m in t.getmembers():
        raw = m.name
        if raw.startswith(("/", "~")) or "\\" in raw or re.match(r"[A-Za-z]:", raw):
            raise Refused(f"{label}: member {raw!r}: an absolute path (or a backslash or drive letter): refused")
        parts = raw.split("/")
        if ".." in parts:
            raise Refused(f"{label}: member {raw!r}: a '..' component (path traversal): refused")
        if m.isdir() and len(parts) > 1 and parts[-1] == "":
            parts = parts[:-1]
        if parts[:1] == ["."]:
            parts = parts[1:]
        if any(p in ("", ".") for p in parts):
            raise Refused(f"{label}: member {raw!r}: an empty or '.' path component: refused")
        if m.issym() or m.islnk():
            raise Refused(f"{label}: member {raw!r}: a symlink or hard link: refused")
        if not (m.isfile() or m.isdir()):
            raise Refused(f"{label}: member {raw!r}: not a regular file or folder: refused")
        if not all(_MEMBER_PART.fullmatch(p) for p in parts):
            raise Refused(f"{label}: a member name holds a character outside letters, digits and . _ + -: refused")
        key = "/".join(parts)
        if key in seen:
            raise Refused(f"{label}: member {key!r} appears twice: refused")
        seen.add(key)
        if parts and sc.hits(key):
            raise Refused(f"{label}: a member name hits the held-out screen: refused (name withheld)")
        if parts and fe.denied(key, raw_ok=raw_ok):
            raise Refused(f"{label}: member {key!r} is a denied type (unfiltered, fragments, peaks ...): refused")
        if m.isdir():
            continue
        if not parts:
            raise Refused(f"{label}: a file member without a name: refused")
        last = parts[-1]
        if last.startswith("._") and len(last) > 2:
            if not 0 <= m.size <= MAC_APPLEDOUBLE_MAX:
                raise Refused(f"{label}: member {key!r}: a '._' member larger than {MAC_APPLEDOUBLE_MAX} bytes, or of "
                              "a negative size: refused")
            companions.append((key, "/".join([*parts[:-1], last[2:]])))
            mac.append(key)
            mac_bytes += m.size
            continue
        if last == ".DS_Store":
            if not 0 <= m.size <= MAC_DS_STORE_MAX:
                raise Refused(f"{label}: member {key!r}: a '.DS_Store' member larger than {MAC_DS_STORE_MAX} bytes "
                              "(or of a negative size): refused")
            mac.append(key)
            mac_bytes += m.size
            continue
        kind = cm.triple_kind(last)
        if kind is None:
            raise Refused(f"{label}: member {key!r} is not part of a 10x triple: refused (a member with no rule)")
        if not 0 <= m.size <= MEMBER_MAX_BYTES:
            raise Refused(f"{label}: member {key!r} is larger than {MEMBER_MAX_BYTES} bytes (or of a negative size): "
                          "refused")
        out.append((parts, kind, m))
    if mac_bytes > MAC_SKIPPED_TOTAL_MAX:
        raise Refused(f"{label}: the skipped macOS members hold more than {MAC_SKIPPED_TOTAL_MAX} bytes together: "
                      "refused")
    for key, beside in companions:
        if beside not in seen:
            raise Refused(f"{label}: member {key!r} names no member beside it ({beside!r} is not in the archive): "
                          "refused")
    if skipped is not None:
        skipped.extend(sorted(mac))
    return out


def _one_triple(items: list, label: str) -> dict:
    by = {}
    for it in items:
        if it[1] in by:
            raise Refused(f"{label}: two '{it[1]}' members in one library: refused")
        by[it[1]] = it
    if sorted(by) != ["barcodes", "features", "matrix"]:
        raise Refused(f"{label}: a library needs exactly one matrix, barcodes and features (or genes) member, found "
                      f"{sorted(by)}")
    return by


def _triple_from(by: dict, read) -> tuple:
    def name(k):
        return "/".join(by[k][0]) if isinstance(by[k][0], list) else by[k][0]

    return _from_parts(read(by["matrix"]), name("matrix"), _read_text(read(by["features"]), name("features")),
                       _read_text(read(by["barcodes"]), name("barcodes")))


def _one_folder_triple(items: list, label: str) -> dict:
    """The single 10x triple of a per-sample archive: its three members in one folder (or none), nothing else
    (archive_members has already refused any member that is not part of a triple)."""
    if len({tuple(it[0][:-1]) for it in items}) > 1:
        raise Refused(f"{label}: members in more than one folder: a per-sample archive holds one library, refused")
    return _one_triple(items, label)


def _read_sample_archive(path: Path):
    """A per-sample archive read member by member in memory, every member checked by archive_members first."""
    with tarfile.open(path, "r:*") as t:
        items = archive_members(t, path.name)
        by = _one_folder_triple(items, path.name)

        def read(it):
            f = t.extractfile(it[2])
            if f is None:
                raise Refused(f"{path.name}: member {'/'.join(it[0])!r} cannot be read")
            return f.read()

        return _triple_from(by, read)


def read_series_unit(paths: list[Path], unit: dict) -> list[tuple[str, str, list[str], tuple]]:
    """The libraries of one pinned series file, split by its pin's rule: [(GSM, folder or '', barcodes kept, (M,
    symbols, types, barcodes))]. An archive is read member by member in memory (never extracted; archive_members
    checks every member first); a loose triple is read from its three files. split none: the file is one unit (its
    samples are one arm and one stage, each pinned for the same genes, checked when the pins load); each folder of an
    archive (a member's parent path) must hold exactly one triple, the file must hold exactly as many libraries as
    the entry names samples (stage 7), and every library is returned under the entry's first sample.
    barcode_suffix: every barcode must end with exactly one key of the map (with something before it); a barcode no
    key fits refuses the file, and so does a key whose sample gets no barcode. folder: every file member lies in a
    top-level folder that the map names; each folder holds exactly one triple. A matrix of more barcodes than
    MAX_BARCODES per library is refused."""
    split = unit["split"]
    label = f"{unit['series']}/{'+'.join(unit['files'])}"
    libs = []
    if len(paths) == 1 and cm.series_file_form(unit["files"]) == "archive":
        with tarfile.open(paths[0], "r:*") as t:
            items = archive_members(t, label)

            def read(it):
                f = t.extractfile(it[2])
                if f is None:
                    raise Refused(f"{label}: member {'/'.join(it[0])!r} cannot be read")
                return f.read()

            if split["by"] == "folder":
                groups: dict[str, list] = {}
                for it in items:
                    if len(it[0]) < 2 or it[0][0] not in split["map"]:
                        raise Refused(f"{label}: member {'/'.join(it[0])!r} lies in no folder the split map names: "
                                      "refused")
                    groups.setdefault(it[0][0], []).append(it)
                if set(groups) != set(split["map"]):
                    raise Refused(f"{label}: folders {sorted(set(split['map']) - set(groups))} of the split map hold "
                                  "no library: refused")
                for folder in sorted(groups):
                    libs.append((split["map"][folder], folder, _triple_from(_one_triple(groups[folder], label), read)))
            elif split["by"] == "none":
                # used as one unit: each folder (a member's parent path) holds one library, all pooled under one label
                groups = {}
                for it in items:
                    groups.setdefault("/".join(it[0][:-1]), []).append(it)
                _none_count(label, len(groups), unit)
                for folder in sorted(groups):
                    libs.append((None, folder, _triple_from(_one_triple(groups[folder], label), read)))
            else:
                libs.append((None, "", _triple_from(_one_triple(items, label), read)))
    else:
        if split["by"] == "none":
            _none_count(label, 1, unit)   # a loose triple is one library
        items = [(p.name, cm.triple_kind(p.name), p) for p in paths]
        libs.append((None, "", _triple_from(_one_triple(items, label), lambda it: Path(it[2]).read_bytes())))
    out = []
    for gsm, folder, (M, symbols, types, bcs) in libs:
        if split["by"] == "barcode_suffix":
            keys = list(split["map"])
            key_of = []
            for b in bcs:
                fit = [k for k in keys if b.endswith(k) and len(b) > len(k)]
                if len(fit) != 1:
                    raise Refused(f"{label}: a barcode that no split key (or more than one) fits: refused")
                key_of.append(fit[0])
            key_of = np.array(key_of, dtype=object)
            for k in keys:
                n_k = int((key_of == k).sum())
                if n_k > cm.MAX_BARCODES:
                    raise Refused(f"{label}: {n_k} barcodes under {k}: an unfiltered matrix, refused")
            owner = np.array([split["map"][k] for k in key_of], dtype=object)
            for g in unit["samples"]:
                rows = np.flatnonzero(owner == g)
                if not len(rows):
                    raise Refused(f"{label}: no barcode carries a split key of {g}: the file does not fit its rule, "
                                  "refused")
                kept = [bcs[i] for i in rows]
                out.append((g, folder, kept, (M[rows], symbols, types, kept)))
        else:
            # stage 7: the cap applies to each library (each folder, or the one loose triple), under every split
            if M.shape[0] > cm.MAX_BARCODES:
                raise Refused(f"{label}{'/' + folder if folder else ''}: {M.shape[0]} barcodes in one library: an "
                              "unfiltered matrix, refused")
            out.append((gsm or unit["samples"][0], folder, bcs, (M, symbols, types, bcs)))
    return out


def _none_count(label: str, n_libs: int, unit: dict) -> None:
    """Stage 7 (the first reviewer of 853b214, blocker 1): under split none a file holds exactly as many libraries
    (folders with a triple, or the one loose or flat triple) as its entry names samples; otherwise refused before any
    count is read."""
    m = len(unit["samples"])
    if n_libs != m:
        raise Refused(f"{label}: split none: the file holds {n_libs} libraries but the entry names {m} samples: "
                      "refused (state a folder or suffix rule from the records, or do not pin)")


# ---------------------------------------------------------------------------------------------------------------------
# Stage 7, piece 2: a per-sample UNFILTERED matrix (a raw_matrices dataset's declared file), read sparsely in two
# streamed passes and through the ONE fixed cell call (t3x_common.raw_cell_call); never densified, never unpacked

RAW_CHUNK_BYTES = 64 * 1024 ** 2   # decompressed matrix bytes parsed per block
RAW_MAX_LINE = 1024 ** 2           # a matrix line longer than this refuses the file
RAW_MAX_COMMENT_LINES = 1000       # stage 8: more comment or empty lines before the size line refuse the file


def _decompressed(fh, name: str):
    """A binary stream of one member or file, decompressed when its name ends in .gz."""
    return gzip.GzipFile(fileobj=fh) if name.endswith(".gz") else fh


def _read_capped(fh, name: str) -> bytes:
    """A whole (small) member read in chunks, refused past MEMBER_MAX_BYTES once decompressed."""
    out, n = io.BytesIO(), 0
    while block := fh.read(1 << 24):
        n += len(block)
        if n > MEMBER_MAX_BYTES:
            raise Refused(f"{name}: decompresses to more than {MEMBER_MAX_BYTES} bytes: refused")
        out.write(block)
    return out.getvalue()


def _mtx_header(fh, label: str) -> tuple[int, int, int]:
    """(features, barcodes, entries) of a MatrixMarket coordinate file, read from its banner and size line."""
    head = fh.readline(RAW_MAX_LINE).decode("ascii", errors="replace").strip().lower().split()
    if (len(head) != 5 or head[0] != "%%matrixmarket" or head[1:3] != ["matrix", "coordinate"]
            or head[3] not in ("integer", "real") or head[4] != "general"):
        raise Refused(f"{label}: the matrix is not a MatrixMarket coordinate integer or real general file: refused")
    skipped = 0
    while True:
        line = fh.readline(RAW_MAX_LINE)
        if not line:
            raise Refused(f"{label}: the matrix has no size line: refused")
        if line.startswith(b"%") or not line.strip():
            # stage 8 (the first reviewer of ec86a1b, lesser point 5): a bounded comment block
            skipped += 1
            if skipped > RAW_MAX_COMMENT_LINES:
                raise Refused(f"{label}: more than {RAW_MAX_COMMENT_LINES} comment or empty lines before the size "
                              "line: refused")
            continue
        parts = line.split()
        if len(parts) != 3 or not all(p.isdigit() for p in parts):
            raise Refused(f"{label}: the matrix size line is not three counts: refused")
        return int(parts[0]), int(parts[1]), int(parts[2])


def _parse_entries(block: bytes, label: str, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(feature, barcode, value) arrays, 0-based, of whole lines; refused for a malformed line, an index out of
    range or a negative or non-finite value."""
    if not block.strip():
        return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0)
    # stage 8 (the first reviewer of ec86a1b, lesser point 1): read with no column names, so a line of four fields is
    # refused rather than its first field taken as an index; ragged rows already raise ParserError
    try:
        df = pd.read_csv(io.BytesIO(block), sep=r"\s+", header=None, engine="c")
    except (ValueError, pd.errors.ParserError):
        raise Refused(f"{label}: a matrix line is not 'feature barcode value': refused") from None
    if df.shape[1] != 3:
        raise Refused(f"{label}: a matrix line is not 'feature barcode value': refused")
    df.columns = ["r", "c", "v"]
    # the two indices must parse as integers (a '1.5' or a missing field is refused, never truncated or filled)
    if df.isna().any().any() or not (pd.api.types.is_integer_dtype(df["r"]) and pd.api.types.is_integer_dtype(df["c"])
            and (pd.api.types.is_integer_dtype(df["v"]) or pd.api.types.is_float_dtype(df["v"]))):
        raise Refused(f"{label}: a matrix line is not 'feature barcode value': refused")
    df = df.astype({"r": np.int64, "c": np.int64, "v": np.float64})
    r, c, v = df["r"].to_numpy() - 1, df["c"].to_numpy() - 1, df["v"].to_numpy()
    if df.isna().any().any() or (r < 0).any() or (r >= shape[0]).any() or (c < 0).any() or (c >= shape[1]).any() \
            or (v < 0).any() or not np.isfinite(v).all():
        raise Refused(f"{label}: a matrix entry out of range, negative or not finite: refused")
    return r, c, v


def _mtx_entries(fh, label: str, shape: tuple[int, int], nnz: int):
    """The entries after the size line, block by block; refused for a line longer than RAW_MAX_LINE and for more or
    fewer entries than the header states."""
    carry, seen = b"", 0
    while True:
        block = fh.read(RAW_CHUNK_BYTES)
        if not block:
            break
        block = carry + block
        cut = block.rfind(b"\n")
        carry = block[cut + 1:]
        if len(carry) > RAW_MAX_LINE:
            raise Refused(f"{label}: a matrix line longer than {RAW_MAX_LINE} bytes: refused")
        if cut < 0:
            continue
        r, c, v = _parse_entries(block[:cut + 1], label, shape)
        seen += len(r)
        if seen > nnz:
            raise Refused(f"{label}: more matrix entries than its header states ({nnz}): refused")
        yield r, c, v
    if carry.strip():
        r, c, v = _parse_entries(carry, label, shape)
        seen += len(r)
        yield r, c, v
    if seen != nnz:
        raise Refused(f"{label}: {seen} matrix entries, not the {nnz} its header states: refused")


def _raw_from(stream, label: str):
    """The called cells of one unfiltered library: (cells x features CSR of the called barcodes only, symbols, feature
    types, called barcodes, the cell call record). stream(kind) opens a fresh decompressed binary stream of the
    matrix, features or barcodes member."""
    def text(kind):
        fh = stream(kind)
        try:
            return _read_capped(fh, f"{label} {kind}").decode("utf-8", errors="replace").splitlines()
        finally:
            fh.close()

    rows = [r.split("\t") for r in text("features") if r.strip()]
    symbols = [r[1] if len(r) > 1 else r[0] for r in rows]
    types = [r[2] if len(r) > 2 else "Gene Expression" for r in rows]
    barcodes = [b.split("\t")[0] for b in text("barcodes") if b.strip()]
    if len(barcodes) > cm.RAW_MAX_BARCODES:
        raise Refused(f"{label}: {len(barcodes)} barcodes in the file, more than {cm.RAW_MAX_BARCODES}: refused")
    gex = np.array([("gene" in t.lower() and "expression" in t.lower()) or t == "" for t in types])
    fh = stream("matrix")
    try:
        n_f, n_b, nnz = _mtx_header(fh, label)
        if (n_f, n_b) != (len(symbols), len(barcodes)):
            raise Refused(f"{label}: the matrix header ({n_f} x {n_b}) does not match {len(symbols)} features x "
                          f"{len(barcodes)} barcodes: refused")
        if nnz > cm.RAW_MAX_ENTRIES:
            raise Refused(f"{label}: {nnz} matrix entries, more than {cm.RAW_MAX_ENTRIES}: refused")
        # pass 1: each barcode's total over gene-expression features, and whether it holds any count at all
        totals, anyc = np.zeros(n_b), np.zeros(n_b, bool)
        for r, c, v in _mtx_entries(fh, label, (n_f, n_b), nnz):
            g = gex[r]
            totals += np.bincount(c[g], weights=v[g], minlength=n_b)
            anyc[c[v > 0]] = True
    finally:
        fh.close()
    called, rec = cm.raw_cell_call(totals)
    n = rec["cells_called"]
    ct = totals[called]
    rec = {"barcodes_in_file": n_b, "barcodes_with_any_count": int(anyc.sum()), **rec,
           "umi_share_called": float(totals[called].sum() / totals.sum()) if totals.sum() > 0 else 0.0,
           "features": n_f, "entries": nnz,
           # stage 8 (the second reviewer of ec86a1b, lesser point 1): records only, read by no guard
           **{f"called_total_p{q}": (float(np.percentile(ct, q)) if len(ct) else None) for q in (5, 25, 50)},
           "called_under_twice_cutoff": int((ct < 2 * rec["cutoff"]).sum())}
    if n < cm.RAW_MIN_CELLS:
        raise Refused(f"{label}: the fixed cell call keeps {n} cells, fewer than {cm.RAW_MIN_CELLS}: refused")
    if n > cm.MAX_BARCODES:
        raise Refused(f"{label}: the fixed cell call keeps {n} cells, more than the cap {cm.MAX_BARCODES}: refused")
    # pass 2: the called barcodes' entries only
    col = np.full(n_b, -1, np.int64)
    idx = np.flatnonzero(called)
    col[idx] = np.arange(len(idx))
    parts_r, parts_c, parts_v = [], [], []
    fh = stream("matrix")
    try:
        _mtx_header(fh, label)
        for r, c, v in _mtx_entries(fh, label, (n_f, n_b), nnz):
            k = called[c]
            parts_r.append(col[c[k]])
            parts_c.append(r[k])
            parts_v.append(v[k])
    finally:
        fh.close()
    M = sparse.coo_matrix((np.concatenate(parts_v) if parts_v else np.zeros(0),
                           (np.concatenate(parts_r) if parts_r else np.zeros(0, np.int64),
                            np.concatenate(parts_c) if parts_c else np.zeros(0, np.int64))),
                          shape=(len(idx), n_f)).tocsr()
    return M, symbols, types, [barcodes[i] for i in idx], rec


def read_raw_sample(paths: list[Path], label: str):
    """A pinned sample's declared unfiltered matrix (one archive alone, or one loose 10x triple) through the fixed cell
    call: (CSR of the called cells, symbols, types, called barcodes, record). An archive is read member by member in
    memory with archive_members (raw_ok: only the unfiltered-matrix deny term is lifted) and must hold exactly one
    triple in one folder and nothing else but the macOS companion members archive_members skips. A file above
    RAW_MAX_FILE_BYTES, a 10x h5, a matrix of more entries than
    RAW_MAX_ENTRIES or more barcodes than RAW_MAX_BARCODES, fewer than RAW_MIN_CELLS or more than MAX_BARCODES called
    cells: refused."""
    paths = [Path(p) for p in paths]
    for p in paths:
        if p.stat().st_size > cm.RAW_MAX_FILE_BYTES:
            raise Refused(f"{label}: {p.name} is larger than {cm.RAW_MAX_FILE_BYTES} bytes: refused")
    if any(p.name.endswith(".h5") for p in paths):
        raise Refused(f"{label}: the unfiltered intake reads a MatrixMarket triple only, not an h5 file: refused")
    tars = [p for p in paths if p.name.endswith((".tar.gz", ".tgz", ".tar"))]
    if tars:
        if len(paths) != 1:
            raise Refused(f"{label}: an unfiltered archive comes alone, not with {len(paths) - 1} other files")
        skipped: list[str] = []
        with tarfile.open(paths[0], "r:*") as t:
            by = _one_folder_triple(archive_members(t, label, raw_ok=True, skipped=skipped), label)

            def stream(kind):
                f = t.extractfile(by[kind][2])
                if f is None:
                    raise Refused(f"{label}: member {'/'.join(by[kind][0])!r} cannot be read")
                return _decompressed(f, by[kind][0][-1])

            res = _raw_from(stream, label)
        if skipped:   # recorded only when present, so a plain archive's record is unchanged
            res[-1]["mac_companions_skipped"] = len(skipped)
        return res
    items = [(p.name, cm.triple_kind(p.name), p) for p in paths]
    if any(it[1] is None for it in items):
        raise Refused(f"{label}: a file of the unfiltered triple is not a matrix, barcodes or features file: refused")
    by = _one_triple(items, label)

    def stream_file(kind):
        p = by[kind][2]
        return gzip.open(p, "rb") if p.name.endswith(".gz") else open(p, "rb")   # noqa: SIM115 (closed by _raw_from)

    return _raw_from(stream_file, label)


def to_panel(M, symbols: list[str], types: list[str], genes: list[str]) -> tuple[sparse.csr_matrix, np.ndarray, dict]:
    """Counts on the panel genes in the reference's order (absent columns zero), the presence mask and a record.
    Only gene-expression features; exact symbol first, then the committed alias table; duplicates summed."""
    gex = np.array([("gene" in t.lower() and "expression" in t.lower()) or t == "" for t in types])
    pos: dict[str, list[int]] = {}
    for i, s in enumerate(symbols):
        if gex[i]:
            pos.setdefault(s, []).append(i)
    rows, cols, how = [], [], {"exact": 0, "alias": 0, "duplicate_symbols": 0}
    for j, g in enumerate(genes):
        idx = pos.get(g)
        if idx:
            how["exact"] += 1
        elif g in cm.ALIASES and pos.get(cm.ALIASES[g]):
            idx = pos[cm.ALIASES[g]]
            how["alias"] += 1
        if idx:
            how["duplicate_symbols"] += int(len(idx) > 1)
            rows += idx
            cols += [j] * len(idx)
    Map = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(symbols), len(genes)))
    present = np.zeros(len(genes), bool)
    present[cols] = True
    C = sparse.csr_matrix(M @ Map).astype(np.int64)
    return C, present, {**how, "matched": int(present.sum()), "of": len(genes)}


def cell_filter(M, types) -> np.ndarray:
    gex = np.array([("gene" in t.lower() and "expression" in t.lower()) or t == "" for t in types])
    return np.asarray((M[:, np.flatnonzero(gex)] > 0).sum(1)).ravel() >= cm.MIN_GENES_PER_CELL


# ---------------------------------------------------------------------------------------------------------------------
# The log ratios (pure numpy; synthetic tests run these)


def _nanmean_rows(A: np.ndarray) -> np.ndarray:
    held = np.isfinite(A)
    n = held.sum(0)
    return np.where(n > 0, np.where(held, A, 0.0).sum(0) / np.maximum(n, 1), np.nan)


def _sum(C, rows) -> np.ndarray:
    return np.asarray(C[rows].sum(0)).ravel().astype(np.float64)


def log_ratio(x_a, x_b, present) -> np.ndarray:
    """log((x_a+a)/N_a) - log((x_b+a)/N_b) on present genes, clipped; NaN for an absent gene and for a gene with
    fewer than MIN_UMI in either arm (not measured: it holds no value in any mean and ends as factor 1 if no group
    or study holds it)."""
    x_a, x_b = np.asarray(x_a, float), np.asarray(x_b, float)
    Na, Nb = x_a[present].sum(), x_b[present].sum()
    out = np.full(len(x_a), np.nan)
    if Na <= 0 or Nb <= 0:
        return out
    L = np.log((x_a + cm.PSEUDOCOUNT) / Na) - np.log((x_b + cm.PSEUDOCOUNT) / Nb)
    L = np.where(np.minimum(x_a, x_b) < cm.MIN_UMI, np.nan, np.clip(L, *cm.LOG_CLIP))
    out[present] = L[present]
    return out


def chimera_log_ratio(x: dict[str, np.ndarray], present) -> np.ndarray:
    """[ko_mut vs ko_host] - [wt_mut vs wt_host], each with the pseudocount; NaN where any arm has < MIN_UMI."""
    a = log_ratio(x["ko_mut"], x["ko_host"], present)
    b = log_ratio(x["wt_mut"], x["wt_host"], present)
    low = np.min(np.vstack([x[k] for k in CHIMERA_ARMS]), 0) < cm.MIN_UMI
    L = np.where(low, np.nan, np.clip(a - b, *cm.LOG_CLIP))
    out = np.full(len(L), np.nan)
    out[present] = L[present]
    return out


def spearman_guard(profile, ref: dict, g: str, present) -> tuple[bool, dict]:
    """The control group's mean-proportion profile must match reference group g best, ahead of the next by the
    fixed margin (computed on control cells and the reference only: not a try)."""
    r = rankdata(np.asarray(profile)[present])
    rho = {}
    for h in ref["groups"]:
        q = rankdata(np.asarray(ref["mu"][h])[present])
        rho[h] = float(np.corrcoef(r, q)[0, 1]) if r.std() > 0 and q.std() > 0 else float("nan")
    order = sorted((v, h) for h, v in rho.items() if np.isfinite(v))[::-1]
    if not order or order[0][1] != g:
        return False, {"best": order[0][1] if order else None, "rho_own": rho.get(g)}
    margin = order[0][0] - (order[1][0] if len(order) > 1 else -1.0)
    return margin >= cm.SPEARMAN_MARGIN, {"best": g, "rho_own": rho[g], "margin": margin}


def unit_log_ratios(C, present, obs: pd.DataFrame, ref: dict, d: float, chimera: bool = False) -> tuple[dict, dict]:
    """One knockout gene's table: {group: L (P,) averaged over passing stages, NaN absent}, and its QC record.

    obs columns: arm (ko/control, or ko_mut/ko_host/wt_mut/wt_host for chimeras), stage (float), sample, call
    (reference group or ''), conf; optional author_group (the committed dictionary's group, or 'none')."""
    arms = CHIMERA_ARMS if chimera else PLAIN_ARMS
    if set(obs["arm"]) - set(arms):
        raise Refused(f"arms {sorted(set(obs['arm']) - set(arms))} are not of this design ({arms})")
    if chimera and set(arms) - set(obs["arm"]):
        raise Refused("a chimera needs knockout and wild-type chimeras (mut and host each): no second difference, "
                      "refused")
    C = sparse.csr_matrix(C)
    conf_ok = obs["conf"].to_numpy(float) >= cm.CONFIDENT
    per_stage: dict[float, dict] = {}
    qc = {"stages": {}, "groups": {}}
    for s in sorted(set(obs["stage"])):
        cm.stage_ok(f"E{s:g}")
        at = obs["stage"].to_numpy() == s
        if not all((at & (obs["arm"] == a).to_numpy()).any() for a in arms):
            qc["stages"][f"E{s:g}"] = {"skipped": "not every arm has cells"}
            continue
        tot = {a: int((at & (obs["arm"] == a).to_numpy()).sum()) for a in arms}
        cshare = {a: float((at & (obs["arm"] == a).to_numpy() & conf_ok).sum() / tot[a]) for a in arms}
        rec = {"cells": tot, "confident_share": cshare, "groups": {}}
        qc["stages"][f"E{s:g}"] = rec
        if min(cshare.values()) < cm.MIN_CONF_SHARE:
            rec["skipped"] = "fewer than half of an arm's cells called confidently"
            continue
        for g in ref["groups"]:
            sel = {a: np.flatnonzero(at & (obs["arm"] == a).to_numpy() & conf_ok & (obs["call"] == g).to_numpy())
                   for a in arms}
            n = {a: len(v) for a, v in sel.items()}
            share = {a: n[a] / tot[a] for a in arms}
            gr = {"cells": n, "call_share": share}
            rec["groups"][g] = gr
            if min(n.values()) < cm.MIN_CELLS:
                gr["skipped"] = f"fewer than {cm.MIN_CELLS} confident cells in an arm"
                continue
            if max(share.values()) > cm.SHARE_RATIO_MAX * min(share.values()):
                gr["skipped"] = "call share differs across arms by more than 2x (relabelling)"
                continue
            ctl = np.concatenate([sel[a] for a in arms if a in CONTROL_ARMS])
            Pc = C[ctl].toarray().astype(float)
            prof = (Pc / np.maximum(Pc.sum(1, keepdims=True), 1)).mean(0)
            ok, info = spearman_guard(prof, ref, g, present)
            gr["profile_guard"] = info
            if not ok:
                gr["skipped"] = "control profile does not match its reference group by the margin"
                continue
            if "author_group" in obs.columns:
                agree = float((obs["author_group"].to_numpy()[ctl] == g).mean())
                gr["author_agreement"] = agree
                if agree < cm.AUTHOR_AGREE:
                    gr["skipped"] = "author-label agreement under one half"
                    continue
            x = {a: _sum(C, sel[a]) for a in arms}
            L = chimera_log_ratio(x, present) if chimera else log_ratio(x["ko"], x["control"], present)
            m = min(len(set(obs["sample"].to_numpy()[sel[a]].astype(str))) for a in arms)
            k = cm.shrink_k(list(n.values()), m, d)
            gr.update(k=k, m=m, passed=True)
            per_stage.setdefault(s, {})[g] = k * L
    out = {}
    for g in ref["groups"]:
        vals = [per_stage[s][g] for s in per_stage if g in per_stage[s]]
        if vals:
            out[g] = _nanmean_rows(np.vstack(vals))
            qc["groups"][g] = {"stages": [f"E{s:g}" for s in per_stage if g in per_stage[s]]}
    qc["status"] = "ok" if out else "left_out"
    if not out:
        qc["reason"] = "no group passes the guards (fail closed: left out, no pooled fallback)"
    return out, qc


# ---------------------------------------------------------------------------------------------------------------------
# Classification (DGX only; synthetic tests pass their own calls)


def classify(wt_path: Path, C, present, genes: list[str], names: list[str], threads: int = 8):
    """Group calls and confidences for outside cells by transfer from the WT E8.75 reference's lineage groups."""
    import anndata as ad

    from vec import harmonise as hm
    from vec.prep import load_stage, log1p_cp10k

    pres = [g for g, p in zip(genes, present, strict=True) if p]
    ref = load_stage(wt_path, genes=genes, drop_duplicates=True)
    ref = ref[:, pres].copy()
    Cq = sparse.csr_matrix(C)[:, np.flatnonzero(present)]
    q = ad.AnnData(X=log1p_cp10k(Cq),
                   obs=pd.DataFrame(index=names), var=pd.DataFrame(index=pres))
    q.layers["counts"] = Cq.astype(np.int32)
    qd = float(np.median(np.asarray(Cq.sum(1)).ravel()))
    rd = hm.median_depth(ref, pres)
    depth = qd if (rd is not None and rd > qd * 1.05) else None
    model = hm.fit_transfer(ref, "group", genes=pres, classifier="logreg", n_pcs=50, C=1.0, clip=10.0,
                            threshold=cm.CONFIDENT, min_class_cells=10, exclude_labels=("Unknown",), depth=depth,
                            depth_tol=0.05, seed=0, n_jobs=threads)
    df = hm.predict(model, q, depth_draws=5, depth_tol=0.05, seed=0, proba=False)
    return df["label_raw"].astype(str).to_numpy(), df["label_conf"].to_numpy(float), {
        "reference_depth": rd, "query_depth": qd, "refit_depth": depth, "genes_used": len(pres)}


# ---------------------------------------------------------------------------------------------------------------------
# Orchestration


def write_unit(out_dir: Path, acc: str, gene: str, Lg: dict, qc: dict, ref: dict, meta: dict) -> dict:
    name = f"{acc}__{gene}".replace("+", "_and_").replace("/", "_")
    out_dir = cm.outside_repo(out_dir)
    unit = {"accession": acc, "gene": gene, "genes": ref["genes"], "refstats_sha256": ref["file_sha256"],
            "L": {g: [None if not np.isfinite(v) else float(v) for v in L] for g, L in Lg.items()}, **meta}
    up = out_dir / f"{name}.unit.json"
    cm.dump_json(unit, up)
    rec = {"unit_file": str(up), "unit_sha256": cm.sha256_file(up), "status": qc["status"]}
    if qc["status"] == "ok":
        phi, used = tb.finish_pooled(Lg, ref)
        rec["tables"] = tb.write_table(out_dir / "tables" / name, name, ref["genes"], phi,
                                       tb.finish_per_group(Lg, ref), {**meta, "groups_pooled": used})
    return rec


def run(acc: str, dest: Path, wt: Path, refstats: Path, out: Path, datasets: Path | None = None,
        threads: int = 8, classifier=None) -> dict:
    ds = cm.dataset(acc, datasets)
    if cm.is_bulk(ds):   # stage 9: a bulk dataset goes through bulk.run_bulk; the single-cell path below is unchanged
        import bulk

        return bulk.run_bulk(acc, dest, wt, refstats, out, datasets)
    floor = cm.qc_floors(datasets)["min_confident_share_per_sample"]
    if ds.get("cell_labels"):
        raise Refused(f"{acc}: per-cell genotype designs (cell-level labels inside one sample) are excluded in this "
                      "stage")
    if not ds.get("pinned"):
        raise Refused(f"{acc}: the pinned block is empty: nothing is processed")
    dest, out = cm.outside_repo(dest), cm.outside_repo(out)
    # stage 5: no stale output. This accession's earlier outputs go before any check, so a run that refuses leaves
    # nothing that combine could read
    adir = out / acc
    if adir.is_symlink():
        raise Refused(f"{adir}: a symlinked output folder: remove it by hand")
    if adir.exists():
        shutil.rmtree(adir)
    pl = fe._load_plan(acc, dest, datasets)
    plan_sha = cm.sha256_file(dest / "meta" / acc / "plan.json")
    rdir = dest / "raw" / acc
    mp = rdir / "MANIFEST.json"
    for d in (dest / "raw", rdir, mp):
        if d.is_symlink():
            raise Refused(f"{d}: a symlink in the raw folder: refused (raw/ and MANIFEST.json are written by fetch "
                          "only)")
    if not mp.is_file():
        raise Refused(f"{mp}: no manifest; fetch first")
    man = json.loads(mp.read_text())
    if man.get("plan_sha256") != plan_sha:
        raise Refused(f"{acc}: the manifest was not fetched from this plan")
    cur = cm.pins_sha256(ds)
    if (man.get("pins_sha256") != cur or man.get("pinned") != ds["pinned"]
            or man.get("pinned_series") != ds["pinned_series"]
            or (man.get("pinned_series_files") or []) != (ds.get("pinned_series_files") or [])):
        raise Refused(f"{acc}: the samples, arms, genes, stages, series or series files are not exactly as pinned when "
                      "the files were fetched: fetch again")
    # 1. the pins against the stored metadata (samples and series), the series membership of every pinned sample, and
    #    the screen again on every pinned sample, before any count is read; stage 6: the pinned series files, their
    #    links derived again from the hashed series fields
    pinned = fe.verify_pins(ds, pl)
    sunits = fe.verify_series_units(ds, pl)
    for su in sunits:
        sdir = rdir / su["series"]
        if sdir.is_symlink():
            raise Refused(f"{su['series']}: a symlinked series folder: refused")
        names = {f["name"] for u in sunits if u["series"] == su["series"] for f in u["files"]}
        for q in (sorted(sdir.iterdir()) if sdir.is_dir() else []):
            if q.is_symlink():
                raise Refused(f"{su['series']}/{q.name}: a symlink in a pinned series folder: refused")
            if q.name not in names or f"{su['series']}/{q.name}" not in man["files"]:
                raise Refused(f"{su['series']}/{q.name}: a file in a pinned series folder that is not a pinned series "
                              "file with a manifest record: refused (fetch again)")
        for f in su["files"]:
            rec = man["files"].get(f"{su['series']}/{f['name']}")
            if rec is None:
                raise Refused(f"{su['series']}/{f['name']}: no manifest record: refused (fetch again)")
            if (rec.get("series_metadata_sha256") != su["series_metadata_sha256"] or rec.get("url") != f["url"]
                    or rec.get("pins_sha256") != cur):
                raise Refused(f"{su['series']}/{f['name']}: fetched under another series text, url or pins: refused "
                              "(fetch again)")
    # 1b. every file of a pinned sample is bound to the pin it was fetched under (stage 4): a manifest record per file
    #     with this pin's metadata_sha256, this url and these pins; no other file in a pinned sample's folder
    for s in pinned:
        gsm = s["accession"]
        names = {f["name"] for f in s["files"]}
        sdir = dest / "raw" / acc / gsm
        if sdir.is_symlink():
            raise Refused(f"{gsm}: a symlinked sample folder: refused")
        on_disk = sorted(q.name for q in sdir.iterdir()) if sdir.is_dir() else []
        for name in on_disk:
            if (sdir / name).is_symlink():
                raise Refused(f"{gsm}/{name}: a symlink in a pinned sample's folder: refused")
            if name not in names or f"{gsm}/{name}" not in man["files"]:
                raise Refused(f"{gsm}/{name}: a file in a pinned sample's folder that is not one of its pinned files "
                              "with a manifest record: refused (fetch again)")
        for f in s["files"]:
            rec = man["files"].get(f"{gsm}/{f['name']}")
            if rec is None:
                raise Refused(f"{gsm}/{f['name']}: no manifest record: refused (fetch again)")
            if rec.get("metadata_sha256") != s["metadata_sha256"]:
                raise Refused(f"{gsm}/{f['name']}: fetched under metadata_sha256 {rec.get('metadata_sha256')}, but "
                              f"the sample is pinned now with {s['metadata_sha256']}: refused (fetch again)")
            if rec.get("url") != f["url"] or rec.get("pins_sha256") != cur:
                raise Refused(f"{gsm}/{f['name']}: fetched from another url or under other pins: refused (fetch "
                              "again)")
    units = cm.pinned_units(ds["pinned"])
    for s in pinned:
        cm.stage_ok(f"E{s['stage']:g}")
    screen = [s["accession"] for s in pinned]
    ref = tb.load_ref_stats(refstats)
    if ref["wt_sha256"] != cm.sha256_file(wt):
        raise Refused("the WT file is not the one the reference statistics were computed from")
    genes = ref["genes"]
    # 2. verify the pinned files, then read counts on the panel (no other sample's file is opened)
    blocks, obs_parts, presence, panel_rec, dropped, count_sha = [], [], [], {}, {}, {}
    by_gsm = {s["accession"]: s for s in pinned}
    # the library label each pinned sample's cells carry (m and the confidence floor count labels): its own accession,
    # or for a series file used as one unit (split none) with several samples, one label for the whole file
    label_of = {g: g for g in by_gsm}
    unit_recs = []

    def checked(key: str, p: Path) -> Path:
        if p.is_symlink():
            raise Refused(f"{key}: a symlink: refused")
        digest = cm.sha256_file(p) if p.is_file() else None
        if key not in man["files"] or digest is None or digest != man["files"][key]["sha256"]:
            raise Refused(f"{key}: missing or its sha256 differs from the manifest")
        count_sha[key] = digest
        return p

    def add(label: str, pin: str, M, symbols, types, bcs, folder: str = "") -> None:
        s = by_gsm[pin]
        prefix = f"{folder}/" if folder else ""
        keep = cell_filter(M, types)
        dropped[label] = dropped.get(label, 0) + int((~keep).sum())
        Cp, pres, prec = to_panel(M[np.flatnonzero(keep)], symbols, types, genes)
        panel_rec[label + (f"/{folder}" if folder else "")] = prec
        presence.append(pres)
        blocks.append(Cp)
        obs_parts.append(pd.DataFrame({"arm": s["arm"], "gene": ";".join(s["genes"]), "stage": s["stage"],
                                       "sample": label, "pin": pin},
                                      index=[f"{label}:{prefix}{b}" for b, k in zip(bcs, keep, strict=True) if k]))

    raw_rec: dict[str, dict] = {}
    for s in pinned:
        gsm = s["accession"]
        if not s["files"]:
            continue   # stage 6: its cells come from a pinned series file (check_pins, check_series_file_pins)
        files = [checked(f"{gsm}/{f['name']}", dest / "raw" / acc / gsm / f["name"]) for f in s["files"]]
        raw = {bool(f.get("raw")) for f in s["files"]}
        if raw == {True}:
            # stage 7: a declared unfiltered matrix: the one fixed cell call, the same for every sample of every arm;
            # then the ordinary per-cell rule (add) on the called cells, as for any filtered library
            M, symbols, types, bcs, raw_rec[gsm] = read_raw_sample(files, gsm)
        elif raw == {False}:
            M, symbols, types, bcs = read_sample(files)
            if M.shape[0] > cm.MAX_BARCODES:
                raise Refused(f"{gsm}: {M.shape[0]} barcodes: an unfiltered matrix, refused")
        else:
            raise Refused(f"{gsm}: an unfiltered matrix mixed with other files: refused")
        add(gsm, gsm, M, symbols, types, bcs)
    for su in sunits:
        paths = [checked(f"{su['series']}/{f['name']}", rdir / su["series"] / f["name"]) for f in su["files"]]
        one = su["split"]["by"] == "none" and len(su["samples"]) > 1
        label = f"series:{su['series']}/{su['files'][0]['name']}" if one else None
        libs = read_series_unit(paths, {**su, "files": [f["name"] for f in su["files"]]})
        rec = {"series": su["series"], "files": [f["name"] for f in su["files"]], "samples": su["samples"],
               "split": su["split"], "libraries": []}
        for gsm, folder, kept, (M, symbols, types, bcs) in libs:
            lab = label or gsm
            for g in su["samples"]:
                if one:
                    label_of[g] = lab
            add(lab, gsm, M, symbols, types, bcs, folder=folder)
            rec["libraries"].append({"label": lab, "folder": folder, "barcodes": len(kept)})
        unit_recs.append(rec)
    # stage 8 (the first reviewer of ec86a1b, blocker 2; the second reviewer's lesser point 2): the second guard
    # behind fetch.check_unit_forms. Both arms of a unit are read the same way, before any group call is made
    for u in units:
        members = u["ko"] + u["control"]
        n_raw = sum(g in raw_rec for g in members)
        if 0 < n_raw < len(members):
            raise Refused(f"{u['gene']}: an unfiltered and a filtered library in one unit: refused")
    present = np.logical_and.reduce(presence)
    C = sparse.vstack(blocks).tocsr()
    obs = pd.concat(obs_parts)
    qc = {"accession": acc, "tier": ds["tier"], "study": ds["study"], "d": ds["d"], "screened_samples": screen,
          "pinned": ds["pinned"], "pinned_series": ds["pinned_series"], "pins_sha256": cur,
          "plan_sha256": plan_sha, "manifest_sha256": cm.sha256_file(mp), "count_files_sha256": count_sha,
          "screen_record": pl["screen"], "cells_dropped_low_genes": dropped, "panel_match": panel_rec,
          "panel_genes_present": int(present.sum()), "panel_genes": len(genes),
          "created": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "git": _git_head()}
    if sunits:
        qc["series_file_units"] = unit_recs
        qc["pinned_series_files"] = ds["pinned_series_files"]
    if cm.control_series(ds):
        qc["control_series"] = cm.control_series(ds)
        qc["sister_series_controls"] = sorted(g for g, s in by_gsm.items() if s.get("sister"))
    if ds.get("raw_matrices"):
        # stage 7: per raw sample, barcodes in the file, barcodes with any count, m, the cutoff, cells called and the
        # share of all gene-expression UMIs the called cells hold; the rule's constants once
        qc["raw_cell_call"] = {"rule": cm.raw_call_rule(), "raw_matrices": ds["raw_matrices"], "samples": raw_rec}
    if present.mean() < cm.MIN_PANEL_MATCH:
        qc.update(status="failed", reason=f"only {present.sum()} of {len(genes)} panel genes matched")
        cm.dump_json(qc, adir / "QC.json")
        return qc
    # 3. group calls
    calls, conf, crec = (classifier or classify)(wt, C, present, genes, list(obs.index), threads)
    obs["call"], obs["conf"] = calls, conf
    qc["classifier"] = crec
    # 3b. the transfer has never run on real cells: each sample's confident share goes into QC.json, and a table
    # with a sample under the floor fixed in datasets.yaml is refused
    by_sample = {}
    for gsm, sub in obs.groupby("sample"):
        share = float((sub["conf"].to_numpy(float) >= cm.CONFIDENT).mean()) if len(sub) else 0.0
        by_sample[str(gsm)] = {"cells": int(len(sub)), "confident_share": share, "unconfident_share": 1.0 - share}
    qc["confidence_by_sample"] = by_sample
    qc["min_confident_share_per_sample"] = floor
    # 4. one table per knockout gene
    qc["units"] = {}
    for u in units:
        low = sorted(g for g in u["ko"] + u["control"]
                     if by_sample.get(label_of[g], {}).get("confident_share", 0.0) < floor)
        if low:
            qc["units"][u["gene"]] = {"status": "refused_low_confidence", "samples_under_floor": low,
                                      "reason": f"confident share under {floor} in {low}: no table written"}
            continue
        rows = obs["pin"].isin(u["ko"] + u["control"]).to_numpy()
        sub = obs[rows].copy()
        sub["arm"] = np.where(sub["pin"].isin(u["ko"]), "ko", "control")
        # stage 6: a control from the sister series fixes d at the weaker value (t3x_common.unit_design_factor)
        sister = sorted(g for g in u["control"] if by_gsm[g].get("sister"))
        d_u = cm.unit_design_factor(ds, bool(sister))
        Lg, uqc = unit_log_ratios(C[np.flatnonzero(rows)], present, sub, ref, d_u)
        cells = sub.groupby(["stage", "arm"]).size()
        uqc["cells_per_stage_and_arm"] = {f"E{s:g}:{a}": int(v) for (s, a), v in cells.items()}
        meta = {"tier": ds["tier"], "study": ds["study"], "d": d_u, "ko_samples": u["ko"],
                "control_samples": u["control"], "present": present.tolist(), "git": qc["git"]}
        if ds.get("raw_matrices"):
            # stage 7: the one flag downstream that says a library was called from an unfiltered matrix
            uqc["raw_called_samples"] = sorted(g for g in u["ko"] + u["control"] if g in raw_rec)
            uqc["raw_called"] = bool(uqc["raw_called_samples"])
            if uqc["raw_called"]:
                uqc["raw_call_arms"] = raw_call_arms(raw_rec, u["ko"], u["control"])
        if cm.control_series(ds):
            uqc["d"] = meta["d"]
            uqc["control_source"] = meta["control_source"] = (
                f"sister series {cm.control_series(ds)['accession']} (control_series): {sister}" if sister
                else "the dataset's own series")
        uqc.update(write_unit(adir, acc, u["gene"], Lg, uqc, ref, meta))
        qc["units"][u["gene"]] = uqc
    qc["status"] = "ok" if any(v["status"] == "ok" for v in qc["units"].values()) else "left_out"
    cm.dump_json(qc, adir / "QC.json")
    return qc


RAW_ARM_KEYS = ("cutoff", "cells_called", "umi_share_called")


def raw_call_arms(raw_rec: dict[str, dict], ko: list[str], control: list[str]) -> dict:
    """Stage 8 (the second reviewer of ec86a1b, lesser point 1): for a unit read through the fixed cell call, the
    per-arm medians over its samples of the cutoff, the cells called and the share of gene-expression UMIs in called
    cells, and their knockout-to-control ratios (None where the control median is 0). A record only: no guard and no
    threshold reads it."""
    med = {arm: {k: float(np.median([raw_rec[g][k] for g in gs])) for k in RAW_ARM_KEYS}
           for arm, gs in (("ko", ko), ("control", control))}
    ratio = {k: (med["ko"][k] / med["control"][k] if med["control"][k] else None) for k in RAW_ARM_KEYS}
    return {"median_ko": med["ko"], "median_control": med["control"], "ko_to_control": ratio,
            "note": "a record only: no guard or threshold reads it"}


def _git_head() -> dict:
    import subprocess

    def g(*a):
        try:
            return subprocess.run(["git", "-C", str(cm.REPO), *a], capture_output=True, text=True, timeout=30,
                                  env={"GIT_OPTIONAL_LOCKS": "0", "PATH": "/usr/bin:/bin"}).stdout.strip()
        except Exception:  # noqa: BLE001
            return "unknown"

    return {"commit": g("rev-parse", "HEAD"), "dirty": bool(g("status", "--porcelain", "--untracked-files=no"))}


def reference(wt: Path, out: Path) -> dict:
    from vec.prep import load_stage

    out = cm.outside_repo(out)
    a = load_stage(wt, drop_duplicates=True)
    st = tb.ref_stats_from_counts(a.layers["counts"], a.obs["group"].astype(str).to_numpy(), list(a.var_names))
    tb.save_ref_stats(st, out, cm.sha256_file(wt))
    return st


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="step", required=True)
    r = sub.add_parser("reference")
    r.add_argument("--wt", type=Path, required=True)
    r.add_argument("--out", type=Path, required=True)
    u = sub.add_parser("run")
    u.add_argument("accession")
    u.add_argument("--dest", type=Path, required=True)
    u.add_argument("--wt", type=Path, required=True)
    u.add_argument("--refstats", type=Path, required=True)
    u.add_argument("--out", type=Path, required=True)
    u.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    try:
        if a.step == "reference":
            st = reference(a.wt, a.out)
            print(f"reference: {len(st['genes'])} genes, {len(st['groups'])} groups -> {a.out}")
            return 0
        qc = run(a.accession, a.dest, a.wt, a.refstats, a.out, threads=a.threads)
    except Refused as e:
        print(f"refused: {e}")
        return 2
    print(f"{a.accession}: status {qc['status']}; panel genes present {qc['panel_genes_present']}/{qc['panel_genes']}")
    for gsm, v in (qc.get("confidence_by_sample") or {}).items():
        print(f"  {gsm}: {v['cells']} cells, confident share {v['confident_share']:.3f}")
    for gname, v in (qc.get("units") or {}).items():
        if qc.get("modality") == "bulk":   # round 2 (the second reviewer's lesser point 6): status and reason only
            print(f"  {gname}: {v['status']}" + (f"; {v['reason']}" if v.get("reason") else ""))
            continue
        passed = sorted(v.get("groups", {}))
        print(f"  {gname}: {v['status']}; groups with a factor: {passed}")
    return 0 if qc["status"] != "failed" else 1


if __name__ == "__main__":
    sys.exit(main())
