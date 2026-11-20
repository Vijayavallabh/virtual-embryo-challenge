#!/usr/bin/env python
"""Write experiments/t3_ext/panel_ids.tsv, the gene-id map of the bulk Ensembl route (stage 9; README "Gene ids").

No annotation file is downloaded (the main session's decision on the design's identifier point). The source is the
features file of an ALREADY PINNED AND FETCHED 10x sample: raw/<ACC>/<GSM>/<FILE> under the fetch destination, a pinned
file of a pinned sample of a cleared dataset whose manifest record holds the current pins and the file's sha256. FILE
is a loose features (or genes) file, or a per-sample archive whose one features member is read in memory after
process.archive_members has checked every member. Only gene-expression rows are read: Ensembl gene id, symbol.

The file holds one row per reference gene, in the reference's gene order (refstats.npz): symbol, ensembl_gene_id
(NA when the source gives the symbol no id or several ids, or the id names two panel symbols), chromosome (NA: a 10x
features file states none), reason. Header lines record the source file's name and sha256, the reference gene list's
sha256, and one '# extra' line per sex gene of the sex call (t3x_common.SEX_GENES) and per symbol named with --extra
(the knocked-out genes, whose own rows G9 reads). The script prints the new file's sha256; the main session puts it
into t3x_common.PANEL_IDS_SHA256 in a reviewed commit together with the file. Until then every dataset that needs the
Ensembl route is refused.

  uv run --frozen --no-sync python -I -B experiments/t3_ext/make_panel_ids.py --dest ~/v/vec-ext/t3 \\
      --accession GSE168768 --gsm GSM5169142 --file GSM5169142_features.tsv.gz \\
      --refstats ~/v/vec-ext/t3/reference/refstats.npz --extra Sox7 --extra Zmym2 [--out PATH]
"""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import t3x_common as cm  # noqa: E402
import t3x_tables as tb  # noqa: E402

Refused = cm.Refused
COLUMNS = ("symbol", "ensembl_gene_id", "chromosome", "reason")


def source_file(dest: Path, acc: str, gsm: str, name: str, datasets: Path | None = None) -> Path:
    """The pinned, fetched file: refused unless the dataset is cleared and pins the sample with this file, and the
    manifest record of raw/<ACC>/<GSM>/<FILE> holds the current pins and the file's sha256."""
    ds = cm.dataset(acc, datasets)
    if cm.is_bulk(ds):
        raise Refused(f"{acc}: the source is a single-cell (10x) sample")
    if not any(p["accession"] == gsm and name in p["files"] for p in ds.get("pinned") or []):
        raise Refused(f"{acc}: {gsm} is not pinned with the file {name}")
    rdir = cm.outside_repo(dest) / "raw" / acc
    man = json.loads((rdir / "MANIFEST.json").read_text()) if (rdir / "MANIFEST.json").is_file() else {}
    rec = (man.get("files") or {}).get(f"{gsm}/{name}")
    p = rdir / gsm / name
    if rec is None or p.is_symlink() or not p.is_file() or rec.get("pins_sha256") != cm.pins_sha256(ds) \
            or cm.sha256_file(p) != rec.get("sha256"):
        raise Refused(f"{acc}/{gsm}/{name}: not a file fetched under the current pins (manifest record and sha256)")
    return p


def features_rows(p: Path) -> tuple[list[tuple[str, str]], str]:
    """(gene-expression (id, symbol) rows, the member read or ''), from a loose features file or a per-sample archive
    (one features member; every member checked by process.archive_members)."""
    if p.name.endswith((".tar.gz", ".tgz", ".tar")):
        import process as pr

        with tarfile.open(p, "r:*") as t:
            items = [it for it in pr.archive_members(t, p.name, raw_ok=True) if it[1] == "features"]
            if len(items) != 1:
                raise Refused(f"{p.name}: {len(items)} features members: one is needed")
            f = t.extractfile(items[0][2])
            raw = f.read() if f is not None else b""
            member = "/".join(items[0][0])
    else:
        if cm.triple_kind(p.name) != "features":
            raise Refused(f"{p.name}: not a features or genes file")
        raw, member = p.read_bytes(), ""
    if (member or p.name).endswith(".gz"):
        raw = gzip.decompress(raw)
    rows = []
    for ln in raw.decode("utf-8").splitlines():
        f = ln.split("\t")
        if len(f) < 2 or not ln.strip():
            continue
        typ = f[2] if len(f) > 2 else "Gene Expression"
        if "gene" in typ.lower() and "expression" in typ.lower():
            rows.append((f[0], f[1]))
    if not rows:
        raise Refused(f"{p.name}: no gene-expression row")
    return rows, member


def build(rows: list[tuple[str, str]], genes: list[str], extra: list[str], source: str, source_sha: str) -> str:
    """The text of panel_ids.tsv (see the module docstring)."""
    ids_of: dict[str, list[str]] = {}
    for i, s in rows:
        ids_of.setdefault(s, [])
        if i not in ids_of[s]:
            ids_of[s].append(i)
    one = {s: v[0] for s, v in ids_of.items() if len(v) == 1}
    owners: dict[str, set] = {}
    for g in genes:
        if g in one:
            owners.setdefault(one[g], set()).add(g)
    out = ["# panel_ids.tsv: the gene-id map of the bulk Ensembl route (experiments/t3_ext/make_panel_ids.py); one row "
           "per reference gene, in the reference's order",
           f"# source_file\t{source}", f"# source_sha256\t{source_sha}",
           f"# reference_genes\t{len(genes)}",
           f"# reference_genes_sha256\t{hashlib.sha256(chr(10).join(genes).encode('utf-8')).hexdigest()}",
           "# chromosome\tNA for every gene: a 10x features file states no chromosome"]
    for s in [*cm.SEX_GENES, *[e for e in extra if e not in cm.SEX_GENES]]:
        v = ids_of.get(s, [])
        out.append(f"# extra\t{s}\t{v[0] if len(v) == 1 else 'NA'}\t"
                   + ("one id in the source" if len(v) == 1 else f"{len(v)} ids in the source"))
    out.append("\t".join(COLUMNS))
    for g in genes:
        v = ids_of.get(g, [])
        if len(v) == 1 and len(owners[v[0]]) > 1:
            out.append(f"{g}\tNA\tNA\tthe id names two panel symbols in the source")
        elif len(v) == 1:
            out.append(f"{g}\t{v[0]}\tNA\tone id in the source")
        else:
            out.append(f"{g}\tNA\tNA\t{'no id' if not v else f'{len(v)} ids'} in the source")
    return "\n".join(out) + "\n"


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dest", type=Path, required=True)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--gsm", required=True)
    ap.add_argument("--file", required=True)
    ap.add_argument("--refstats", type=Path, required=True)
    ap.add_argument("--extra", action="append", default=[])
    ap.add_argument("--extra-file", type=Path)
    ap.add_argument("--out", type=Path, default=cm.PANEL_IDS)
    a = ap.parse_args(argv)
    try:
        p = source_file(a.dest, a.accession, a.gsm, a.file)
        rows, member = features_rows(p)
        extra = list(a.extra)
        if a.extra_file:
            extra += [ln.strip() for ln in a.extra_file.read_text().splitlines() if ln.strip()]
        genes = tb.load_ref_stats(a.refstats)["genes"]
        src = f"{a.accession}/{a.gsm}/{a.file}" + (f":{member}" if member else "")
        text = build(rows, genes, extra, src, cm.sha256_file(p))
    except Refused as e:
        print(f"refused: {e}")
        return 2
    a.out.write_text(text, encoding="utf-8")
    print(f"wrote {a.out}: {len(genes)} rows; sha256 {cm.sha256_file(a.out)}")
    print("put that value into t3x_common.PANEL_IDS_SHA256 in the same reviewed commit as the file")
    return 0


if __name__ == "__main__":
    sys.exit(main())
