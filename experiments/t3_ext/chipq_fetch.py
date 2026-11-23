"""CHIPQ-S, the database part's fetch (T3 candidate of the outside-data round of 2026-11-23; register row X018).

Standard library only. Every file lives under one data root outside every checkout, by default ~/v/vec-ext/t3/chipatlas
on whichever machine runs it (DATA_ROOT; the functions take the root as an argument, which is how the tests reach a
temporary one, and no command-line option takes another path). Requests go to two hosts only, dbarchive.biosciencedbc.jp
and chip-atlas.dbcls.jp, over https, with a generic client string that names no team, project, person or machine.
Genome mm10 only; of the per-antigen files only the target-gene files within 5 kb of the TSS (<antigen>.5.tsv), and
only through the URL the analysis list gives for that antigen. A signal, bigWig, bigBed or BED file is never fetched.

Steps (python -I -B experiments/t3_ext/chipq_fetch.py <step>):
  lists    the analysis list, the experiment list and the archive's README (it carries the licence; saved and hashed,
           nothing of it printed), into <root>/lists/.
  rows     for each query antigen held for mm10 (Ctnnb1, Gata4; Mab21l2 only if held) its experiment-list rows go to
           <root>/screen/<antigen>.rows.tsv for the two screening readers, and only the number of rows is printed.
           For every other mm10 antigen with a target-gene file the committed keyword filter (DROP_PATTERNS) runs on
           its experiment rows, and <root>/screen/placebo_kept.json maps each such antigen to its kept experiment ids,
           sorted. Only the number of antigens that keep at least two experiments is printed.
  targets  the target-gene file of each query antigen and of each placebo antigen with at least two kept experiments,
           into <root>/targets/. Prints counts and total bytes only.
  verify   every file of the manifest re-hashed against it; exit code only (0 all match, 2 not).

The manifest (<root>/MANIFEST.json) records for every file its URL (None for a file the rows step derives), size,
sha256 and fetch time. Nothing is deleted or overwritten: a file already on disk is kept when it matches its manifest
entry and refused otherwise (also a file on disk that the manifest does not know, and a leftover .part). A file whose
bytes differ from its manifest entry is refused by every step that reads it. Only the manifest itself is rewritten
(through a .part file and a rename) when a step adds entries; an entry once written is never changed.

Nothing here opens a target-gene file's content: targets downloads and hashes bytes. The only reader of such a file is
chipq.open_target. Exit codes: 0 done, 2 refused.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

DATA_ROOT = Path("~/v/vec-ext/t3/chipatlas")      # outside every checkout, on whichever machine runs it (X018)
HOSTS = ("dbarchive.biosciencedbc.jp", "chip-atlas.dbcls.jp")
USER_AGENT = "python-urllib (public metadata and data download)"   # generic, as experiments/t3_ext/fetch.py
GENOME = "mm10"
ARCHIVE = "https://dbarchive.biosciencedbc.jp/data/chip-atlas/LATEST/"
LISTS = {  # name in <root>/lists/ -> URL (the DBCLS archive's download page)
    "chip_atlas_analysis_list.zip": ARCHIVE + "chip_atlas_analysis_list.zip",
    "chip_atlas_experiment_list.zip": ARCHIVE + "chip_atlas_experiment_list.zip",
    "README_e.html": ARCHIVE + "README_e.html",   # the archive's README: it states the licence (CC BY-SA 4.0)
}
ANALYSIS = "chip_atlas_analysis_list.zip"
EXPERIMENTS = "chip_atlas_experiment_list.zip"
QUERY_ANTIGENS = ("Ctnnb1", "Gata4", "Mab21l2")   # their rows go to the two readers; never through the keyword filter
PLACEBO_MIN_KEPT = 2                              # the specification: "at least 2 screened-clean experiments"
MAX_BYTES = {"list": 512 * 1024 ** 2, "target": 256 * 1024 ** 2}   # a download past this is refused (.part kept)

# the analysis list's +/-5 kb column: its header cell and the one URL form taken from it (the archive serves the
# same path over https). The antigen named in the URL must be the row's antigen.
ANALYSIS_COLUMNS = ("Antigen", "Cell type class in Colocalization", "Recorded (+) or not (-) in Target Genes",
                    "Genome assembly")
TARGET5_COL = 5
TARGET5_HEADER = re.compile(r"Target Genes \(TSS .{1,3}5k\)")
TARGET5_URL = re.compile(r"https?://dbarchive\.biosciencedbc\.jp/kyushu-u/mm10/target/([A-Za-z0-9._+-]+)\.5\.tsv")
NAME = re.compile(r"[A-Za-z0-9._+-]+")
EXP_ID = re.compile(r"[SDE]RX\d+")
DENY = re.compile(r"bigwig|\.bw$|\.bed(\.gz)?$|bigbed|\.bb$|/eachdata/|/assembled/", re.IGNORECASE)

# ---------------------------------------------------------------------------------------------------------------------
# The keyword filter for placebo antigens (the specification's SCREEN paragraph: "Placebo antigens never reach a reader:
# they pass the same drops as a committed keyword filter in code"). One pattern per word of the drop classes, quoted
# from the specification; each pattern also takes that word's spelling, hyphenation and inflection variants and its
# common abbreviation, and nothing beyond. A row is dropped when ANY pattern matches ANY metadata field of the row (a
# doubt drops). The plain words "untreated" and "wild type" match no pattern. Case-insensitive, except the
# sh/si-prefix form, which is case-sensitive inside its own group so that a gene symbol in capitals is not read as one.

_B = r"(?<![a-z0-9])"   # no letter or digit before
_E = r"(?![a-z])"       # no letter after
DROP_CLASSES = {
    "i": "any mutant, flox, Cre, knockout, knockdown, sh/si/shRNA, CRISPR, transgene, reporter or overexpression",
    "ii": ("any stabilised, activated, exon-3-deleted or constitutively active allele, and any sample described as "
           "Wnt-stimulated, Wnt3a/Wnt-ligand treated, GSK3-inhibitor (CHIR/LiCl) treated, or cultured in a medium "
           "naming a pathway drug"),
    "iii": "any drug, inhibitor, agonist, ligand or growth-factor treatment",
    "iv": "for Gata4, the same list",
}
DROP_PATTERNS: tuple[tuple[str, str, str], ...] = (
    # (class, the quoted word, pattern)
    ("i", "mutant", r"mutant|mutat(?:ed|ion)"),
    ("i", "flox", r"flox|" + _B + r"fl\s*/\s*fl" + _E),
    ("i", "Cre", _B + r"cre(?:ert2?|er)?(?![a-z])"),
    ("i", "knockout", r"knock[\s_-]*out|" + _B + r"ko" + _E),
    ("i", "knockdown", r"knock[\s_-]*down|" + _B + r"kd" + _E),
    ("i", "sh/si/shRNA", _B + r"(?:sh|si)[\s_-]?rna|short[\s_-]*hairpin|small[\s_-]*interfering|"
                         r"(?-i:(?<![A-Za-z0-9])(?:sh|si)(?=[A-Z0-9_-]))"),
    ("i", "CRISPR", r"crispr"),
    ("i", "transgene", r"transgen"),
    ("i", "reporter", r"reporter"),
    ("i", "overexpression", r"over[\s_-]*express|" + _B + r"oe" + _E),
    ("i", "mutant (genotype notation)", r"[+-]\s*/\s*-"),
    ("ii", "stabilised", r"stabili[sz]"),
    ("ii", "activated", r"activated"),
    ("ii", "exon-3-deleted", r"exon[\s_-]*3|delet"),
    ("ii", "constitutively active", r"constitutive"),
    ("ii", "Wnt-stimulated / Wnt3a / Wnt-ligand", r"wnt"),
    ("ii", "stimulated", r"stimulat"),
    ("ii", "GSK3-inhibitor", r"gsk[\s_-]*3"),
    ("ii", "CHIR", _B + r"chir"),
    ("ii", "LiCl", _B + r"licl" + _E + r"|lithium"),
    ("ii", "a medium naming a pathway drug (2i)", _B + r"2i" + _E),
    ("iii", "drug", r"drug"),
    ("iii", "inhibitor", r"inhibit"),
    ("iii", "agonist", r"agonist"),
    ("iii", "ligand", r"ligand"),
    ("iii", "growth-factor", r"growth[\s_-]*factor"),
    ("iii", "treatment", r"(?<!un)(?<!non)(?<!non-)(?<!non )treat"),
)
_DROP = tuple((c, w, re.compile(p, re.IGNORECASE)) for c, w, p in DROP_PATTERNS)


def dropped_by(fields) -> list[str]:
    """The drop classes (with the word) whose pattern matches any field; empty means the row is kept."""
    hits = []
    for c, w, rx in _DROP:
        if any(rx.search(str(f)) for f in fields):
            hits.append(f"{c}:{w}")
    return hits


def keep_row(fields) -> bool:
    return not dropped_by(fields)


# ---------------------------------------------------------------------------------------------------------------------
# Helpers


class Refused(Exception):
    """A rule of this step refuses the input; the message says which (never a gene symbol, a score or a column)."""


def sha256_file(path, chunk: int = 1 << 23) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


REFERENCE_CHECKOUT = Path("~/v/virtual-embryo-challenge")


def outside_repo(path) -> Path:
    """The resolved path, refused inside this checkout, the main checkout or the DGX's reference checkout."""
    p = Path(path).expanduser().resolve()
    roots = {REPO.resolve(), REFERENCE_CHECKOUT.expanduser().resolve()}
    if ".claude" in REPO.resolve().parts:
        parts = REPO.resolve().parts
        roots.add(Path(*parts[:parts.index(".claude")]))
    for r in roots:
        if p == r or r in p.parents:
            raise Refused(f"{p} lies inside a checkout ({r}): the data root is outside every checkout")
    return p


def root_of(root=None) -> Path:
    return outside_repo(DATA_ROOT if root is None else root)


def https_url(url: str) -> str:
    """The URL on an allowed host over https (an http URL of an allowed host is taken over https), no port, user,
    query or fragment, and never a denied file type."""
    p = urllib.parse.urlparse(url.strip())
    if p.scheme not in ("http", "https") or p.netloc not in HOSTS or p.params or p.query or p.fragment:
        raise Refused(f"{url}: not an allowed URL (hosts {HOSTS}, no port, user, query or fragment)")
    if DENY.search(p.path):
        raise Refused(f"{url}: a signal, bigWig or BED path is never fetched")
    return urllib.parse.urlunparse(("https", p.netloc, p.path, "", "", ""))


def _download(url: str, dst: Path, max_bytes: int) -> tuple[int, str]:
    """Stream the URL into dst via dst.part; refuses a URL off the allowed hosts and any existing dst or .part."""
    url = https_url(url)
    part = dst.with_name(dst.name + ".part")
    if dst.exists() or dst.is_symlink() or part.exists() or part.is_symlink():
        raise Refused(f"{dst.name}: already on disk (or a leftover .part): nothing is overwritten; remove it by hand")
    dst.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    h = hashlib.sha256()
    n = 0
    with urllib.request.urlopen(req, timeout=300) as r, open(part, "xb") as f:
        while block := r.read(1 << 20):
            n += len(block)
            if n > max_bytes:
                raise Refused(f"{dst.name}: over the byte limit ({max_bytes}); the .part is left for removal by hand")
            h.update(block)
            f.write(block)
    if dst.exists():
        raise Refused(f"{dst.name}: appeared while downloading: refused")
    part.rename(dst)
    return n, h.hexdigest()


def _rel(root: Path, p: Path) -> str:
    return p.relative_to(root).as_posix()


def load_manifest(root: Path) -> dict:
    mf = root / "MANIFEST.json"
    if mf.is_symlink():
        raise Refused("MANIFEST.json is a symlink: refused")
    if not mf.exists():
        return {"root": "data root of the CHIPQ-S fetch (X018)", "files": {}}
    m = json.loads(mf.read_bytes())
    if not isinstance(m, dict) or not isinstance(m.get("files"), dict):
        raise Refused("MANIFEST.json is not of the fixed form: refused")
    return m


def save_manifest(root: Path, man: dict) -> None:
    mf = root / "MANIFEST.json"
    part = mf.with_name(mf.name + ".part")
    if part.exists() or part.is_symlink():
        raise Refused("MANIFEST.json.part is on disk: remove it by hand")
    root.mkdir(parents=True, exist_ok=True)
    with open(part, "x") as f:
        f.write(json.dumps(man, indent=1, sort_keys=True) + "\n")
    part.replace(mf)


def checked_bytes(root: Path, rel: str, man: dict | None = None) -> bytes:
    """The file's bytes, read once, refused unless they match its manifest entry (size and sha256)."""
    man = load_manifest(root) if man is None else man
    rec = man["files"].get(rel)
    p = root / rel
    if rec is None:
        raise Refused(f"{rel}: not in the manifest: refused")
    if p.is_symlink() or not p.is_file():
        raise Refused(f"{rel}: missing or a symlink: refused")
    raw = p.read_bytes()
    if len(raw) != rec["size"] or sha256_bytes(raw) != rec["sha256"]:
        raise Refused(f"{rel}: its bytes differ from the manifest entry: refused")
    return raw


def _fetch_one(root: Path, man: dict, rel: str, url: str, kind: str) -> int:
    """Fetch url to root/rel once. An entry already in the manifest is kept when the file on disk matches it (no
    request) and refused otherwise; a file on disk the manifest does not know is refused."""
    p = root / rel
    rec = man["files"].get(rel)
    if rec is not None:
        if rec.get("url") != url:
            raise Refused(f"{rel}: the manifest records another URL: refused")
        checked_bytes(root, rel, man)
        return rec["size"]
    size, digest = _download(url, p, MAX_BYTES[kind])
    man["files"][rel] = {"url": url, "size": size, "sha256": digest,
                         "fetched": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())}
    save_manifest(root, man)
    return size


def _write_derived(root: Path, man: dict, rel: str, raw: bytes, source: list[str]) -> None:
    """A file the rows step derives: written once; a rerun that would give other bytes is refused."""
    p = root / rel
    rec = man["files"].get(rel)
    if rec is not None or p.exists() or p.is_symlink():
        if rec is None:
            raise Refused(f"{rel}: on disk but not in the manifest: refused")
        if checked_bytes(root, rel, man) != raw:
            raise Refused(f"{rel}: a rerun gives other bytes than the manifest's file: nothing is overwritten")
        return
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "xb") as f:
        f.write(raw)
    man["files"][rel] = {"url": None, "derived_from": source, "size": len(raw), "sha256": sha256_bytes(raw),
                         "fetched": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())}
    save_manifest(root, man)


# ---------------------------------------------------------------------------------------------------------------------
# The two lists


def _zip_member(raw: bytes, what: str) -> bytes:
    try:
        z = zipfile.ZipFile(io.BytesIO(raw))
    except zipfile.BadZipFile:
        raise Refused(f"{what}: not a zip archive: refused") from None
    names = [n for n in z.namelist() if not n.endswith("/")]
    if len(names) != 1:
        raise Refused(f"{what}: {len(names)} members, not one: refused")
    return z.read(names[0])


def _text(raw: bytes) -> str:
    # the archive's lists are ASCII but for a few bytes in a header (the +/- sign in a legacy encoding) and some
    # authors' metadata; undecodable bytes become U+FFFD, which no antigen name or experiment id can contain
    return raw.decode("utf-8", errors="replace")


def analysis_rows(raw_zip: bytes) -> dict[str, str]:
    """{antigen: its +/-5 kb target-file URL} for the mm10 rows of the analysis list that record target genes."""
    rows = list(csv.reader(io.StringIO(_text(_zip_member(raw_zip, ANALYSIS)), newline="")))
    if not rows:
        raise Refused("the analysis list is empty: refused")
    head = rows[0]
    if tuple(head[:4]) != ANALYSIS_COLUMNS or len(head) <= TARGET5_COL or not TARGET5_HEADER.fullmatch(
            head[TARGET5_COL]):
        raise Refused("the analysis list's header is not of the fixed form: refused")
    out: dict[str, str] = {}
    for r in rows[1:]:
        if len(r) != len(head):
            raise Refused("an analysis-list row has another number of fields than the header: refused")
        if r[3] != GENOME or r[2] != "+":
            continue
        m = TARGET5_URL.fullmatch(r[TARGET5_COL])
        if not NAME.fullmatch(r[0]) or m is None or m.group(1) != r[0]:
            raise Refused("an mm10 analysis-list row has a name or +/-5 kb URL outside the fixed form: refused")
        if r[0] in out:
            raise Refused("an mm10 antigen appears twice in the analysis list: refused")
        out[r[0]] = r[TARGET5_COL]
    return out


# The experiment list's layout (README_e.html, "Experiment list", and the file itself: one CSV member with this header
# and 17 fields per row). The keyword filter reads the metadata fields SCREENED_COLS (antigen class, antigen, cell-type
# class, cell type, cell-type description, title, the authors' metadata): not the id, the genome, the processing logs
# (read counts) or the seven file links.
EXPERIMENT_COLUMNS = ("Experimental ID", "Genome assembly", "Antigen class", "Antigen", "Cell type class", "Cell type",
                      "Cell type description", "Processing logs", "Title", "Meta data", "BigWig",
                      "Peak-call (BED) (q < 1E-05)", "Peak-call (BED) (q < 1E-10)", "Peak-call (BED) (q < 1E-20)",
                      "Peak-call (BigBed) (q < 1E-05)", "Peak-call (BigBed) (q < 1E-10)",
                      "Peak-call (BigBed) (q < 1E-20)")
EXP_ID_COL, EXP_GENOME_COL, EXP_ANTIGEN_COL = 0, 1, 3
SCREENED_COLS = (2, 3, 4, 5, 6, 8, 9)


def experiment_rows(raw_zip: bytes):
    """The data rows of the experiment list, each a list of its 17 fields; refuses any other form."""
    r = csv.reader(io.StringIO(_text(_zip_member(raw_zip, EXPERIMENTS)), newline=""))
    head = next(r, None)
    if head is None or tuple(head) != EXPERIMENT_COLUMNS:
        raise Refused("the experiment list's header is not of the fixed form: refused")
    for f in r:
        if len(f) != len(EXPERIMENT_COLUMNS) or not EXP_ID.fullmatch(f[EXP_ID_COL]):
            raise Refused("an experiment-list row is not of the fixed form: refused")
        yield f


def screened_fields(f: list[str]) -> list[str]:
    return [f[i] for i in SCREENED_COLS]


def rows_tsv(rows: list[list[str]]) -> bytes:
    """The readers' file: the header and the rows, every field as the list has it, tab-separated (a tab or line
    break inside a field becomes one space)."""
    def clean(x: str) -> str:
        return re.sub(r"[\t\r\n]", " ", x)

    lines = ["\t".join(EXPERIMENT_COLUMNS)] + ["\t".join(clean(x) for x in f) for f in rows]
    return ("\n".join(lines) + "\n").encode("utf-8")


# ---------------------------------------------------------------------------------------------------------------------
# Steps


def step_lists(root=None) -> dict:
    root = root_of(root)
    man = load_manifest(root)
    sizes = {n: _fetch_one(root, man, f"lists/{n}", url, "list") for n, url in LISTS.items()}
    return {"files": len(sizes)}


def _lists(root: Path, man: dict) -> tuple[dict[str, str], bytes]:
    targets = analysis_rows(checked_bytes(root, f"lists/{ANALYSIS}", man))
    return targets, checked_bytes(root, f"lists/{EXPERIMENTS}", man)


def step_rows(root=None) -> dict:
    root = root_of(root)
    man = load_manifest(root)
    targets, exp_raw = _lists(root, man)
    src = [f"lists/{ANALYSIS}", f"lists/{EXPERIMENTS}"]
    query_rows: dict[str, list[list[str]]] = {q: [] for q in QUERY_ANTIGENS}
    placebo: dict[str, list[str]] = {a: [] for a in targets if a not in QUERY_ANTIGENS}
    seen: set[str] = set()
    for f in experiment_rows(exp_raw):
        if f[EXP_GENOME_COL] != GENOME:
            continue
        if f[EXP_ID_COL] in seen:
            raise Refused("an mm10 experiment id appears twice in the experiment list: refused")
        seen.add(f[EXP_ID_COL])
        a = f[EXP_ANTIGEN_COL]
        if a in query_rows:
            query_rows[a].append(f)
        elif a in placebo and keep_row(screened_fields(f)):
            placebo[a].append(f[EXP_ID_COL])
    counts = {}
    for q in QUERY_ANTIGENS:
        held = q in targets
        counts[q] = len(query_rows[q]) if held else 0
        if held and query_rows[q]:
            _write_derived(root, man, f"screen/{q}.rows.tsv", rows_tsv(query_rows[q]), src)
    kept = {a: sorted(v) for a, v in sorted(placebo.items())}
    raw = (json.dumps({"filter": [list(p) for p in DROP_PATTERNS], "kept": kept}, indent=1, sort_keys=True)
           + "\n").encode("utf-8")
    _write_derived(root, man, "screen/placebo_kept.json", raw, src)
    return {"rows": counts, "placebo_at_least_two": sum(1 for v in kept.values() if len(v) >= PLACEBO_MIN_KEPT)}


def placebo_kept(root: Path, man: dict | None = None) -> dict[str, list[str]]:
    """The keyword filter's kept ids per placebo antigen, read from its manifest-checked file."""
    raw = checked_bytes(root, "screen/placebo_kept.json", man)
    d = json.loads(raw)
    if [list(p) for p in DROP_PATTERNS] != d.get("filter"):
        raise Refused("placebo_kept.json was written under another keyword filter: refused")
    return d["kept"]


def target_rel(antigen: str) -> str:
    if not NAME.fullmatch(antigen):
        raise Refused("an antigen name outside the fixed form: refused")
    return f"targets/{antigen}.5.tsv"


def step_targets(root=None) -> dict:
    root = root_of(root)
    man = load_manifest(root)
    targets, _ = _lists(root, man)
    want = [q for q in QUERY_ANTIGENS if q in targets]
    kept = placebo_kept(root, man)
    want += sorted(a for a, v in kept.items() if len(v) >= PLACEBO_MIN_KEPT)
    total = 0
    for a in want:
        if a not in targets:
            raise Refused("a placebo antigen has no target-gene file in the analysis list: refused")
        total += _fetch_one(root, man, target_rel(a), targets[a], "target")
    return {"query_files": sum(1 for a in want if a in QUERY_ANTIGENS),
            "placebo_files": sum(1 for a in want if a not in QUERY_ANTIGENS), "total_bytes": total}


def step_verify(root=None) -> bool:
    root = root_of(root)
    man = load_manifest(root)
    try:
        for rel in sorted(man["files"]):
            checked_bytes(root, rel, man)
    except Refused:
        return False
    return bool(man["files"])


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("step", choices=("lists", "rows", "targets", "verify"))
    a = ap.parse_args(argv)
    try:
        if a.step == "lists":
            r = step_lists()
            print(f"lists: {r['files']} files in the manifest")
        elif a.step == "rows":
            r = step_rows()
            for q, n in r["rows"].items():
                print(f"{q}: {n} rows")
            print(f"placebo antigens keeping at least {PLACEBO_MIN_KEPT} experiments: {r['placebo_at_least_two']}")
        elif a.step == "targets":
            r = step_targets()
            print(f"query files {r['query_files']}, placebo files {r['placebo_files']}, total bytes {r['total_bytes']}")
        else:
            return 0 if step_verify() else 2
    except Refused as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
