"""CHIPQ-S, the database part's fetch (T3 candidate of the outside-data round of 2026-11-23; register row X018).

Standard library only. Every file lives under one data root outside every checkout, by default ~/v/vec-ext/t3/chipatlas
on whichever machine runs it (DATA_ROOT; the functions take the root as an argument, which is how the tests reach a
temporary one, and no command-line option takes another path). Requests go to two hosts only, dbarchive.biosciencedbc.jp
and chip-atlas.dbcls.jp, over https, with a generic client string that names no team, project, person or machine. No
redirect is followed: a redirect, or an answer from any other URL than the one asked for, is refused before a byte is
written.
Genome mm10 only; of the per-antigen files only the target-gene files within 5 kb of the TSS (<antigen>.5.tsv), and
only through the URL the analysis list gives for that antigen. A signal, bigWig, bigBed or BED file is never fetched.

Steps (python -I -B experiments/t3_ext/chipq_fetch.py <step>):
  lists    the analysis list and the experiment list, into <root>/lists/. With the target-gene files these are the
           specification's three kinds of file and the only ones fetched (the licence words are in the specification
           and go into register row X018 as the archive states them).
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
entry and refused otherwise (also a file on disk that the manifest does not know). A partial download (.part),
left by a failed attempt or an earlier run, is moved aside as <name>.failed.<n>: never deleted, never read, refused
nowhere. A failed download is tried again up to RETRIES times after a short pause, and every request waits
REQUEST_PAUSE first. A file whose
bytes differ from its manifest entry is refused by every step that reads it. Only the manifest itself is rewritten
(through a .part file and a rename) when a step adds entries; an entry once written is never changed.

No message names a target-gene file by its path (the path is its antigen's name): file_label names it by its manifest
entry's sha256 (first 12 hex) instead, and any other error prints its type only.

Nothing here opens a target-gene file's content: targets downloads and hashes bytes. The only reader of such a file is
chipq.open_target. Exit codes: 0 done, 2 refused.
"""

from __future__ import annotations

import csv
import hashlib
import http.client
import io
import json
import re
import sys
import time
import urllib.error
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
}
ANALYSIS = "chip_atlas_analysis_list.zip"
EXPERIMENTS = "chip_atlas_experiment_list.zip"
QUERY_ANTIGENS = ("Ctnnb1", "Gata4", "Mab21l2")   # their rows go to the two readers; never through the keyword filter
PLACEBO_MIN_KEPT = 2                              # the specification: "at least 2 screened-clean experiments"
MAX_BYTES = {"list": 512 * 1024 ** 2, "target": 256 * 1024 ** 2}   # a download past this is refused (.part set aside)
RETRIES = 3            # a failed download (a network or OS error, not a refusal) is tried again up to three times
RETRY_PAUSE = 10.0     # seconds before each retry
REQUEST_PAUSE = 1.0    # seconds before every request, so hundreds of target files are fetched at a gentle pace
_sleep = time.sleep    # (the tests replace it)

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
# doubt drops). The plain words "untreated" and "wild type" match no pattern, nor does a field whose whole value
# states no treatment (below). Case-insensitive, except the groups marked (?-i:...), which are case-sensitive so that a
# gene symbol in capitals is not read as one. The forms marked kwdecide are the rule agent's A13 addendum: spellings,
# notations or abbreviations of the paragraph's own words; no class and no word is added.

_B = r"(?<![a-z0-9])"   # no letter or digit before
_E = r"(?![a-z])"       # no letter after
_L = r"(?<![a-z])"      # no letter before (a digit may come before: a gene symbol usually ends in one, "Gata4OE")
_MINUS = "\\-\u2212\u2013"   # hyphen-minus, the minus sign and the en dash (text copied from papers)
_AA = "ACDEFGHIKLMNPQRSTVWY"   # the one-letter residue codes
# The reverse case of "treat" (rule agent kwdecide, A13 addendum): a field whose entire value states no treatment
# ("treatment: none", "treatment=untreated", "treatment: no treatment", or the whole field "no treatment") states Tier
# U's "untreated" and names no class (iii) treatment. Only that whole value is exempt: the key must open the field or
# follow a "|", ";" or newline, the value must end the field or meet one of those separators, and "vehicle", "DMSO",
# "control", "n/a" or an empty value are not exempt (a doubt drops).
_VALUE_END = r"\s*(?:\Z|[|;\n])"
_KEY_START = r"(?:\A|(?<=[|;\n])|(?<=[|;\n] ))"
_NO_TREATMENT_KEY = (r"(?!" + _KEY_START + r"treatments?\s?[:=]\s?(?:none|untreated|no treatment)" + _VALUE_END + ")")
_NO_TREATMENT_VALUE = ("(?!(?:(?<=^no )|"
                       + "|".join(f"(?<=treatment{k}{a}{sep}{b}no )" for k in ("", "s") for sep in ":="
                                  for a in ("", " ") for b in ("", " "))
                       + r")treatment" + _VALUE_END + ")")
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
    # "mutagenized" shares the stem (kwdecide). A point-mutation name, residue-position-residue ("S33Y", "V600E", a stop
    # as "*" or "X"), is the notation of a mutation, as +/- is of a mutant: case-sensitive, with no letter or digit on
    # either side and two to four position digits. A one-digit position is not read as one: it would drop strain,
    # antigen, receptor and tag names and named agents that spell no mutation (C3H, N2A, H2A.Z, E2A, G9A, S5P, D1R,
    # P2A, G2M, S1P, I3C). Look-alikes with two or more digits (T47D, S100A, R26R) still drop and only thin the pool.
    ("i", "mutant", r"mutant|mutat(?:ed|ion)|mutagen|(?-i:(?<![A-Za-z0-9])[" + _AA + r"][1-9][0-9]{1,3}[" + _AA
     + r"X*](?![A-Za-z0-9]))"),
    # fl/fl, fl/+ and fl/- (floxed over null), and f/f, the short form of fl/fl (kwdecide)
    ("i", "flox", r"flox|" + _L + r"fl\s*/\s*(?:fl|\+|[" + _MINUS + r"])(?![a-z])|" + _L + r"f\s*/\s*f" + _E),
    # Cre joined to a driver name is the standard allele notation ("AlbCre", "Vil1cre", "LysMcre", "iCre"). The
    # case-sensitive groups keep "secreted", symbols that begin with Cre+letter ("Creb1", "Crem") and symbols that
    # begin with a capital and "cre" undropped: "cre" counts in any case after a non-letter ("cRE", "CRe"), "cre" after
    # a digit or a lowercase-capital pair, and "Cre" after any letter, a capital driver name included ("CMVCre",
    # "GFAPCre", "CAGGCreER", "UBCCreERT2").
    ("i", "Cre", r"(?-i:(?<![A-Za-z])(?i:cre)|(?<=[0-9])cre|(?<=[a-z][A-Z])[Cc]re|(?<=[a-z])Cre|(?<=[A-Z])Cre)"
     r"(?:ert2?|er)?" + _E),
    # "knocked out", "knock-outs", "KO" and "KOs", the KO abbreviation with a c/d/i prefix ("cKO", "dKO", "iKO") and KO
    # straight after a symbol ("Gata4cKO", "SmadKO"). The prefixed and plural forms are case-sensitive, so symbols
    # such as "Cko1", "Kos2" or "Pkd1" are kept. Knockdown the same with down and KD.
    ("i", "knockout", r"knock(?:ed|ing)?[\s_-]*outs?|" + _L + r"ko" + _E
     + r"|(?-i:(?<![A-Za-z])(?:[cdi]?KOs|[cdi]KO|[cdi]kos?|[CDI]KOs?)|(?<=[a-z])KOs?)" + _E),
    ("i", "knockdown", r"knock(?:ed|ing)?[\s_-]*downs?|" + _L + r"kd" + _E
     + r"|(?-i:(?<![A-Za-z])(?:[cdi]?KDs|[cdi]KD|[cdi]kds?|[CDI]KDs?)|(?<=[a-z])KDs?)" + _E),
    # esiRNA (endoribonuclease-prepared siRNA) and the lower-case control names (shctrl, sictrl, shscr) are notations
    # of the word (kwdecide); B8's refusal of a broad case-insensitive sh/si prefix stands, so only sh/si followed by a
    # control word is matched in any case.
    ("i", "sh/si/shRNA", _B + r"(?:sh|e?si)[\s_-]?rna|short[\s_-]*hairpin|small[\s_-]*interfering|"
                         r"(?-i:(?<![A-Za-z0-9])(?:sh|si)(?=[A-Z0-9_-]))|"
                         + _B + r"(?:sh|si)[\s_-]?(?:ctrl|ctl|control|scr(?:ambled?)?)"),
    # Cas9 and dCas9: Cas stands for CRISPR-associated (kwdecide). sgRNA alone is not matched: single guide RNA has no
    # form of the word, and the filter adds no word (A13).
    ("i", "CRISPR", r"crispr|cas[\s_-]?9(?![0-9])"),
    # Tg, the allele-notation abbreviation of transgene/transgenic ("Tg(...)", "Tg/+"), as a whole token: case-sensitive
    # Tg or tg, with no letter or digit before and no letter after (kwdecide)
    ("i", "transgene", r"transgen|(?-i:(?<![A-Za-z0-9])[Tt]g(?![a-z]))"),
    ("i", "reporter", r"reporter"),
    ("i", "overexpression", r"over[\s_-]*express|" + _L + r"oe" + _E),
    # +/-, -/- and the reversed heterozygous -/+ (kwdecide)
    ("i", "mutant (genotype notation)", "[+" + _MINUS + r"]\s*/\s*[" + _MINUS + "]|[" + _MINUS + r"]\s*/\s*\+"),
    ("ii", "stabilised", r"stabili[sz]"),
    ("ii", "activated", r"activated"),
    # ex3, delta-ex3, Δex3, dex3 and del(ex3) are notations of an exon-3 deletion (kwdecide; a bare ex3 drops as a
    # doubt). "dex" takes no space before the 3, so "Dex 3uM" (a named drug) is not read as one.
    ("ii", "exon-3-deleted", r"exon[\s_-]*3|delet|(?<![a-z0-9])(?:ex[\s_-]?3|d(?:elta)?[_-]?ex3)(?![0-9])"),
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
    ("iii", "treatment", r"(?<!un)(?<!non)(?<!non-)(?<!non )" + _NO_TREATMENT_KEY + _NO_TREATMENT_VALUE + r"treat"),
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


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    """Every redirect is refused: no request leaves the URL it was given, so none can leave the two hosts."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise Refused("the server answered with a redirect: refused (no redirect is followed)")


def _opener(*extra):
    """The one opener every request goes through: urllib's default handlers with the redirect handler replaced (the
    tests add a fake https handler through `extra`)."""
    return urllib.request.build_opener(_NoRedirect, *extra)


def file_label(rel: str, rec: dict | None = None, path: Path | None = None) -> str:
    """How a message names a file. A target-gene file's path is its antigen's name, and no message may carry a placebo
    antigen's name, so such a file is named by a non-identifying key: its manifest entry's sha256 (first 12 hex), or,
    for a file the manifest does not know, the sha256 of its bytes on disk, or nothing. Other files keep their path."""
    if not rel.startswith("targets/"):
        return rel
    if rec is not None and isinstance(rec.get("sha256"), str):
        return f"a target-gene file (manifest entry sha256 {rec['sha256'][:12]})"
    if path is not None and path.is_file() and not path.is_symlink():
        return f"a target-gene file not in the manifest (its bytes' sha256 {sha256_file(path)[:12]})"
    return "a target-gene file not in the manifest"


def _download(url: str, dst: Path, max_bytes: int, label: str | None = None) -> tuple[int, str]:
    """Stream the URL into dst via dst.part; refuses a URL off the allowed hosts, any redirect (no .part is made) and
    any existing dst or .part. Messages name the file by `label` (default: its name), never by its URL."""
    label = dst.name if label is None else label
    try:
        url = https_url(url)
    except Refused:
        raise Refused(f"{label}: its URL is not an allowed URL: refused") from None
    part = dst.with_name(dst.name + ".part")
    if dst.exists() or dst.is_symlink() or part.exists() or part.is_symlink():
        raise Refused(f"{label}: already on disk but not in the manifest (or a .part beside it): nothing is "
                      "overwritten: refused")
    dst.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    h = hashlib.sha256()
    n = 0
    with _opener().open(req, timeout=300) as r:
        if r.geturl() != url:
            raise Refused("the answer came from another URL than the one asked for: refused")
        with open(part, "xb") as f:
            while block := r.read(1 << 20):
                n += len(block)
                if n > max_bytes:
                    raise Refused(f"{label}: over the byte limit ({max_bytes}): refused")
                h.update(block)
                f.write(block)
    if dst.exists():
        raise Refused(f"{label}: appeared while downloading: refused")
    part.rename(dst)
    return n, h.hexdigest()


def set_aside(part: Path) -> Path | None:
    """Move a partial download aside under <name>.failed.<n> (the first n not on disk). Such a file is never deleted,
    never read and refused nowhere; a rerun then needs no removal by hand. None when there is no .part."""
    if not (part.exists() or part.is_symlink()):
        return None
    base = part.name[:-len(".part")] if part.name.endswith(".part") else part.name
    n = 1
    while (part.with_name(f"{base}.failed.{n}")).exists() or part.with_name(f"{base}.failed.{n}").is_symlink():
        n += 1
    dst = part.with_name(f"{base}.failed.{n}")
    part.rename(dst)
    return dst


_TRANSIENT = (OSError, urllib.error.URLError, http.client.HTTPException)


def _download_retried(url: str, dst: Path, max_bytes: int, label: str) -> tuple[int, str]:
    """_download after a short pause, tried again up to RETRIES times on a network or OS error (never on a refusal).
    A leftover .part from an earlier run, and the .part of every failed or refused attempt, is set aside first."""
    part = dst.with_name(dst.name + ".part")
    set_aside(part)
    last = ""
    for attempt in range(1 + RETRIES):
        _sleep(REQUEST_PAUSE if attempt == 0 else RETRY_PAUSE)
        try:
            return _download(url, dst, max_bytes, label)
        except Refused:
            set_aside(part)
            raise
        except _TRANSIENT as e:
            set_aside(part)
            last = type(e).__name__
    raise Refused(f"{label}: the download failed {1 + RETRIES} times (last: {last}): refused; partial files were "
                  "set aside as .failed.<n>")


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
        raise Refused(f"{file_label(rel)}: not in the manifest: refused")
    if p.is_symlink() or not p.is_file():
        raise Refused(f"{file_label(rel, rec)}: missing or a symlink: refused")
    raw = p.read_bytes()
    if len(raw) != rec["size"] or sha256_bytes(raw) != rec["sha256"]:
        raise Refused(f"{file_label(rel, rec)}: its bytes differ from the manifest entry: refused")
    return raw


def _fetch_one(root: Path, man: dict, rel: str, url: str, kind: str) -> int:
    """Fetch url to root/rel once. An entry already in the manifest is kept when the file on disk matches it (no
    request) and refused otherwise; a file on disk the manifest does not know is refused; a leftover .part is set
    aside (set_aside) and the file fetched again."""
    p = root / rel
    rec = man["files"].get(rel)
    if rec is not None:
        if rec.get("url") != url:
            raise Refused(f"{file_label(rel, rec)}: the manifest records another URL: refused")
        checked_bytes(root, rel, man)
        return rec["size"]
    size, digest = _download_retried(url, p, MAX_BYTES[kind], file_label(rel, None, p))
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


# The experiment list's layout (the archive's documentation, "Experiment list", and the file itself: one CSV member
# with this header and 17 fields per row). The keyword filter reads the metadata fields SCREENED_COLS (antigen class,
# antigen, cell-type class, cell type, cell-type description, title, the authors' metadata): not the id, the genome,
# the processing logs (read counts) or the seven file links.
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
    except Exception as e:  # noqa: BLE001  (an OS or network message could carry a file's path or URL: type only)
        print(f"REFUSED: {type(e).__name__} during the step", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
