"""Shared pieces of the outside-knockout factor tables (experiments/t3_ext): the fixed constants, the dataset list's
loader and checks, the pinned sample list, the stage cap, sha256, and the guard that keeps every output outside the
repository.

Every constant below was fixed on 2026-11-20, before any outside count was read, from design numbers only (the plan
of the T3 outside-knockout campaign and the challenger's required changes). None comes from a score.

Two pins live here (stage 3, 2026-11-20):
  NEVER_FETCH       the held and excluded accessions, a held dataset's second screened series included (GSE280587
                    of GSE280543, stage 4). dataset() refuses them whatever a YAML says, and load_datasets() refuses
                    a YAML that gives any of them a status other than held_phenocopy or excluded_design. Moving one
                    out needs a reviewed commit that cites the captain's ruling. The held set the code checks is
                    NEVER_FETCH together with the accession and every screened series of every dataset that is not
                    cleared (held_series()), so a held series missing from this constant is still held.
  DATASETS_SHA256   the sha256 of the committed datasets.yaml. Reading the committed file refuses when its bytes no
                    longer hash to this value, so the pin-and-commit step (README, "The pinned sample list") updates
                    both in the same commit. Tests reach other YAML content only through a function argument; no CLI
                    takes another path.

Stage 6 (2026-11-20, campaign/t3-ext-archives) adds two declarations, checked here for form:
  series_files      per named series, the patterns a file of that series' OWN supplementary list must match to be
                    offered; pinned_series_files maps each pinned series file (one archive, or a loose 10x triple) to
                    the pinned samples whose cells it holds, with a split rule (check_series_file_pins).
  control_series    a sister series of the same study that supplies controls only; it is screened, shown, hashed and
                    pinned like the dataset's own series, with its relation and a reader note stated
                    (check_control_series). A knockout unit whose controls come from it takes d = SISTER_D.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DATASETS = HERE / "datasets.yaml"
EXTERNAL = REPO / "experiments" / "external"
if str(EXTERNAL) not in sys.path:
    sys.path.insert(0, str(EXTERNAL))
import windows  # noqa: E402  (experiments/external/windows.py: parse_stage, in_window)

# ---------------------------------------------------------------------------------------------------------------------
# Fixed constants (2026-11-20; no score in any of them)

CONFIDENT = 0.6          # a cell's group call counts only with label_conf >= 0.6 (rules.yaml confident_call)
MIN_CONF_SHARE = 0.5     # a (stage) drops if fewer than half of either arm's cells are called confidently
SHARE_RATIO_MAX = 2.0    # a (group, stage) drops if its share of an arm's cells differs by more than 2x across arms
SPEARMAN_MARGIN = 0.05   # the control group's best-matching reference group must be its own, ahead by 0.05
AUTHOR_AGREE = 0.5       # where author labels exist: >= half of the group's control calls agree with the dictionary
MIN_CELLS = 100          # confident cells of the group in EACH arm at a stage (the frozen transfer's minimum)
MIN_UMI = 20             # a gene with fewer summed UMIs in any arm of (group, stage) is NaN there (not measured)
PSEUDOCOUNT = 10.0       # per-gene count pseudocount on the summed UMIs (the challenger's change 3)
N_HALF = 200             # the cell term n_h / (n_h + 200); 200 is K7's MIN_CELLS
CLIP = (0.25, 4.0)       # the frozen clip of gene_factor
STAGE_MAX = 9.5          # every outside cell is capped at E9.5 inclusive (the T1 window, applied to every cell)
MIN_GENES_PER_CELL = 200  # a cell with fewer detected genes (all features) is not a cell
MAX_BARCODES = 60000     # a per-sample matrix with more barcodes is an unfiltered matrix: refused
MIN_PANEL_MATCH = 0.5    # a dataset matching fewer than half of the panel genes fails its QC
DESIGNS = ("germline", "conditional", "f0_crispr", "chimera")
RELATIONS = ("littermate", "same_capture", "same_series", "heterozygous")
TIERS = ("wnt", "cardiac", "unrelated")
STATUSES = ("cleared", "held_phenocopy", "excluded_design")
# Held and excluded datasets: never fetched by any flag or YAML (the main session's decision of 2026-11-20).
NEVER_FETCH = frozenset({"GSE165300", "GSE126128", "GSE280543", "GSE280587", "GSE259342", "E-MTAB-7325"})
# sha256 of the committed datasets.yaml; the pin-and-commit step updates it in the same commit as the YAML.
DATASETS_SHA256 = "2b18574c6d26fd57eeda84fc29a336024104cce328240d72711d628c6b46cf6e"
PIN_KEYS = ("accession", "arm", "gene", "stage", "files", "metadata_sha256")
SERIES_PIN_KEYS = ("accession", "series_metadata_sha256")
PIN_ARMS = ("ko", "control")
# stage 6: a series-level file of a dataset (README, "Series-level files"), and a sister-series control
SERIES_FILE_UNIT_KEYS = ("series", "files", "samples", "split")
SPLITS = ("none", "barcode_suffix", "folder")
CONTROL_SERIES_KEYS = ("accession", "screen_series", "relation", "reader_note")
RELATION_KEYS = ("same_lab", "same_protocol", "same_publication")
# d for a knockout unit whose controls come from a sister series: the weaker value design_factor gives controls that
# are not littermates or same-capture (fixed 2026-11-20, stage 6; no score)
SISTER_D = 0.5
PLACEHOLDER = "TO FILL"
# Aliases matched only after an exact symbol fails (both directions); nothing else is renamed.
ALIASES = {"T": "Tbxt", "Tbxt": "T", "Hbb-y": "Hbb-bh0", "Hbb-bh0": "Hbb-y"}


class Refused(Exception):
    """A rule of this folder refuses the input; the message says which."""


def sha256_file(path, chunk: int = 1 << 23) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


# the DGX's reference checkout stays exactly as it was rehearsed: nothing of this folder writes into it
REFERENCE_CHECKOUT = Path.home() / "v" / "virtual-embryo-challenge"


def outside_repo(path, repo: Path | None = None) -> Path:
    """The resolved path, refused if it lies inside the repository (this checkout or the main checkout) or inside the
    DGX's reference checkout (~/v/virtual-embryo-challenge)."""
    p = Path(path).expanduser().resolve()
    roots = {(repo or REPO).resolve(), REFERENCE_CHECKOUT.expanduser().resolve()}
    main = (repo or REPO).resolve()
    if ".claude" in main.parts:  # a worktree under <main>/.claude/worktrees/<name>
        roots.add(Path(*main.parts[:main.parts.index(".claude")]))
    for r in roots:
        if p == r or r in p.parents:
            raise Refused(f"{p} lies inside the repository or the reference checkout ({r}): outputs and downloads "
                          "go outside it")
    return p


def dump_json(obj, path: Path) -> None:
    path = outside_repo(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(obj, indent=1, sort_keys=True, default=_json_default) + "\n")
    tmp.replace(path)


def _json_default(o):
    import numpy as np

    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


# ---------------------------------------------------------------------------------------------------------------------
# Stage cap


_RANGE = re.compile(r"e?\s*\d+(?:\.\d+)?\s*(?:-|–|to)\s*e?\s*\d", re.IGNORECASE)


def stage_ok(label) -> float:
    """E-day of a stage label, refused unless it parses as one stage at or before E9.5 and outside the T1 window.
    A range ('E9.5-E10.5') or anything unparseable fails closed."""
    s = str(label).strip()
    if _RANGE.search(s):
        raise Refused(f"stage {label!r} is a range: refused (fail closed)")
    try:
        e = windows.parse_stage(s)
    except ValueError:
        raise Refused(f"stage {label!r} does not parse: refused (fail closed)") from None
    if e > STAGE_MAX or windows.in_window(e, "T1"):
        raise Refused(f"stage {label!r} is after E9.5: refused (the T1 and heart window applies to every cell)")
    return e


# ---------------------------------------------------------------------------------------------------------------------
# The dataset list

_TOP = {"version", "fixed", "notes", "datasets", "qc_floors"}
_KEYS = {"accession", "status", "status_reason", "tier", "study", "kind", "design", "control_relation", "chimera",
         "repository", "screen_series", "stages_used", "control_definition", "arms", "exclude_match",
         "exclude_match_on",
         "design_refuse", "chemistry", "expect", "sample_file_patterns", "files_note", "licence", "why", "risks",
         "verified", "screen_extra_refuse", "review_cleared", "cell_labels", "citation", "pinned",
         "pinned_series", "stage_fields", "series_files", "pinned_series_files", "control_series"}
_SHA = re.compile(r"[0-9a-f]{64}")
_GSM = re.compile(r"GSM\d+")
_GSE = re.compile(r"GSE\d+")
_PLAIN_NAME = re.compile(r"[A-Za-z0-9._+-]+")
_SUFFIX_KEY = re.compile(r"[-_.][A-Za-z0-9]+")
_ARCHIVE = re.compile(r"\.(tar\.gz|tgz|tar)$", re.IGNORECASE)
_FIELD = re.compile(r"[A-Za-z0-9_]+")
# the sample fields a stage is read from when a dataset lists no stage_fields: the characteristics alone (with the
# title fallback of fetch.sample_stage's default rule), exactly the rule from before stage_fields existed
DEFAULT_STAGE_FIELDS = ("characteristics_ch1",)


def check_stage_fields(acc, ds: dict, where: str = "datasets.yaml") -> tuple[str, ...]:
    """A dataset's stage_fields, checked: absent means DEFAULT_STAGE_FIELDS; else a non-empty list of distinct SAMPLE
    field names, each a plain repository key (ASCII letters, digits, underscores). fetch.sample_stage reads them."""
    sf = ds.get("stage_fields")
    if sf is None:
        return DEFAULT_STAGE_FIELDS
    if (not isinstance(sf, list) or not sf or not all(isinstance(f, str) and _FIELD.fullmatch(f) for f in sf)
            or len(set(sf)) != len(sf)):
        raise Refused(f"{where}: {acc}: stage_fields must be a non-empty list of distinct plain field names "
                      "(ASCII letters, digits, underscores)")
    return tuple(sf)


def stage_fields(ds: dict) -> tuple[str, ...]:
    """The sample fields one dataset's stage is read from (check_stage_fields)."""
    return check_stage_fields(ds.get("accession"), ds)


def unit_design_factor(ds: dict, sister_control: bool) -> float:
    """d of one knockout unit: the dataset's d, but at most SISTER_D when any of its controls comes from a sister
    series (stage 6)."""
    return min(float(ds["d"]), SISTER_D) if sister_control else float(ds["d"])


def design_factor(ds: dict) -> float:
    """d of the shrinkage: 1 for a germline or conditional knockout against littermate or same-capture controls, and
    for a chimera with a wild-type-chimera second difference; 0.5 for everything else (F0 CRISPR, heterozygous or
    same-series controls of unverified relation)."""
    if ds["design"] == "chimera":
        return 1.0 if ds.get("chimera", {}).get("second_difference") else 0.5
    if ds["design"] in ("germline", "conditional") and ds["control_relation"] in ("littermate", "same_capture"):
        return 1.0
    return 0.5


def read_datasets_bytes(path: Path | None = None) -> tuple[Path, bytes]:
    """The YAML's path and bytes. The committed file (path None, or a path that resolves to it) is refused unless its
    sha256 is DATASETS_SHA256; another path is read as given (tests only: no CLI passes one)."""
    p = Path(path) if path is not None else DATASETS
    raw = p.read_bytes()
    if p.resolve() == DATASETS.resolve():
        got = hashlib.sha256(raw).hexdigest()
        if got != DATASETS_SHA256:
            raise Refused(f"{p}: sha256 {got} differs from t3x_common.DATASETS_SHA256 {DATASETS_SHA256}: the "
                          "pin-and-commit step must update both in the same commit (README, 'The pinned sample "
                          "list')")
    return p, raw


def datasets_sha256(path: Path | None = None) -> str:
    """sha256 of the dataset list that is read (the committed one is checked against its pin first)."""
    return hashlib.sha256(read_datasets_bytes(path)[1]).hexdigest()


def check_pins(acc: str, ds: dict, where: str = "datasets.yaml") -> list[dict]:
    """The pinned block of one dataset, checked for form only (its hashes are checked against the plan's stored
    metadata by fetch.py and process.py): a list of {accession, arm, gene, stage, files, metadata_sha256}, exact keys,
    arm exactly 'ko' or 'control', a GSM accession, a stage in stages_used, a non-empty list of distinct file names,
    a 64-hex sha256. One sample has one arm, one stage, one file list and one hash however many genes it is pinned
    for, and a (sample, gene) pair is pinned once. Pins on a dataset that is not cleared, or that needs per-cell
    labels, are refused."""
    pins = ds.get("pinned")
    if pins is None:
        pins = []
    if not isinstance(pins, list):
        raise Refused(f"{where}: {acc}: pinned must be a list")
    if pins and (ds.get("status") != "cleared" or acc in NEVER_FETCH):
        raise Refused(f"{where}: {acc}: only a cleared dataset may hold pins")
    if pins and ds.get("cell_labels"):
        raise Refused(f"{where}: {acc}: per-cell genotype designs are excluded in this stage: no pins")
    used = {stage_ok(s) for s in ds.get("stages_used") or []}
    seen, per_sample = set(), {}
    from_series = {g for u in (ds.get("pinned_series_files") or []) if isinstance(u, dict)
                   for g in (u.get("samples") or []) if isinstance(g, str)}
    for p in pins:
        if not isinstance(p, dict) or sorted(p) != sorted(PIN_KEYS):
            raise Refused(f"{where}: {acc}: pin {p!r} must have exactly the keys {PIN_KEYS}")
        a = p["accession"]
        if not isinstance(a, str) or not _GSM.fullmatch(a):
            raise Refused(f"{where}: {acc}: pinned accession {a!r} is not a GSM accession")
        if p["arm"] not in PIN_ARMS:
            raise Refused(f"{where}: {acc}: {a}: arm {p['arm']!r} is not exactly one of {PIN_ARMS}")
        if not isinstance(p["gene"], str) or not p["gene"].strip():
            raise Refused(f"{where}: {acc}: {a}: gene must be named (a control names the knockout gene it controls)")
        if stage_ok(p["stage"]) not in used:
            raise Refused(f"{where}: {acc}: {a}: stage {p['stage']!r} is not in stages_used")
        if (not isinstance(p["files"], list)
                or not all(isinstance(f, str) and f and "/" not in f for f in p["files"])
                or len(set(p["files"])) != len(p["files"])):
            raise Refused(f"{where}: {acc}: {a}: files must be a non-empty list of distinct file names")
        # stage 6: a sample takes its cells from exactly one source, its own per-sample files or one pinned series file
        if a in from_series and p["files"]:
            raise Refused(f"{where}: {acc}: {a}: a sample named by a pinned series file is pinned with files: [] (its "
                          "cells come from that file only)")
        if a not in from_series and not p["files"]:
            raise Refused(f"{where}: {acc}: {a}: files must be a non-empty list of distinct file names (only a sample "
                          "named by a pinned series file has none)")
        if not isinstance(p["metadata_sha256"], str) or not _SHA.fullmatch(p["metadata_sha256"]):
            raise Refused(f"{where}: {acc}: {a}: metadata_sha256 must be 64 lower-case hex digits")
        if (a, p["gene"]) in seen:
            raise Refused(f"{where}: {acc}: {a} is pinned twice for gene {p['gene']}")
        seen.add((a, p["gene"]))
        same = (p["arm"], stage_ok(p["stage"]), tuple(p["files"]), p["metadata_sha256"])
        if per_sample.setdefault(a, same) != same:
            raise Refused(f"{where}: {acc}: {a} is pinned with two different arms, stages, file lists or hashes")
    check_series_pins(acc, ds, bool(pins), where)
    check_series_file_pins(acc, ds, pins, where)
    return pins


def check_series_pins(acc: str, ds: dict, required: bool, where: str = "datasets.yaml") -> list[dict]:
    """The pinned_series block of one dataset, checked for form (stage 4). With sample pins it is required: a list of
    {accession, series_metadata_sha256}, exact keys, one entry for each series of screen_series and no other, each hash
    64 lower-case hex digits. Without sample pins it must be empty. Its hashes are checked against the plan's series
    text by fetch.verify_pins (live at plan time, stored at fetch and process time)."""
    sp = ds.get("pinned_series")
    if sp is None:
        sp = []
    if not isinstance(sp, list):
        raise Refused(f"{where}: {acc}: pinned_series must be a list")
    if not required:
        if sp:
            raise Refused(f"{where}: {acc}: pinned_series without sample pins: the readers pin both together")
        return sp
    for e in sp:
        if not isinstance(e, dict) or sorted(e) != sorted(SERIES_PIN_KEYS):
            raise Refused(f"{where}: {acc}: series pin {e!r} must have exactly the keys {SERIES_PIN_KEYS}")
        if not isinstance(e["accession"], str) or not _GSE.fullmatch(e["accession"]):
            raise Refused(f"{where}: {acc}: series pin accession {e['accession']!r} is not a GEO series")
        if not isinstance(e["series_metadata_sha256"], str) or not _SHA.fullmatch(e["series_metadata_sha256"]):
            raise Refused(f"{where}: {acc}: {e['accession']}: series_metadata_sha256 must be 64 lower-case hex digits")
    got = [e["accession"] for e in sp]
    want = all_screened(ds)
    if len(got) != len(set(got)) or set(got) != set(want):
        raise Refused(f"{where}: {acc}: pinned_series must name every screened series exactly once ({want}), not "
                      f"{got}: the readers pin the series text together with the samples")
    return sp


def control_series(ds: dict) -> dict | None:
    """The dataset's sister control series when it is declared AND named (accession set), else None (absent, or the
    placeholder with accession null, which is inert)."""
    c = ds.get("control_series")
    return c if isinstance(c, dict) and c.get("accession") else None


def all_screened(ds: dict) -> list[str]:
    """Every series the dataset screens, shows, hashes and pins: screen_series, then a named control series'
    screen_series (stage 6)."""
    out = list(ds.get("screen_series") or [])
    c = control_series(ds)
    for s in (c.get("screen_series") or []) if c else []:
        if s not in out:
            out.append(s)
    return out


def sister_series(ds: dict) -> list[str]:
    """The series screened for the named control series only (its screen_series), or []."""
    c = control_series(ds)
    return list(c.get("screen_series") or []) if c else []


def check_control_series(acc: str, ds: dict, where: str = "datasets.yaml") -> dict | None:
    """The control_series declaration, checked for form (stage 6). Absent or null: none. Else a mapping with exactly
    {accession, screen_series, relation, reader_note}. With accession null it is a placeholder: inert, its
    screen_series empty, and no sample can be pinned from it. Named: a GEO series other than the dataset's own, not
    among screen_series, listed first in its own screen_series (which shares no series with the dataset's
    screen_series); relation exactly {same_lab, same_protocol, same_publication}, each a statement; and a reader
    note; none of them may still hold the placeholder text."""
    c = ds.get("control_series")
    if c is None:
        return None
    if not isinstance(c, dict) or sorted(c) != sorted(CONTROL_SERIES_KEYS):
        raise Refused(f"{where}: {acc}: control_series must have exactly the keys {CONTROL_SERIES_KEYS}")
    rel = c["relation"]
    if not isinstance(rel, dict) or sorted(rel) != sorted(RELATION_KEYS):
        raise Refused(f"{where}: {acc}: control_series.relation must state exactly {RELATION_KEYS}")
    if not isinstance(c["screen_series"], list):
        raise Refused(f"{where}: {acc}: control_series.screen_series must be a list")
    if c["accession"] is None:
        if c["screen_series"]:
            raise Refused(f"{where}: {acc}: a control_series placeholder (accession null) screens nothing")
        return None
    s = c["accession"]
    if not isinstance(s, str) or not _GSE.fullmatch(s) or s == acc or s in (ds.get("screen_series") or []):
        raise Refused(f"{where}: {acc}: control_series.accession {s!r} must be a GEO series other than the dataset's "
                      "own screened series")
    if not c["screen_series"] or c["screen_series"][0] != s or not all(
            isinstance(x, str) and _GSE.fullmatch(x) for x in c["screen_series"]):
        raise Refused(f"{where}: {acc}: control_series.screen_series must list GEO series, its own accession first")
    if len(set(c["screen_series"])) != len(c["screen_series"]) or set(c["screen_series"]) & set(
            ds.get("screen_series") or []):
        raise Refused(f"{where}: {acc}: control_series.screen_series must be distinct and share no series with the "
                      "dataset's screen_series")
    for k in RELATION_KEYS:
        v = rel[k]
        if not isinstance(v, str) or not v.strip() or PLACEHOLDER in v:
            raise Refused(f"{where}: {acc}: control_series.relation.{k} must be a statement the readers checked")
    if not isinstance(c["reader_note"], str) or not c["reader_note"].strip() or PLACEHOLDER in c["reader_note"]:
        raise Refused(f"{where}: {acc}: control_series.reader_note is required (what both readers checked)")
    return c


def check_series_file_patterns(acc: str, ds: dict, where: str = "datasets.yaml") -> list[dict]:
    """The series_files declaration (stage 6): a list of {series, patterns}, each series one the dataset screens (its
    own or a named control series'), named once, each pattern a regular expression that compiles."""
    sf = ds.get("series_files")
    if sf is None:
        return []
    if not isinstance(sf, list):
        raise Refused(f"{where}: {acc}: series_files must be a list")
    seen = set()
    screened = set(all_screened(ds))
    for e in sf:
        if not isinstance(e, dict) or sorted(e) != ["patterns", "series"]:
            raise Refused(f"{where}: {acc}: series_files entry {e!r} must have exactly series and patterns")
        if not isinstance(e["series"], str) or e["series"] not in screened or e["series"] in seen:
            raise Refused(f"{where}: {acc}: series_files names {e['series']!r}: it must be a series the dataset "
                          "screens, named once")
        seen.add(e["series"])
        if not isinstance(e["patterns"], list) or not e["patterns"]:
            raise Refused(f"{where}: {acc}: {e['series']}: series_files patterns must be a non-empty list")
        for pat in e["patterns"]:
            try:
                re.compile(pat)
            except (re.error, TypeError):
                raise Refused(f"{where}: {acc}: {e['series']}: pattern {pat!r} does not compile") from None
    return sf


def triple_kind(name: str) -> str | None:
    """'matrix', 'barcodes' or 'features' when a file or member name is part of a 10x triple (the kind's name alone
    or after '_', '-' or '.', optionally gzipped; genes.tsv counts as features), else None."""
    low = name.lower()
    for kind, keys in (("matrix", ("matrix.mtx",)), ("barcodes", ("barcodes.tsv",)),
                       ("features", ("features.tsv", "genes.tsv"))):
        for k in keys:
            for end in (k, k + ".gz"):
                if low == end or (low.endswith(end) and low[-len(end) - 1] in "_-."):
                    return kind
    return None


def series_file_form(files) -> str:
    """'archive' (one .tar.gz, .tgz or .tar) or 'triple' (exactly a matrix, a barcodes and a features or genes file),
    else Refused."""
    if not isinstance(files, list) or not files or not all(
            isinstance(f, str) and _PLAIN_NAME.fullmatch(f) and ".." not in f for f in files):
        raise Refused(f"series file {files!r}: a non-empty list of plain file names")
    if len(files) == 1 and _ARCHIVE.search(files[0]):
        return "archive"
    if len(files) == 3 and sorted(triple_kind(f) or "" for f in files) == ["barcodes", "features", "matrix"]:
        return "triple"
    raise Refused(f"series file {files}: one archive (.tar.gz, .tgz, .tar) or one 10x triple (matrix, barcodes, "
                  "features or genes), nothing else")


def check_series_file_pins(acc: str, ds: dict, pins: list[dict], where: str = "datasets.yaml") -> list[dict]:
    """The pinned_series_files block (stage 6), checked for form. Each entry {series, files, samples, split}: a series
    that series_files declares; files one archive or one loose 10x triple, no file in two entries; samples pinned
    samples of this dataset, none in two entries (each is pinned with files: [], check_pins); split {by: none} or
    {by: barcode_suffix | folder, map: {key: GSM}} whose map names every sample of the entry and no other (a suffix
    is one separator then letters or digits, '-1'; a folder is one plain path component; a folder split needs an
    archive). With by none every sample of the entry is pinned with ONE arm, ONE gene and ONE stage: the file is then
    used as one unit and counts as one library."""
    units = ds.get("pinned_series_files")
    if units is None:
        units = []
    if not isinstance(units, list):
        raise Refused(f"{where}: {acc}: pinned_series_files must be a list")
    if units and not pins:
        raise Refused(f"{where}: {acc}: pinned_series_files without sample pins")
    declared = {e["series"] for e in check_series_file_patterns(acc, ds, where)}
    by_gsm: dict[str, list[dict]] = {}
    for p in pins:
        by_gsm.setdefault(p["accession"], []).append(p)
    seen_files, seen_samples = set(), set()
    for u in units:
        if not isinstance(u, dict) or sorted(u) != sorted(SERIES_FILE_UNIT_KEYS):
            raise Refused(f"{where}: {acc}: series file pin {u!r} must have exactly the keys {SERIES_FILE_UNIT_KEYS}")
        if u["series"] not in declared:
            raise Refused(f"{where}: {acc}: a series file pin names series {u['series']!r}, which series_files does "
                          "not declare")
        form = series_file_form(u["files"])
        if set(u["files"]) & seen_files:
            raise Refused(f"{where}: {acc}: {sorted(set(u['files']) & seen_files)} pinned in two series file entries")
        seen_files |= set(u["files"])
        sm = u["samples"]
        if (not isinstance(sm, list) or not sm or not all(isinstance(g, str) and _GSM.fullmatch(g) for g in sm)
                or len(set(sm)) != len(sm)):
            raise Refused(f"{where}: {acc}: {u['files']}: samples must be a non-empty list of distinct GSM accessions")
        for g in sm:
            if g not in by_gsm:
                raise Refused(f"{where}: {acc}: {u['files']}: {g} is not a pinned sample")
            if g in seen_samples:
                raise Refused(f"{where}: {acc}: {g} is named by two series file entries")
            seen_samples.add(g)
        sp = u["split"]
        if not isinstance(sp, dict) or sp.get("by") not in SPLITS:
            raise Refused(f"{where}: {acc}: {u['files']}: split.by must be one of {SPLITS}")
        if sp["by"] == "none":
            if sorted(sp) != ["by"]:
                raise Refused(f"{where}: {acc}: {u['files']}: split by none takes no map")
            kinds = {(p["arm"], p["gene"], stage_ok(p["stage"])) for g in sm for p in by_gsm[g]}
            if len(kinds) != 1:
                raise Refused(f"{where}: {acc}: {u['files']}: a file used as one unit holds samples of ONE arm, gene "
                              f"and stage only, not {sorted(kinds, key=str)}: state a split rule or do not pin it")
            continue
        if sorted(sp) != ["by", "map"] or not isinstance(sp["map"], dict) or not sp["map"]:
            raise Refused(f"{where}: {acc}: {u['files']}: split by {sp['by']} needs a non-empty map {{key: GSM}}")
        if sp["by"] == "folder" and form != "archive":
            raise Refused(f"{where}: {acc}: {u['files']}: a folder split needs an archive")
        for k in sp["map"]:
            ok = isinstance(k, str) and (_SUFFIX_KEY.fullmatch(k) if sp["by"] == "barcode_suffix"
                                         else (_PLAIN_NAME.fullmatch(k) and k not in (".", "..")))
            if not ok:
                raise Refused(f"{where}: {acc}: {u['files']}: split key {k!r} is not a plain "
                              f"{'barcode suffix (-1)' if sp['by'] == 'barcode_suffix' else 'folder name'}")
        if set(sp["map"].values()) != set(sm):
            raise Refused(f"{where}: {acc}: {u['files']}: the split map must name every sample of the entry and no "
                          "other")
    return units


def pinned_units(pins: list[dict]) -> list[dict]:
    """The knockout units the pins define, in pin order: {gene, ko: [GSM], control: [GSM], stages: [E-day]}. A
    knockout gene needs, at each stage of its knockout pins, a control pinned for that gene at that stage, and every
    control pin needs a knockout pin of its gene at its stage; otherwise refused."""
    units: dict[str, dict] = {}
    for p in pins:
        u = units.setdefault(p["gene"], {"gene": p["gene"], "ko": [], "control": [], "_st": {"ko": set(),
                                                                                             "control": set()}})
        u[p["arm"]].append(p["accession"])
        u["_st"][p["arm"]].add(stage_ok(p["stage"]))
    out = []
    for g, u in units.items():
        st = u.pop("_st")
        if not u["ko"]:
            raise Refused(f"gene {g}: control pins without a knockout pin")
        if st["ko"] != st["control"]:
            raise Refused(f"gene {g}: knockout stages {sorted(st['ko'])} and control stages {sorted(st['control'])} "
                          "differ: every stage needs both arms")
        out.append({**u, "stages": sorted(st["ko"])})
    return out


def load_datasets(path: Path | None = None) -> dict[str, dict]:
    """datasets.yaml, checked: known keys only (nothing can switch the screen off), known statuses, tiers, designs
    and relations, stages that pass the cap, the NEVER_FETCH accessions held or excluded, the pinned blocks' form, no
    sample pinned under two datasets. The committed file must hash to DATASETS_SHA256. Returns {accession: record}
    with 'd' added."""
    path, raw = read_datasets_bytes(path)
    doc = yaml.safe_load(raw.decode("utf-8"))
    if not isinstance(doc, dict) or set(doc) - _TOP:
        raise Refused(f"{path}: unknown top-level keys {sorted(set(doc or {}) - _TOP)}")
    out = {}
    pinned_in: dict[str, str] = {}
    for ds in doc.get("datasets") or []:
        acc = ds.get("accession")
        bad = set(ds) - _KEYS
        if bad:
            raise Refused(f"{path}: {acc}: unknown keys {sorted(bad)} (the screen has no switch)")
        if ds.get("status") not in STATUSES:
            raise Refused(f"{path}: {acc}: status {ds.get('status')!r} not in {STATUSES}")
        if acc in NEVER_FETCH and ds.get("status") not in ("held_phenocopy", "excluded_design"):
            raise Refused(f"{path}: {acc} is in t3x_common.NEVER_FETCH: its status must stay held_phenocopy or "
                          f"excluded_design, not {ds.get('status')!r}")
        if ds.get("tier") not in TIERS:
            raise Refused(f"{path}: {acc}: tier {ds.get('tier')!r} not in {TIERS}")
        if ds.get("design") not in DESIGNS or ds.get("control_relation") not in RELATIONS:
            raise Refused(f"{path}: {acc}: design/control_relation not in {DESIGNS}/{RELATIONS}")
        if ds.get("exclude_match_on", "text") not in ("text", "title"):
            raise Refused(f"{path}: {acc}: exclude_match_on must be text or title")
        if ds["design"] == "chimera" and ds["status"] == "cleared" and not ds.get("chimera", {}).get(
                "second_difference"):
            raise Refused(f"{path}: {acc}: a chimera without a wild-type-chimera second difference is not cleared")
        for s in ds.get("stages_used") or []:
            stage_ok(s)
        check_stage_fields(acc, ds, str(path))
        if acc in out:
            raise Refused(f"{path}: {acc} listed twice")
        if not isinstance(ds.get("screen_series"), list) or acc not in ds["screen_series"]:
            raise Refused(f"{path}: {acc}: screen_series must be a list that holds the dataset's own accession")
        from t3x_screen import check_review_cleared

        check_review_cleared(ds.get("review_cleared") or [])
        check_control_series(acc, ds, str(path))
        check_series_file_patterns(acc, ds, str(path))
        pins = check_pins(acc, ds, str(path))
        if pins:
            pinned_units(pins)
        for p in pins:
            other = pinned_in.setdefault(p["accession"], acc)
            if other != acc:
                raise Refused(f"{path}: {p['accession']} is pinned under {other} and {acc}: a sample belongs to one "
                              "series")
        out[acc] = {**ds, "pinned": list(pins), "pinned_series": list(ds.get("pinned_series") or []),
                    "pinned_series_files": list(ds.get("pinned_series_files") or []),
                    "series_files": list(ds.get("series_files") or []), "d": design_factor(ds)}
    held = held_series(out)
    for acc, ds in out.items():
        bad = sorted(set(all_screened(ds)) & held)
        if ds["status"] == "cleared" and bad:
            raise Refused(f"{path}: {acc}: screen_series includes the held or excluded series {bad}: refused (a "
                          "cleared dataset never screens, fetches or reads a held series)")
        ds["held_series"] = sorted(held)
    return out


def held_series(allds: dict[str, dict]) -> set[str]:
    """The held set: NEVER_FETCH plus the accession and every screened series of every dataset that is not cleared.
    A pinned sample whose record lists any of them is refused (fetch.verify_pins), and a plan whose screened series
    hold such a sample is refused with its text withheld."""
    held = set(NEVER_FETCH)
    for acc, ds in allds.items():
        if ds.get("status") != "cleared":
            held.add(acc)
            held |= set(ds.get("screen_series") or [])
            c = ds.get("control_series")   # stage 6: a held dataset's sister series is held too
            if isinstance(c, dict):
                held |= {x for x in [c.get("accession"), *(c.get("screen_series") or [])] if isinstance(x, str)}
    return held


def qc_floors(path: Path | None = None) -> dict:
    """datasets.yaml's qc_floors, checked; refused when missing (a table is never built without its floor)."""
    path, raw = read_datasets_bytes(path)
    doc = yaml.safe_load(raw.decode("utf-8")) or {}
    fl = doc.get("qc_floors")
    if not isinstance(fl, dict) or set(fl) != {"min_confident_share_per_sample"}:
        raise Refused(f"{path}: qc_floors must hold exactly min_confident_share_per_sample")
    v = fl["min_confident_share_per_sample"]
    if not isinstance(v, int | float) or not 0 < float(v) <= 1:
        raise Refused(f"{path}: min_confident_share_per_sample must be in (0, 1]")
    return {"min_confident_share_per_sample": float(v)}


def dataset(acc: str, path: Path | None = None, need_cleared: bool = True) -> dict:
    """One record; refused for a NEVER_FETCH accession whatever the YAML says, if the accession is not listed, or if
    it is not cleared (when need_cleared)."""
    if acc in NEVER_FETCH:
        raise Refused(f"{acc} is in t3x_common.NEVER_FETCH (held or excluded): refused whatever a YAML says")
    allds = load_datasets(path)
    if acc not in allds:
        raise Refused(f"{acc} is not in datasets.yaml: refused")
    ds = allds[acc]
    if need_cleared and ds["status"] != "cleared":
        raise Refused(f"{acc} has status {ds['status']} ({ds.get('status_reason', '')}): refused")
    return ds


# ---------------------------------------------------------------------------------------------------------------------
# The canonical serialisation of one sample's screened metadata (README, "The pinned sample list")


def canonical_metadata(gsm: str, fields: dict) -> bytes:
    """The bytes whose sha256 is a pin's metadata_sha256: the UTF-8 encoding of
        json.dumps({"accession": gsm, "fields": F}, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    where F maps every field name of the sample's GEO family SOFT record (the 'Sample_' prefix removed) to the list of
    its values, each a string exactly as fetch.parse_soft gives it (surrounding whitespace stripped, nothing else
    changed), in the order of the SOFT file, plus the key '_other_lines' when the record holds lines that are not
    'key = value' fields (parse_soft keeps them). Every field is included because the screen reads every field. A
    series is serialised the same way (series_metadata_sha256): its accession, and F its fields with the 'Series_'
    prefix removed, contact fields included."""
    if not isinstance(fields, dict) or not all(
            isinstance(k, str) and isinstance(v, list) and all(isinstance(x, str) for x in v)
            for k, v in fields.items()):
        raise Refused(f"{gsm}: metadata fields must map names to lists of strings")
    return json.dumps({"accession": gsm, "fields": fields}, sort_keys=True, ensure_ascii=False,
                      separators=(",", ":")).encode("utf-8")


def metadata_sha256(gsm: str, fields: dict) -> str:
    return hashlib.sha256(canonical_metadata(gsm, fields)).hexdigest()


def series_metadata_sha256(gse: str, fields: dict) -> str:
    """A series' hash: the same canonical serialisation as a sample's (canonical_metadata)."""
    return metadata_sha256(gse, fields)


def pins_sha256(ds: dict) -> str:
    """sha256 of one dataset's pins, canonically serialised: {"samples": pinned, "series": pinned_series}, plus (stage
    6, only when present, so a dataset without them keeps its hash) "series_files": pinned_series_files and
    "control_series": the named control_series declaration. The fetch manifest and each of its file records keep it;
    fetch downloads again and process refuses on any change."""
    rec = {"samples": ds.get("pinned") or [], "series": ds.get("pinned_series") or []}
    if ds.get("pinned_series_files"):
        rec["series_files"] = ds["pinned_series_files"]
    if control_series(ds):
        rec["control_series"] = control_series(ds)
    return hashlib.sha256(json.dumps(rec, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
                          .encode("utf-8")).hexdigest()


def shrink_k(n_arms, m: int, d: float) -> float:
    """k = min(n_h / (n_h + 200), m / (m + 1)) * m / (m + 1) * d, with n_h the harmonic mean of the arms' cell
    numbers and m the smallest number of independent libraries (or pools) of any arm. The cell term is capped by the
    unit term because cells of one library are not independent (the challenger's change 3)."""
    n = [float(x) for x in n_arms]
    if min(n) <= 0 or m < 1:
        return 0.0
    nh = len(n) / sum(1.0 / x for x in n)
    unit = m / (m + 1.0)
    return min(nh / (nh + N_HALF), unit) * unit * float(d)


LOG_CLIP = (math.log(CLIP[0]), math.log(CLIP[1]))
