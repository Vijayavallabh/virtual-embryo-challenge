"""Shared pieces of the Task 1 outside-atlas program (experiments/t1_ext): the atlas list, the pinned block, the
canonical record text, the stage rule, paths.

Third version (2026-11-21): the PINS are the first line of defence. A sample is used only when atlases.yaml pins it with
the sha256 of its full verbatim record (record_text), after two readers have read that record line by line. The stage
rule is the second line and can only refuse. The rule is written and tested before any record is parsed (the
challenger's change on plug-day conventions, integer labels, somite counts and Theiler stages; the second version after
the two reviews of 2026-11-20 reads the whole record, fails closed on unparsed tokens, handles plug = E1 and widens
readings when no convention is stated; the third closes the second verifier's three blockers). A reading is an interval
in E-days. Anything any reading could put inside the T1 window (E9.5, E13.5] counts as inside. The canonical
processing text and its clearance live here too. Everything here is pure: no network, no data file.
"""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import re
import sys
import unicodedata
from collections.abc import Iterable, Mapping
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ATLASES = HERE / "atlases.yaml"
EXTERNAL = REPO / "experiments" / "external"
for _p in (EXTERNAL,):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
import windows  # noqa: E402

TASKS = ("T1",)
WINDOW = (9.5, 13.5)  # T1: low exclusive, high inclusive (windows.WINDOWS["T1"])
BIRTH = windows.BIRTH_EDAY
PENDING = "PENDING"
MARKS = ("USED", "WINDOW", "NOT_USED", "UNKNOWN")
STATUSES = ("in_use", "use", "second", "held", "dropped")  # dropped: by the readers' screens (row T72, cell 3)
FETCHABLE = ("use", "second")
# Row T72 (candidate MX2s), cells 3 and 8, as literal lists: no step requests any of these by any route, whatever
# atlases.yaml says (the yaml's never_fetch lists must name exactly these too; load_atlases checks it).
NEVER_FETCHED = {
    "GSE244783": ("GSM7829701", "GSM7829700", "GSM7829702", "GSM7829704", "GSE244783_RAW.tar"),
    "GSE230531": ("GSM7226279", "GSM7226280", "GSM7226270", "GSM7226271", "GSM7226272", "GSM7226273",
                  "GSE230531_RAW.tar"),
}
NEVER_FETCHED_ALL = frozenset(x for v in NEVER_FETCHED.values() for x in v)
# The six samples whose barcodes and features files alone `fetch.py fetch-meta` requests, before any matrix (cell 3).
SCREEN_SAMPLES = {"GSE244783": {"GSM7829699": "early", "GSM7829703": "late"},
                  "GSE230531": {"GSM7226268": "early", "GSM7226269": "early", "GSM7226274": "late",
                                "GSM7226276": "late"}}
# The register rows of the three members (cell 3: X002, X019 GSE244783, X020 GSE230531).
REGISTER_ROWS = {"GSE193346": "X002", "GSE244783": "X019", "GSE230531": "X020"}
# The screens' folder (cell 3: fetch-meta "writes only under ~/v/scratch/t72_screen, never into the row's data
# folders"). The plan, pin-template, fetch-meta and screens.py all use it.
SCREEN_ROOT = Path.home() / "v" / "scratch" / "t72_screen"
SCREEN_POINTS = ("1", "2", "3", "4", "5", "6", "7", "a", "b", "c", "d")
SCREEN_ROLES = ("barcodes", "features")
# The DGX's reference checkout stays exactly as rehearsed: nothing here writes into it.
REFERENCE_CHECKOUT = Path.home() / "v" / "virtual-embryo-challenge"


class Refused(RuntimeError):
    """A rule of atlases.yaml or of this program stops the step."""


# ------------------------------------------------------------------------------------------------ modules and files


def load_module(name: str, path: Path):
    """Import a module from a file under a unique name (this folder's fetch.py would otherwise collide with another
    folder's fetch.py on sys.path)."""
    key = f"t1_ext_mod_{name}"
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[key] = mod
    spec.loader.exec_module(mod)
    return mod


def sibling(name: str):
    return load_module(name, HERE / f"{name}.py")


def sha256_file(path, chunk: int = 1 << 23) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def sha256_json(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=_json_default).encode()).hexdigest()


def outside_repo(path, repo: Path | None = None) -> Path:
    """The resolved path, refused if it lies inside this checkout, the main checkout (when this is a worktree) or the
    DGX's reference checkout."""
    p = Path(path).expanduser().resolve()
    main = (repo or REPO).resolve()
    roots = {main, REFERENCE_CHECKOUT.expanduser().resolve()}
    if ".claude" in main.parts:
        roots.add(Path(*main.parts[:main.parts.index(".claude")]))
    for r in roots:
        if p == r or r in p.parents:
            raise Refused(f"{p} lies inside the repository or the reference checkout ({r}): outputs go outside it")
    return p


ROW_DATA_RELS = ("data/external/t1_ext_t72", "data/derived/t1_ext_t72")


def writable(path, repo: Path | None = None) -> Path:
    """The resolved path where a step of this folder may write: outside every checkout (outside_repo), or inside row
    T72's two output folders under the checkout's own data/ (a real directory; row_t72.guards checks the rest)."""
    p = Path(path).expanduser().resolve()
    r = Path(repo or REPO)
    data = r / "data"
    if data.is_dir() and not data.is_symlink():
        for rel in ROW_DATA_RELS:
            root = data.resolve() / rel.split("/", 1)[1]
            if p == root or root in p.parents:
                return p
    return outside_repo(path, repo)


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
    if isinstance(o, set | frozenset):
        return sorted(o)
    raise TypeError(type(o))


def dump_json(obj, path) -> Path:
    path = writable(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(obj, indent=1, sort_keys=True, default=_json_default) + "\n")
    tmp.replace(path)
    return path


# ------------------------------------------------------------------------------------------------ the atlas list


def load_atlases(path: Path | None = None) -> dict:
    cfg = yaml.safe_load(Path(path or ATLASES).read_text())
    if tuple(cfg.get("tasks") or ()) != TASKS:
        raise Refused("atlases.yaml must serve T1 only (tasks: [T1])")
    for acc, a in cfg["atlases"].items():
        if a.get("status") not in STATUSES:
            raise Refused(f"{acc}: status {a.get('status')!r} is not one of {STATUSES}")
        if a["status"] != "in_use":
            for k in ("series", "form", "early", "late", "groups_formed", "qc", "early_tissue"):
                if k not in a:
                    raise Refused(f"{acc}: {k} missing")
            if a["form"] not in cfg["forms"]:
                raise Refused(f"{acc}: unknown form {a['form']!r}")
            for s in a["series"]:
                if not re.fullmatch(r"GSE\d+", s):
                    raise Refused(f"{acc}: series {s!r} is not a GEO series accession")
        if a["status"] in FETCHABLE:
            # the strict form (ruling row R27): every fetchable atlas names its samples per end, its never-fetch list,
            # and, for 10x files, the screens run on barcodes and features before any matrix is fetched
            for e in ("early", "late"):
                if not a[e].get("gsms"):
                    raise Refused(f"{acc}: the {e} end names no samples (gsms): the strict form names every sample")
            if not isinstance(a.get("never_fetch"), Mapping):
                raise Refused(f"{acc}: never_fetch missing")
            if a["form"] == "tenx_mtx" and not a.get("screens"):
                raise Refused(f"{acc}: a 10x atlas without screens is refused")
        if acc in REGISTER_ROWS and a.get("register") != REGISTER_ROWS[acc]:
            raise Refused(f"{acc}: register must be {REGISTER_ROWS[acc]} (row T72's members)")
        if acc in NEVER_FETCHED and a.get("status") != "in_use":
            nf = a.get("never_fetch") or {}
            if set(nf.get("samples") or []) | set(nf.get("files") or []) != set(NEVER_FETCHED[acc]):
                raise Refused(f"{acc}: never_fetch must name exactly row T72's list {list(NEVER_FETCHED[acc])}")
            if a.get("status") in FETCHABLE and dict(named_gsms({**a, "accession": acc})) != SCREEN_SAMPLES[acc]:
                raise Refused(f"{acc}: the named samples must be exactly row T72's {SCREEN_SAMPLES[acc]}")
        check_never_fetch(acc, a)
        check_screen_block(acc, a)
        if a.get("screens"):
            unknown = set(a["screens"]) - set(SCREENS)
            if unknown:
                raise Refused(f"{acc}: unknown screens {sorted(unknown)} (known: {SCREENS})")
        check_pins(acc, a)
    return cfg


# ------------------------------------------------------------------------------------------------ the strict form
#
# Ruling row R27 reopened the line in its strict form: per atlas, the samples of each end are named (`gsms`), a pin may
# name only those, and the pins must cover every named sample; the never-fetch list names samples and files that no
# step may fetch (every sample of the series that is not named is also never fetched, since only pinned samples are).
SCREENS = ("suffix", "own_peaks", "no_hashing", "whitelist_sized")


def named_gsms(a: Mapping) -> dict[str, str]:
    """GSM -> end, for the samples the atlas entry names (empty when it names none: synthetic test entries)."""
    out: dict[str, str] = {}
    for e in ("early", "late"):
        for g in a.get(e, {}).get("gsms") or []:
            if g in out:
                raise Refused(f"{a.get('accession')}: {g} is named at both ends")
            out[g] = e
    return out


def check_screen_block(acc: str, a: Mapping) -> None:
    """The readers' screens as the frozen yaml records them (row T72, cell 3: "the frozen yaml records each atlas's
    status and its pins and the sha256 of each screened barcodes and features file"). Per atlas with screens:
      screen: {verdict: PENDING}                                    before the screens
      screen: {verdict: kept, files: {GSM: {barcodes: {name, sha256}, features: {name, sha256}}, ...},
               read_by: [two readers], date: "YYYY-MM-DD"}           status use: both readers said yes to every point
      screen: {verdict: dropped, failed: [points], files: {...}, read_by: [...], date: ...}   status dropped
    A kept atlas's files name exactly its named samples and the two roles, each with a plain file name and a sha256
    (as `screens.py` printed them); a dropped atlas names the points that failed (1 to 7, a to d) and is never
    fetched, loaded or tabled. No exception is granted on the spot."""
    sc = a.get("screen")
    if a.get("status") == "dropped":
        if not isinstance(sc, Mapping) or sc.get("verdict") != "dropped" or not sc.get("failed") or \
                not set(sc["failed"]) <= set(SCREEN_POINTS):
            raise Refused(f"{acc}: a dropped atlas records screen.verdict dropped and the failed points "
                          f"(of {SCREEN_POINTS})")
        return
    if not a.get("screens"):
        if sc is not None:
            raise Refused(f"{acc}: a screen record on an atlas without screens")
        return
    if a.get("status") not in FETCHABLE:
        return
    if not isinstance(sc, Mapping) or sc.get("verdict") not in ("PENDING", "kept"):
        raise Refused(f"{acc}: screen.verdict must be PENDING or kept (a failed screen sets status: dropped)")
    if sc["verdict"] == "PENDING":
        if set(sc) - {"verdict", "files"} or sc.get("files"):
            raise Refused(f"{acc}: a PENDING screen holds no files, readers or date")
        return
    if set(sc) != {"verdict", "files", "read_by", "date"}:
        raise Refused(f"{acc}: a kept screen holds exactly verdict, files, read_by and date")
    _check_readers(f"{acc}: screen", {**sc, "record_sha256": "0" * 64})
    files = sc["files"]
    named = named_gsms({**a, "accession": acc})
    if not isinstance(files, Mapping) or set(files) != set(named):
        raise Refused(f"{acc}: the screen's files must name exactly the atlas's samples {sorted(named)}")
    for g, f in files.items():
        if not isinstance(f, Mapping) or set(f) != set(SCREEN_ROLES):
            raise Refused(f"{acc}: {g}: the screen records exactly the roles {SCREEN_ROLES}")
        for role, x in f.items():
            if not isinstance(x, Mapping) or set(x) != {"name", "sha256"} or not isinstance(x["sha256"], str) or \
                    not _SHA.fullmatch(x["sha256"]) or not isinstance(x["name"], str) or \
                    not re.fullmatch(rf"{g}_[A-Za-z0-9._+-]+", x["name"]):
                raise Refused(f"{acc}: {g} {role}: the screen records a plain file name named after the sample and "
                              f"its sha256 (64 lower-case hex digits)")


def screen_kept(acc: str, a: Mapping) -> dict | None:
    """For an atlas with screens: its screen block, which must say kept (else refused): {GSM: {role: {name,
    sha256}}}. None for an atlas without screens (a synthetic test entry; load_atlases refuses a fetchable 10x atlas
    without them)."""
    if not a.get("screens"):
        return None
    sc = a.get("screen") or {}
    if sc.get("verdict") != "kept":
        raise Refused(f"{acc}: atlases.yaml records the screen as {sc.get('verdict')!r}: its matrices are not "
                      f"fetched, loaded or tabled")
    return {g: {r: dict(x) for r, x in f.items()} for g, f in sc["files"].items()}


def never_fetched_refusal(*names: str | None) -> str | None:
    """Why a request naming any of `names` (accessions, file names, URLs) is refused under row T72's literal list."""
    for n in names:
        if not n:
            continue
        for x in sorted(NEVER_FETCHED_ALL):
            if x.endswith(".tar"):
                if x.lower() in str(n).lower():
                    return f"{x} is on row T72's never-fetched list"
            elif re.search(rf"(?<!\d){x}(?!\d)", str(n)):
                return f"{x} is on row T72's never-fetched list"
        if re.search(r"_raw\.tar", str(n), re.IGNORECASE):
            return f"{n}: a series archive (RAW.tar) is never fetched"
    return None


def check_never_fetch(acc: str, a: Mapping) -> None:
    nf = a.get("never_fetch")
    if nf is None:
        return
    if not isinstance(nf, Mapping) or set(nf) - {"samples", "files"}:
        raise Refused(f"{acc}: never_fetch holds only samples and files")
    for g in nf.get("samples") or []:
        if not isinstance(g, str) or not re.fullmatch(r"GSM\d+", g):
            raise Refused(f"{acc}: never_fetch sample {g!r} is not a GSM accession")
    clash = set(nf.get("samples") or []) & set(named_gsms({**a, "accession": acc}))
    if clash:
        raise Refused(f"{acc}: {sorted(clash)} are both named and never fetched")


def never_fetch_refusal(a: Mapping, gsm: str | None = None, url: str | None = None) -> str | None:
    """Why a sample or a file may not be fetched under the atlas's never-fetch list (None when it may): a listed sample
    by accession, a listed file by name, any file whose name or path names a listed sample, and every series archive
    (<GSE>_RAW.tar) of the atlas's series by name."""
    lit = never_fetched_refusal(gsm, url)
    if lit:
        return lit
    nf = a.get("never_fetch") or {}
    samples = set(nf.get("samples") or [])
    files = set(nf.get("files") or []) | {f"{s}_RAW.tar" for s in a.get("series", [])}
    if gsm is not None and gsm in samples:
        return f"{gsm} is on the never-fetch list"
    if url is not None:
        name = file_name(url)
        if name in files:
            return f"{name} is on the never-fetch list (a series archive or a named file)"
        for g in samples:
            if re.search(rf"(?<!\d){g}(?!\d)", url):
                return f"{name} names {g}, which is on the never-fetch list"
    return None


def pins_complete(a: Mapping, pins: Mapping) -> list[str]:
    """Refusals when the atlas names its samples and the pins do not pin exactly those samples at their ends."""
    named = named_gsms(a)
    if not named:
        return []
    out = []
    missing = sorted(set(named) - set(pins["samples"]))
    if missing:
        out.append(f"named samples not pinned: {missing} (the strict form uses every named sample or none)")
    for g, p in pins["samples"].items():
        if named.get(g) != p["end"]:
            out.append(f"{g} is pinned at the {p['end']} end but the atlas names it "
                       f"{'at the ' + named[g] + ' end' if g in named else 'nowhere'}")
    return out


def atlas(acc: str, cfg: Mapping | None = None) -> dict:
    """One atlas entry; an accession not under `atlases` is refused."""
    cfg = cfg or load_atlases()
    if acc not in cfg["atlases"]:
        raise Refused(f"{acc} is not in experiments/t1_ext/atlases.yaml: refused")
    a = copy.deepcopy(cfg["atlases"][acc])
    a["accession"] = acc
    return a


def entry_sha256(acc: str, cfg: Mapping | None = None) -> str:
    """The sha256 of one atlas entry and the file-wide rules it uses (forms, protocol and enrichment terms, the stage
    tables, the plug-day conventions, the canonical forms), so that pinning the derived extra classes later does not
    invalidate a plan. The atlas's pinned block is part of the entry: a plan made under other pins is not fetched. The
    atlas's processing_cleared entry is left out: fetch and load check it against the text itself."""
    cfg = cfg or load_atlases()
    a = atlas(acc, cfg)
    a.pop("processing_cleared", None)  # checked on its own, against the live text, by fetch and load
    return sha256_json({"atlas": a, "forms": cfg["forms"], "protocol_terms": cfg["protocol_terms"],
                        "enrichment_terms": cfg["enrichment_terms"], "theiler": THEILER, "somites": SOMITES,
                        "conventions": CONVENTIONS, "unstated": UNSTATED, "window": WINDOW,
                        "processing_canonical": PROCESSING_FORMAT, "record_canonical": RECORD_FORMAT,
                        "screen_rules": cfg.get("screen_rules"), "series_text_rule": SERIES_TEXT_RULE,
                        "pooling_words": POOLING_WORDS.pattern,
                        "pooling_words_staged": POOLING_WORDS_STAGED.pattern,
                        "pooling_stage_words": POOLING_STAGE_WORDS.pattern})


def require_frozen(cfg: Mapping, repo: Path | None = None) -> None:
    """Load and tables run only on a frozen atlases.yaml whose derived extra classes are pinned and committed."""
    if cfg.get("status") != "frozen":
        raise Refused("atlases.yaml is not frozen: the main session freezes it after reading every plan table")
    ex = cfg["extra_competing"]
    if ex.get("derived_sha256") in (None, PENDING):
        raise Refused("extra_competing.derived_sha256 is PENDING: run `groups.py derive` and pin its output first")
    f = Path(repo or REPO) / ex["derived_file"]
    if not f.exists() or sha256_file(f) != ex["derived_sha256"]:
        raise Refused(f"{ex['derived_file']} is missing or does not match extra_competing.derived_sha256")


# ------------------------------------------------------------------------------------------------ the stage rule
#
# The stage rule is the SECOND line of defence (third version, 2026-11-21). The first line is the pinned block of
# atlases.yaml: a sample is used only when two readers have read its full record and pinned it (README, "Pinned
# samples"). The rule below can only refuse: it never makes a sample usable that is not pinned.
#
# A sample's stage comes from the WHOLE of its own record (record_stage): every field of its repository record except
# the few listed in NOT_STAGE (dates, contacts, the platform, the row count, unparsed lines), read in one of two modes,
# plus every supplementary file name. Every reading is an interval in E-days. The rule, in order:
#   1. The plug-day convention is read from every sentence of the series' summary and design and of the sample's own
#      protocol, processing and description text that names a plug, mating, coitus, copulation, conception or
#      fertilisation (convention). E0.5, E0 and E1 (day 1) are recognised; any other such sentence naming a day, or two
#      different conventions, gives UNKNOWN. Each reading is widened by its convention's shift interval (CONVENTIONS;
#      UNSTATED when none is stated).
#   2. Descriptive fields (title, source name, characteristics, description, relations, labels, ...) are read whole,
#      value by value; the sample's own protocol and processing fields sentence by sentence, with the convention's own
#      tokens blanked in convention sentences. A reading in protocol text counts exactly as a reading anywhere else.
#   3. The readings of one value (or one protocol sentence) join into the span they cover, so pooled or listed stages
#      ("hearts pooled from E12.5 and E14.5 embryos", "12.5 and 14.5 dpc") read as the interval between them.
#   4. A stage-like token left over after the patterns have consumed theirs (LEFTOVERS, a half-day number no pattern
#      read, and, in a value or sentence that gave a reading, a companion integer from 5 to 21) gives UNKNOWN: an
#      unparsed form never drops silently. A bare number in a stage-keyed field is read as an E-day.
#   5. No reading at all gives UNKNOWN; readings that do not overlap give UNKNOWN.
#   6. Otherwise the union of the widened readings decides: WINDOW when it touches (E9.5, E13.5], else PERMITTED.
#   window_reading is True whenever any widened reading touches the window, whatever the mark; such a sample keeps only
#   the fields that gave a reading (fetch.classify).
# All text is NFKC-normalised before it is read (a full-width E reads as E); the record itself is hashed verbatim.

# Theiler stage -> E-day interval: the EMAP ranges as the AI assistant recalls them, widened by up to half a day; NOT
# checked against a source. They can only widen a sample's interval.
THEILER: dict[int, tuple[float, float]] = {
    9: (6.25, 7.25), 10: (6.5, 7.75), 11: (7.25, 8.0), 12: (7.5, 8.75), 13: (8.0, 9.25), 14: (8.5, 9.75),
    15: (9.0, 10.25), 16: (9.5, 10.75), 17: (10.0, 11.25), 18: (10.5, 11.25), 19: (11.0, 12.25), 20: (11.5, 13.0),
    21: (12.5, 14.0), 22: (13.5, 15.0), 23: (14.5, 15.5), 24: (15.5, 16.5), 25: (16.5, 17.5), 26: (17.5, 18.5),
    27: (18.5, 19.5), 28: (19.5, 30.0)}
# Somite count ranges -> E-day interval, general knowledge widened, NOT checked against a source. The last row ends at
# 15.0 so that a high count can never read narrower than the window (the first review's point 3).
SOMITES: tuple[tuple[int, int, float, float], ...] = (
    (0, 7, 7.5, 8.5), (8, 12, 8.0, 9.0), (13, 20, 8.5, 9.5), (21, 29, 9.0, 10.25), (30, 34, 9.75, 10.75),
    (35, 39, 10.25, 11.25), (40, 99, 10.75, 15.0))

# Plug-day conventions: the shift interval, in days, from an author's embryonic label onto the plug = E0.5 scale (the
# noon of the plug day is E0.5) on which the window is written.
#   E0.5  labels as written.
#   E0    the plug day is E0: a label counts from the plug day's midnight (as E0.5 does: no shift) or from its noon
#         (+0.5); both are possible, so [0, +0.5].
#   E1    the plug day is E1 (day 1): one day later than E0's count, so [-1.0, -0.5]. E14.5 then reads
#         [13.5, 14.0] and meets the window's closed upper edge.
#   UNSTATED  no sentence states one: any of the three is possible, so the union [-1.0, +0.5]. A half-day label n.5
#         then reads [n - 0.5, n + 1.0]: E14.5 meets 13.5 and E9.5 passes 9.5, so both are WINDOW unless the record
#         states E0.5 (or, for E14.5, E0). E8.5 ([7.5, 9.0]) and E16.5 ([15.5, 17.0]) stay outside.
CONVENTIONS: dict[str, tuple[float, float]] = {"E0.5": (0.0, 0.0), "E0": (0.0, 0.5), "E1": (-1.0, -0.5)}
UNSTATED = (-1.0, 0.5)
# Characteristics keys whose bare numbers are read as E-days (brief: "a bare number in a stage-keyed field").
STAGE_KEY = re.compile(r"(?<![a-z])(?:stage|age|time|day|develop|gestation|dpc|e-?day|embryonic|theiler|somite|"
                       r"week|month|hpf|pnd|postnatal)", re.IGNORECASE)
# A stage key naming a unit that is not an embryonic day: its bare number cannot be read, so it is a leftover.
_KEY_OTHER_UNIT = re.compile(r"week|month|year|postnatal|pnd|hour|hpf|\bp\b", re.IGNORECASE)

_NUM = r"(\d{1,2}(?:\.\d+)?)"
_END = r"(?!\d|\.\d)"   # a number ends here (a sentence's full stop may follow)
_L = r"(?:^|(?<=[^a-z0-9]))"  # a token starts the text or follows a non-alphanumeric character
_RX_THEILER = re.compile(_L + r"(?:ts|theiler stage|theiler)\s?(\d{1,2})(?!\d)")
_RX_SOMITE_RANGE = re.compile(_L + r"(\d{1,2})\s*(?:-|–|to)\s*(\d{1,2})\s*(?:somites?|ss)(?![a-z])")
_RX_SOMITE = re.compile(_L + r"(\d{1,2})\s*(?:somites?|ss)(?![a-z])|somites?\s*[:=]?\s*(\d{1,2})(?!\d)")
_RX_WEEKS = re.compile(_L + _NUM + r"(?:\s*(?:-|–|to)\s*" + _NUM + r")?\s*(weeks?|wks?|months?|mo)(?![a-z])")
_RX_POSTNATAL_WORD = re.compile(r"(?<![a-z])(?:old|postnatal|post-natal|after birth|of age|adult|pnd)(?![a-z])")
_RX_ADULT = re.compile(_L + r"adult(?![a-z])")
_RX_NEONATAL = re.compile(_L + r"(?:neonatal|newborn|neonate)s?(?![a-z])")
_RX_POSTNATAL = re.compile(_L + r"(?:postnatal day|postnatal|pnd|p)\s?(\d{1,2})" + _END)
# Embryonic day prefixes: E, ED, GD, e.d., embryonic day (also glued: EmbryonicDay12.5), dpc, pc (the verifier's
# blocker 2, (ii)).
_E = r"(?:e|ed|gd|e\.d\.?|embryonic ?day|dpc|d\.p\.c\.?|pc|p\.c\.?)"
_DPC = r"(?:dpc|d\.p\.c\.?|pc|p\.c\.?|days? post[- ]?coit(?:um|us|al))"
_RX_E_RANGE = re.compile(_L + _E + r"\s?" + _NUM + r"\s*(?:-|–|—|to|~)\s*" + _E + r"?\s?" + _NUM + _END)
# A range or list whose first number has no prefix, closed by a dpc unit: "13.5-14.5 dpc", "12.5 and 14.5 days
# post-coitum" (the verifier's blocker 1, (i)); read as the span [min, max], an embryonic range.
_RX_DPC_RANGE = re.compile(_L + _NUM + r"\s*(?:-|–|—|to|and|&|/|,)\s*" + _NUM + r"\s?" + _DPC + r"(?![a-z])")
_RX_DPC_AFTER = re.compile(_L + _NUM + r"\s?" + _DPC + r"(?![a-z])")
_RX_E_POINT = re.compile(_L + _E + r"\s?" + _NUM + _END)
# "day 12.5" also glued to letters (Heartday12.5); "12 days".
_RX_DAYS = re.compile(_L + _NUM + r"\s*days?(?![a-z])|day\s?" + _NUM + _END)
_RX_BARE = re.compile(_L + _NUM + _END)
# Blanked before any reading: GEO and SRA accessions, read lengths (PE150, SE75), scientific notation (1E6, 2.5e5).
_RX_BLANK = re.compile(r"(?:GS[EM]|SR[RXSP]|SAM[NDE]|PRJ[A-Z]{2})\d+|(?<![A-Za-z])[PS]E\d{2,3}(?!\d)|"
                       r"\d(?:\.\d+)?[Ee][+-]?\d+", re.IGNORECASE)
# A stage token spelt with an underscore, a p or a decimal comma (E12_5, e14p5, E12,5): normalised to E12.5 (same
# length). A capital E may be glued to the preceding letters (HeartE12_5); a lower-case e must follow a non-letter.
_RX_SPELT = re.compile(r"(?:(?<![A-Za-z])[Ee]|E)(\d{1,2})[,_p](\d)(?!\d)")
_RX_GLUED = re.compile(r"(?<=[A-Za-z0-9])(?=E\d|TS\d)")
# A unit after a number: the number is a quantity, not a stage (0.5 mg/ml, 15 min, 10x, 1.5 ml, 37 C is not in range).
_UNIT = (r"(?!\s?(?:%|°|×|x(?![a-z])|[mnµμu]?(?:l|g|m|mol)(?![a-z])|min(?:ute)?s?(?![a-z])|h(?:rs?|ours?)?(?![a-z])|"
         r"s(?:ec(?:ond)?s?)?(?![a-z])|rpm(?![a-z])|bp(?![a-z])|kb(?![a-z])|k(?![a-z])|cycles?(?![a-z])|"
         r"cells?(?![a-z])))")
# A half-day number no pattern read (the verifier's blocker 1, (ii)): any, in a descriptive field; in protocol and
# processing prose only one that could touch the window under any convention (9.5 to 14.5), because quantities such
# as 0.5, 1.5 or 2.5 are everywhere in protocols. Not after a v (a version number).
_RX_HALF = re.compile(r"(?<![\d.v])\d{1,2}\.5(?!\d)" + _UNIT)
_RX_HALF_PROTOCOL = re.compile(r"(?<![\d.v])(?:9|1[0-4])\.5(?!\d)" + _UNIT)
# A companion integer from 5 to 21, standing alone, in a value or sentence that already gave a reading ("E14, 13
# hearts"); not glued to letters (rep12, 10x) and not a quantity with a unit.
_RX_COMPANION = re.compile(r"(?<![a-z0-9.])(?:[5-9]|1\d|2[01])(?!\d|[.,]\d|[a-z_])" + _UNIT)
# Stage-like tokens that may not be left over (rule 4). On the lower-case residual unless marked case-sensitive.
LEFTOVERS: tuple[tuple[str, re.Pattern, bool], ...] = (
    ("E followed by digits", re.compile(r"(?<![a-z])e[\s._:,-]{0,2}\d"), False),
    ("a capital E followed by digits (glued to letters too)", re.compile(r"E[\s._:,-]{0,2}\d"), True),
    ("a Theiler stage", re.compile(r"theiler|(?<![a-z])ts[\s._-]?\d"), False),
    ("dpc or pc", re.compile(r"dpc|d\.p\.c|(?<![a-z])p\.?c\.?[\s._-]?\d|\dp\.?c(?![a-z])|post[- ]?coit"), False),
    ("a somite count", re.compile(r"somit|\d\s?ss(?![a-z])|(?<![a-z])ss[\s._-]?\d"), False),
    ("a day prefix", re.compile(r"(?<![a-z])(?:[eg]d|e\.d\.?)\s?\d"), False),
    ("a number with a day suffix", re.compile(r"\d\s?d(?![a-z])"), False),
    ("embryonic, stage or day glued to a number", re.compile(r"(?:embryonic|stage|day)\s?\d"), False),
)
LEFTOVER_WORDS = re.compile(r"(?<![a-z])(?:days?|weeks?|wks?|months?)(?![a-z])|gestation|"
                            r"(?<![a-z])d\d{1,2}(?!\d)")
SHIFTABLE = ("embryonic day", "embryonic range", "dpc", "bare number (stage-keyed field)", "bare days (either reading)")
MODES = ("field", "protocol")


def nfkc(text) -> str:
    return unicodedata.normalize("NFKC", str(text))


def _lower(text: str) -> str:
    """Lower case, one character for one (a character whose lower case is longer becomes '?')."""
    return "".join(c.lower() if len(c.lower()) == 1 else "?" for c in text)


def _norm(text: str) -> tuple[list[str], list[str]]:
    """(N, L) as character lists of one length: N the text with accessions, read lengths and scientific numbers
    blanked and spelt stage tokens normalised; L its lower case, with a boundary (a space) made before a capital E or
    TS glued to the character before it."""
    n = _RX_BLANK.sub(lambda m: " " * len(m.group(0)), str(text))
    n = _RX_SPELT.sub(lambda m: f"{m.group(0)[0]}{m.group(1)}.{m.group(2)}", n)
    low = list(_lower(n))
    for m in _RX_GLUED.finditer(n):
        low[m.start() - 1] = " "
    return list(n), low


def readings(text: str, stage_keyed: bool = False, words: bool = True,
             mode: str = "field") -> tuple[list[dict], list[dict]]:
    """(readings, leftovers) of one value (mode 'field') or one protocol sentence (mode 'protocol'), after NFKC.
    Readings are {low, high, kind} intervals as written (no plug-day shift); tokens are consumed in a fixed order
    (Theiler, somites, then in field mode weeks or months, adult, neonatal and postnatal day, then embryonic range, dpc
    range, dpc after the number, embryonic day, days, then, in a stage-keyed field, a bare number as an E-day), so no
    part of the text is read twice. Protocol mode reads no postnatal or adult form: none can meet the window, and in
    protocol prose they describe the parents or reagents. Leftovers are the stage-like tokens no pattern consumed
    (LEFTOVERS; a half-day number; a companion integer when the text gave a reading; with `words`, also the words day,
    week, gestation and month)."""
    if mode not in MODES:
        raise ValueError(mode)
    N, Lc = _norm(nfkc(text))
    out: list[dict] = []
    left: list[dict] = []

    def blank(a: int, b: int) -> None:
        for i in range(a, b):
            N[i] = " "
            Lc[i] = " "

    def take(rx, fn):
        for m in list(rx.finditer("".join(Lc))):
            r = fn(m)
            if isinstance(r, str):
                left.append({"token": m.group(0).strip(), "what": r})
            elif r is not None:
                out.append(r)
            blank(m.start(), m.end())

    def theiler(m):
        ts = int(m.group(1))
        if ts not in THEILER:
            return {"low": 0.0, "high": 1e4, "kind": "theiler stage outside the table"}
        return {"low": THEILER[ts][0], "high": THEILER[ts][1], "kind": f"TS{ts}"}

    def somite_iv(a: int, b: int):
        rows = [r for r in SOMITES if not (b < r[0] or a > r[1])]
        if not rows:
            return {"low": 0.0, "high": 1e4, "kind": "somite count outside the table"}
        return {"low": min(r[2] for r in rows), "high": max(r[3] for r in rows), "kind": f"{a}-{b} somites"}

    full = "".join(Lc)

    def weeks(m):
        n = float(m.group(1))
        n2 = float(m.group(2)) if m.group(2) else n
        if m.group(3).startswith("m"):
            return {"low": BIRTH + 28.0 * n, "high": BIRTH + 31.0 * n2 + 31.0, "kind": "months of age"}
        if not _RX_POSTNATAL_WORD.search(full):
            return "weeks without a postnatal word (gestational or postnatal cannot be told apart)"
        return {"low": BIRTH + 7.0 * n, "high": BIRTH + 7.0 * (n2 + 1), "kind": "weeks after birth"}

    def emb(value: str, kind: str):
        x = float(value)
        return {"low": x, "high": x + 0.5 if "." not in value else x, "kind": kind}

    def e_range(m):
        a, b = emb(m.group(1), ""), emb(m.group(2), "")
        return {"low": min(a["low"], b["low"]), "high": max(a["high"], b["high"]), "kind": "embryonic range"}

    def days(m):
        n = float(m.group(1) or m.group(2))
        return {"low": n, "high": BIRTH + n, "kind": "bare days (either reading)"}

    take(_RX_THEILER, theiler)
    take(_RX_SOMITE_RANGE, lambda m: somite_iv(int(m.group(1)), int(m.group(2))))
    take(_RX_SOMITE, lambda m: somite_iv(*(2 * [int(m.group(1) or m.group(2))])))
    if mode == "field":
        take(_RX_WEEKS, weeks)
        take(_RX_ADULT, lambda m: {"low": BIRTH + 28.0, "high": 1e4, "kind": "adult"})
        take(_RX_NEONATAL, lambda m: {"low": BIRTH, "high": BIRTH + 7.0, "kind": "neonatal"})
        take(_RX_POSTNATAL, lambda m: {"low": BIRTH + float(m.group(1)), "high": BIRTH + float(m.group(1)),
                                       "kind": "postnatal day"})
    take(_RX_E_RANGE, e_range)
    take(_RX_DPC_RANGE, e_range)
    take(_RX_DPC_AFTER, lambda m: emb(m.group(1), "dpc"))
    take(_RX_E_POINT, lambda m: emb(m.group(1), "embryonic day"))
    take(_RX_DAYS, days)
    if stage_keyed:
        take(_RX_BARE, lambda m: emb(m.group(1), "bare number (stage-keyed field)"))
    rl, rn = "".join(Lc), "".join(N)
    for what, rx, case in LEFTOVERS:
        for m in rx.finditer(rn if case else rl):
            left.append({"token": m.group(0).strip(), "what": what})
    for m in (_RX_HALF if mode == "field" else _RX_HALF_PROTOCOL).finditer(rl):
        left.append({"token": m.group(0), "what": "a half-day number no pattern read"})
    if out:
        for m in _RX_COMPANION.finditer(rl):
            left.append({"token": m.group(0), "what": "an integer beside a stage reading (a companion stage?)"})
    if words:
        for m in LEFTOVER_WORDS.finditer(rl):
            left.append({"token": m.group(0), "what": "the word day, week, gestation or month"})
    return out, left


def file_readings(name: str) -> tuple[list[dict], list[dict]]:
    """The stage tokens of a supplementary file name (E8_5, E14.5, HeartE12_5, TS21, E9.5DayEmbryo), read like a field
    after its GSM prefix; every kind counts, and leftovers count."""
    return readings(re.sub(r"^GSM\d+_", " ", name, flags=re.IGNORECASE), stage_keyed=False, words=True)


def touches_window(lo: float, hi: float) -> bool:
    """True if [lo, hi] meets the T1 window (9.5, 13.5]; checked against windows.in_window at both ends."""
    inside = hi > WINDOW[0] and lo <= WINDOW[1]
    if windows.in_window(lo, "T1") or windows.in_window(hi, "T1"):
        assert inside
    return inside


def _sentences(texts: Iterable[str]) -> list[str]:
    return [s.strip() for t in texts for s in re.split(r"(?<=[.;!?])\s+", nfkc(t)) if s.strip()]


_RX_CONV = re.compile(r"(?:(?<![a-z0-9])e|embryonic day|(?<![a-z])day|dpc|d\.p\.c\.?|(?<![a-z])pc)\s?"
                      r"(0\.5|0|1|zero|one)(?!\d|\.\d|[a-z])"
                      r"|(?<![\d.])(0\.5|0|1)\s?(?:dpc|d\.p\.c|pc(?![a-z])|days? post[- ]?coit)"
                      r"|(?<![a-z])(first|1st) day")
_RX_CONV_DAY = re.compile(r"(?<![a-z])(?:days?|dpc|d\.p\.c|e\s?\d|embryonic|gestation(?:al)?|pc\s?\d)|post[- ]?coit")
_CONV_OF = {"0.5": "E0.5", "0": "E0", "zero": "E0", "1": "E1", "one": "E1", "first": "E1", "1st": "E1"}
# A sentence that may state the convention: it names a plug, mating, coitus (not "post coitus", a dpc unit),
# copulation, conception or fertilisation (the verifier's lesser point 2).
CONV_SENTENCE = re.compile(r"plug|(?<![a-z])mat(?:ing|ings|ed)(?![a-z])|copulat|conception|fertili[sz]|"
                           r"(?<!post)(?<!post-)(?<!post )coitus")


def convention(texts: Iterable[str]) -> dict:
    """The plug-day convention stated in protocol text: {name, shift (low, high) or None, sentences, reason}. Each
    sentence that names a plug, mating, coitus, copulation, conception or fertilisation is read: one of E0.5, E0 or E1
    (day 1) names its convention; such a sentence naming a day in any other form, or two different conventions over all
    sentences, gives shift None (UNKNOWN). No such sentence: UNSTATED."""
    found: dict[str, list[str]] = {}
    other: list[str] = []
    for s in _sentences(texts):
        t = _lower(s)
        if not CONV_SENTENCE.search(t):
            continue
        vals = {_CONV_OF[next(g for g in m.groups() if g)] for m in _RX_CONV.finditer(t)}
        if vals:
            for v in vals:
                found.setdefault(v, []).append(s)
        elif _RX_CONV_DAY.search(t):
            other.append(s)
    sents = sorted({x for v in found.values() for x in v} | set(other))
    if other:
        return {"name": "unrecognised", "shift": None, "sentences": sents,
                "reason": f"a protocol sentence names a plug or mating and a day in a form the rule does not read: "
                          f"{other[0]!r}"}
    if len(found) > 1:
        return {"name": "conflict", "shift": None, "sentences": sents,
                "reason": f"the protocol text states more than one plug-day convention ({sorted(found)})"}
    if found:
        name = next(iter(found))
        return {"name": name, "shift": CONVENTIONS[name], "sentences": sents,
                "reason": f"plug = {name} stated: readings widened by {CONVENTIONS[name]}"}
    return {"name": "unstated", "shift": UNSTATED, "sentences": [],
            "reason": f"no plug-day convention stated: readings widened by {UNSTATED} (E0.5, E0 or E1)"}


def protocol_sentences(value: str) -> list[str]:
    """One protocol or processing value as sentences (NFKC), each with its convention tokens (E0.5, day 1, 0.5 dpc ...)
    blanked when the sentence is a convention sentence: the convention is read by convention(); the rest of the sentence
    is read for the stage like any other text."""
    out = []
    for s in _sentences([value]):
        low = _lower(s)
        if CONV_SENTENCE.search(low):
            chars = list(s)
            for m in _RX_CONV.finditer(low):
                for i in range(m.start(), m.end()):
                    chars[i] = " "
            s = "".join(chars)
        out.append(s)
    return out


def _values(name: str, values: Iterable[str]):
    """(text, stage_keyed, words) pieces of one stored field. A characteristics line 'key: value' is read whole; its
    value's bare numbers count when the key names a stage; the key's own words (age, day, ...) are not leftovers."""
    for v in values:
        v = str(v)
        key = None
        if name.startswith("characteristics: "):
            key = name[len("characteristics: "):]
        elif name.startswith("characteristics"):
            k, sep, rest = v.partition(":")
            if sep:
                key, v = k, rest
        if key is None:
            yield [(v, False, True)]
            continue
        keyed = bool(STAGE_KEY.search(key))
        if keyed and _KEY_OTHER_UNIT.search(key) and re.search(r"\d", v) and not re.search(r"[a-z]", v.lower()):
            yield [(key, False, False), (v, False, True), ("unit " + key, None, None)]
            continue
        yield [(key, False, False), (v, keyed, True)]


def stage_reading(fields: Mapping[str, Iterable[str]], files: Iterable[str] = (),
                  protocol: Iterable[str] = (), protocol_fields: Mapping[str, Iterable[str]] | None = None,
                  series_fields: Mapping[str, Iterable[str]] | None = None) -> dict:
    """A sample's stage from its own record (see the rule above). `fields`: descriptive field name -> verbatim values;
    `files`: its supplementary file names; `protocol`: its and its series' text, read for the plug-day convention;
    `protocol_fields`: the sample's own protocol and processing fields, read sentence by sentence (protocol mode); a
    reading there counts as one anywhere else. `series_fields`: the sample's protocol and processing fields that are
    byte-identical across every sample record of its series (the stage-parser exception of ruling row R27, see
    series_wide_fields): they give no stage reading, a pooling word in them refuses (UNKNOWN; see
    pooling_in_series_text), and when there is one the sample's stage must come from its own title and its own
    characteristics, which must both give a reading and agree (else UNKNOWN). The exception can only turn a sample
    that would otherwise be refused into one that is read from its title and characteristics; it never makes an
    unpinned sample usable. Returns {mark: WINDOW | UNKNOWN | PERMITTED, low, high (the widened union), nominal_low,
    nominal_high (as written), readings, leftovers, convention, window_reading, reason, fields_read, files_read,
    series_text_fields}."""
    series_fields = {k: list(v) for k, v in (series_fields or {}).items()}
    pooling = pooling_in_series_text(series_fields)
    conv = convention(protocol)
    sh = conv["shift"] if conv["shift"] is not None else UNSTATED
    rs: list[dict] = []
    left: list[dict] = []
    fields_read, files_read = set(), set()

    def add(where: str, pieces, is_file: bool, mode: str = "field") -> None:
        got: list[dict] = []
        for text, keyed, words in pieces:
            if keyed is None:  # a stage key in a unit other than embryonic days, with a bare number
                left.append({"field": where, "token": text, "what": "a bare number under a key naming another unit"})
                continue
            r, lo = readings(text, keyed, words, mode)
            got += r
            left.extend({"field": where, **x} for x in lo)
            if lo:
                (files_read if is_file else fields_read).add(where)
        if not got:
            return
        (files_read if is_file else fields_read).add(where)
        lo_n, hi_n = min(x["low"] for x in got), max(x["high"] for x in got)
        sa = sh[0] if any(x["kind"] in SHIFTABLE for x in got) else 0.0
        sb = sh[1] if any(x["kind"] in SHIFTABLE for x in got) else 0.0
        kind = got[0]["kind"] if len(got) == 1 else f"span of {len(got)} readings in one value"
        rs.append({"field": where, "kind": kind, "low": lo_n, "high": hi_n, "wlow": lo_n + sa, "whigh": hi_n + sb})

    for k, vals in fields.items():
        for pieces in _values(k, vals):
            add(k, pieces, False)
    for k, vals in (protocol_fields or {}).items():
        for v in vals:
            for sent in protocol_sentences(v):
                add(k, [(sent, False, False)], False, "protocol")
    for f in files:
        add(f"file {f}", [(re.sub(r"^GSM\d+_", " ", f, flags=re.IGNORECASE), False, True)], True)
    window_reading = any(touches_window(r["wlow"], r["whigh"]) for r in rs)
    base = {"readings": rs, "leftovers": left, "convention": conv, "window_reading": window_reading,
            "series_text_fields": sorted(series_fields),
            "fields_read": sorted(fields_read), "files_read": sorted(x[len("file "):] for x in files_read),
            "low": min((r["wlow"] for r in rs), default=None), "high": max((r["whigh"] for r in rs), default=None),
            "nominal_low": min((r["low"] for r in rs), default=None),
            "nominal_high": max((r["high"] for r in rs), default=None)}

    def out(mark: str, reason: str) -> dict:
        return {"mark": mark, "reason": reason, **base}

    if conv["shift"] is None:
        return out("UNKNOWN", conv["reason"])
    if pooling:
        x = pooling[0]
        return out("UNKNOWN", f"a pooling word in a field read as series text ({x['field']}: {x['token']!r})")
    if left:
        x = left[0]
        return out("UNKNOWN", f"an unparsed stage-like token ({x['field']}: {x['token']!r}, {x['what']})")
    if not rs:
        return out("UNKNOWN", "no stage reading in the sample's own record")
    for a in rs:
        for b in rs:
            if a["high"] < b["low"] or b["high"] < a["low"]:
                return out("UNKNOWN", f"conflicting stage readings ({a['field']}: {a['kind']} "
                                      f"[{a['low']:g}, {a['high']:g}]; {b['field']}: {b['kind']} "
                                      f"[{b['low']:g}, {b['high']:g}])")
    if window_reading:
        return out("WINDOW", f"a reading, widened by the plug-day convention {conv['name']} {tuple(sh)}, meets the "
                             f"T1 window (E9.5, E13.5]")
    if series_fields:
        def span(sel):
            got = [r for r in rs if sel(r["field"])]
            return (min(r["low"] for r in got), max(r["high"] for r in got)) if got else None

        t, c = span(lambda f: f == "title"), span(lambda f: f.startswith("characteristics"))
        if t is None or c is None or t != c:
            return out("UNKNOWN", f"fields read as series text ({sorted(series_fields)}): the stage must come from "
                                  f"the sample's own title and characteristics, which must agree (title "
                                  f"{t if t is not None else 'no reading'}, characteristics "
                                  f"{c if c is not None else 'no reading'})")
    return out("PERMITTED", "")


# ------------------------------------------------------------------------------------------------ one sample's record
#
# A record is a GEO family SOFT block as geo.parse_soft gives it: field name (prefix removed) -> values, keys in the
# order they first appear, a key's values in file order. Stored (plan.json, the fetch manifest) as a list of
# [key, [values]] pairs, because the JSON writer sorts object keys.

# Sample fields NOT read for the stage (the verifier's lesser point 3, with status added: it holds the release date).
NOT_STAGE = re.compile(r"supplementary_file.*|contact_.*|platform_id|submission_date|last_update_date|status|"
                       r"data_row_count|_other_lines")
# Sample fields read in protocol mode (sentence by sentence); every other read field is descriptive.
PROTOCOL_FIELD = re.compile(r".*protocol.*|data_processing")
# Sample fields that also feed the plug-day convention (with the series' summary and overall design).
SERIES_PROTOCOL = ("summary", "overall_design")


def record_pairs(rec) -> list[tuple[str, list[str]]]:
    """A record as (key, values) pairs in its order, from a dict or from stored pairs; checked for form."""
    items = rec.items() if isinstance(rec, Mapping) else rec
    out = []
    for kv in items:
        k, v = kv
        if not isinstance(k, str) or not isinstance(v, list | tuple) or not all(isinstance(x, str) for x in v):
            raise Refused("a record must map field names to lists of strings")
        out.append((k, list(v)))
    if len({k for k, _ in out}) != len(out):
        raise Refused("a record holds a field name twice")
    return out


def record_dict(rec) -> dict[str, list[str]]:
    return dict(record_pairs(rec))


def stored_record(rec) -> list[list]:
    """The form a record is stored in (order kept through a sorting JSON writer)."""
    return [[k, v] for k, v in record_pairs(rec)]


def file_name(url: str) -> str:
    return url.rstrip("/").rsplit("/", 1)[-1]


def supplementary(rec) -> list[str]:
    """Supplementary file URLs of a record, in its order (keys supplementary_file*, 'NONE' dropped)."""
    return [v for k, vals in record_pairs(rec) if k.startswith("supplementary_file") for v in vals
            if v.strip().upper() != "NONE"]


def stage_fields(rec) -> dict[str, list[str]]:
    """The descriptive fields a stage is read from: every field of the record except NOT_STAGE and the protocol and
    processing fields (read by protocol_fields)."""
    return {k: v for k, v in record_pairs(rec) if v and not NOT_STAGE.fullmatch(k) and not PROTOCOL_FIELD.fullmatch(k)}


def protocol_fields(rec) -> dict[str, list[str]]:
    """The sample's own protocol and processing fields, read for the stage in protocol mode."""
    return {k: v for k, v in record_pairs(rec) if v and PROTOCOL_FIELD.fullmatch(k) and not NOT_STAGE.fullmatch(k)}


def convention_text(rec, series_protocol: Iterable[str]) -> list[str]:
    """The text read for the plug-day convention: the series' summary and design, then the sample's own protocol,
    processing and description values."""
    d = record_dict(rec)
    own = [v for k, vals in protocol_fields(rec).items() for v in vals] + list(d.get("description", []))
    return [*series_protocol, *own]


# The stage-parser exception (ruling row R27, third cell, in its words): "a protocol or processing field that is
# byte-identical across every sample record of a series is read as series text: it gives no stage reading, and it
# still refuses on pooling words; the sample's stage must then come from its own title and characteristics, which must
# agree". Read here as: a FIELD is one key of the record with all its values (GEO's repeated `!Sample_<key>` lines,
# in order); it is series text when every sample record of the series holds that key with the same list of values,
# compared as the stored text (geo.parse_soft strips surrounding whitespace); a series with one sample has none. A
# series-text field still feeds the plug-day convention (as the series' own summary and design do) and is still printed
# in the processing text; it is never read for a stage, a leftover or a companion integer (pooling_in_series_text
# reads its sentences for readings and leftovers only to decide a refusal, never to give the sample a stage).
SERIES_TEXT_RULE = "t1x-series-text-v1 (ruling row R27: whole field, every sample of the series, at least two)"
# Pooling words: refused anywhere in series text (the exception's own clause): word forms of pool and merge. Words
# such as "combined" and "mixed", which 10x protocol prose also uses for reagents ("master mix", "combined with"), are
# in POOLING_WORDS_STAGED below and refuse in a sentence about stages.
POOLING_WORDS = re.compile(r"(?<![a-z])(?:pool(?:s|ed|ing)?|merg(?:e|es|ed|ing))(?![a-z])")
# Pooling words refused where the value names stages (the reviews of 290da3f and of 75356ec): a value of series text
# with one of these words in any of its sentences refuses when the value AS A WHOLE (any of its sentences, read in
# protocol mode) gives a stage reading or a leftover, or names a stage, an age, a time point or a time course. So
# "Embryos from E8.5 and E10.5 were combined", "Samples from all stages were mixed" and "Embryos were collected at E8.5
# and E10.5; all were then mixed." refuse, while reagent prose in a value that names no stage ("cDNA was combined with
# RT Master Mix.") does not. The word set takes noun forms and glued prefixes (mixture, combination, admixed,
# intermixed, premixed, co-mixed, recombined): the review's pattern, with the Unicode hyphens U+2010 and U+2011 also
# taken wherever it allows a hyphen (refuses more).
_HYPHEN = "[-‐‑]"
POOLING_WORDS_STAGED = re.compile(r"(?<![a-z])(?:pre|inter|ad|co|re)?" + _HYPHEN + r"?(?:combin(?:e|es|ed|ing|ations?)|"
                                  r"mix(?:es|ed|ing|tures?)?|multiplex\w*|superload\w*|aggregat\w*)(?![a-z])")
POOLING_STAGE_WORDS = re.compile(r"(?<![a-z])(?:stages?|ages?|time[ ‐‑-]?points?|"
                                 r"time[ ‐‑-]?courses?)(?![a-z])")


def strip_format_chars(text: str) -> str:
    """The text with every character of Unicode category Cf (zero-width spaces and joiners, soft hyphens ...) removed:
    NFKC keeps them, and one inside a word would defeat a word match."""
    return "".join(c for c in str(text) if unicodedata.category(c) != "Cf")


def pooling_in_series_text(series_fields: Mapping[str, Iterable[str]]) -> list[dict]:
    """The pooling words of the fields read as series text: {field, token}, in field and value order. A word of
    POOLING_WORDS anywhere; a word of POOLING_WORDS_STAGED in any sentence (protocol_sentences) of a value when the
    value AS A WHOLE names stages: any of its sentences gives a stage reading or a leftover in protocol mode, or the
    value matches POOLING_STAGE_WORDS (the review of 75356ec, blocker 3: stages in one sentence and the pooling word in
    the next refuse too). Each value is read both as written and with its Cf characters removed, and a match in either
    counts: removing a character can only add a match that way, never lose one (a Cf character between a letter and the
    word is a boundary as written)."""
    out = []
    for k, vals in series_fields.items():
        for v in vals:
            for text in dict.fromkeys((str(v), strip_format_chars(v))):
                out.extend({"field": k, "token": m.group(0)} for m in POOLING_WORDS.finditer(_lower(nfkc(text))))
                sents = protocol_sentences(text)
                staged = [m.group(0) for sent in sents for m in POOLING_WORDS_STAGED.finditer(_lower(sent))]
                if not staged:
                    continue
                names_stages = bool(POOLING_STAGE_WORDS.search(_lower(nfkc(text)))) or any(
                    any(readings(sent, False, True, "protocol")) for sent in sents)
                if names_stages:
                    out.extend({"field": k, "token": t} for t in staged)
    return out


def series_wide_fields(sample_records: Mapping[str, object]) -> dict[str, list[str]]:
    """The protocol and processing fields (PROTOCOL_FIELD) that are byte-identical across EVERY sample record of one
    series: key -> values. Empty for a series with fewer than two sample records; a field missing from one sample, or
    differing in one value, one character or the order of its values, is not series text."""
    recs = [record_dict(r) for r in sample_records.values()]
    if len(recs) < 2:
        return {}
    out = {}
    for k, v in recs[0].items():
        if not v or not PROTOCOL_FIELD.fullmatch(k) or NOT_STAGE.fullmatch(k):
            continue
        if all(r.get(k) == v for r in recs[1:]):
            out[k] = list(v)
    return out


def series_text_of(rec, series_text: Mapping[str, Iterable[str]] | None) -> dict[str, list[str]]:
    """The sample's own protocol and processing fields that equal the series text given for its series."""
    if not series_text:
        return {}
    own = protocol_fields(rec)
    return {k: v for k, v in own.items() if k in series_text and list(series_text[k]) == v}


def record_stage(rec, series_protocol: Iterable[str] = (), series_text: Mapping[str, Iterable[str]] | None = None
                 ) -> dict:
    """The stage rule on one sample's whole record (the plan, the fetch, the loader and read_reduced all call this).
    `series_text`: the series-wide fields of the sample's series (series_wide_fields over every sample record of the
    series, as the plan computed it); without it no field is read as series text."""
    names = [file_name(u) for u in supplementary(rec)]
    st = series_text_of(rec, series_text)
    own = {k: v for k, v in protocol_fields(rec).items() if k not in st}
    return stage_reading(stage_fields(rec), names, convention_text(rec, list(series_protocol)), own, st)


# Canonical record text (pinned by tests/test_t1_ext_stage.py::test_record_text_canonical_form; README, "Canonical
# record text"):
#   line 1: "t1x-record-v1 SAMPLE <GSM>" or "t1x-record-v1 SERIES <GSE>";
#   then every field of the record in its order (keys in the order they first appear in the SOFT block, a key's values
#     together in file order), one line "<key>[<i>]: <value>" per value (i from 1), or "<key>: (none)" for a key with
#     no value; every field is included: protocol and processing text, relations, every supplementary file URL, the
#     contact fields and the lines the parser kept under _other_lines;
#   under any value holding a character outside printable ASCII (32-126; a tab counts), one more line
#     "<key>[<i>] escaped: <the value with Python's unicode_escape, decoded as ASCII>";
#   lines joined by "\n", one "\n" at the end, UTF-8. A pin's record_sha256 is the sha256 of those bytes.
RECORD_FORMAT = "t1x-record-v1"


def record_text(kind: str, acc: str, rec) -> str:
    if kind not in ("SAMPLE", "SERIES"):
        raise ValueError(kind)
    want = r"GSM\d+" if kind == "SAMPLE" else r"GSE\d+"
    if not re.fullmatch(want, acc):
        raise Refused(f"{acc!r} is not a {kind} accession")
    lines = [f"{RECORD_FORMAT} {kind} {acc}"]
    for k, vals in record_pairs(rec):
        lines += _text_lines(k, vals)
    return "\n".join(lines) + "\n"


def record_sha256(kind: str, acc: str, rec) -> str:
    return text_sha256(record_text(kind, acc, rec))


# ------------------------------------------------------------------------------------------------ the pinned block
#
# atlases.yaml, per atlas (README, "Pinned samples"):
#   pinned:
#     samples:
#       - {accession: GSMn, end: early|late, stage: E8.5, record_sha256: <hex>, read_by: [A, B], date: YYYY-MM-DD}
#     series:
#       - {accession: GSEn, record_sha256: <hex>, read_by: [A, B], date: YYYY-MM-DD}
# A sample is USED only when it is pinned AND the stage rule does not refuse it; nothing is fetched or loaded for an
# atlas whose pinned samples do not cover both ends. fetch.py pin-template prints the block (no hash typed by hand).
SAMPLE_PIN_KEYS = ("accession", "end", "stage", "record_sha256", "read_by", "date")
SERIES_PIN_KEYS = ("accession", "record_sha256", "read_by", "date")
_SHA = re.compile(r"[0-9a-f]{64}")
_PLACEHOLDER = re.compile(r"READER_\d+|PENDING|\?+|", re.IGNORECASE)


def _check_readers(where: str, p: Mapping) -> None:
    rb = p.get("read_by")
    if (not isinstance(rb, list) or len(rb) < 2 or not all(isinstance(x, str) and x.strip() for x in rb)
            or len({x.strip().lower() for x in rb}) != len(rb) or any(_PLACEHOLDER.fullmatch(x.strip()) for x in rb)):
        raise Refused(f"{where}: read_by must name at least two different readers (placeholders refused)")
    if not isinstance(p.get("date"), str) or not re.fullmatch(r"20\d\d-\d\d-\d\d", p["date"]):
        raise Refused(f"{where}: date must be YYYY-MM-DD (quoted as text)")
    if not isinstance(p.get("record_sha256"), str) or not _SHA.fullmatch(p["record_sha256"]):
        raise Refused(f"{where}: record_sha256 must be 64 lower-case hex digits")


def check_pins(acc: str, a: Mapping) -> dict:
    """The pinned block of one atlas, checked for form (its hashes are checked against the records by the plan, the
    fetch and the loader). Returns {samples: {GSM: pin}, series: {GSE: pin}}. An atlas in use or held holds no pins;
    a sample pin's stage is one of its end's labels, and all pins of one end name one stage; series pins name only the
    atlas's series, each once."""
    pb = a.get("pinned")
    if pb in (None, {}):
        return {"samples": {}, "series": {}}
    if not isinstance(pb, Mapping) or set(pb) - {"samples", "series"}:
        raise Refused(f"{acc}: pinned must hold only samples and series")
    if a.get("status") in ("in_use", "held") and (pb.get("samples") or pb.get("series")):
        raise Refused(f"{acc}: an atlas with status {a.get('status')!r} holds no pins")
    samples, series = {}, {}
    named = named_gsms({**a, "accession": acc})
    for p in pb.get("samples") or []:
        if not isinstance(p, Mapping) or sorted(p) != sorted(SAMPLE_PIN_KEYS):
            raise Refused(f"{acc}: sample pin {p!r} must have exactly the keys {SAMPLE_PIN_KEYS}")
        g = p["accession"]
        if not isinstance(g, str) or not re.fullmatch(r"GSM\d+", g):
            raise Refused(f"{acc}: pinned accession {g!r} is not a GSM accession")
        if g in samples:
            raise Refused(f"{acc}: {g} is pinned twice")
        why = never_fetch_refusal(a, gsm=g)
        if why:
            raise Refused(f"{acc}: {why}: it is never pinned")
        if named and named.get(g) != p["end"]:
            raise Refused(f"{acc}: {g} is not a sample the atlas names at the {p['end']} end (the strict form pins "
                          f"only the named samples)")
        if p["end"] not in ("early", "late"):
            raise Refused(f"{acc}: {g}: end must be early or late")
        if p["stage"] not in a[p["end"]]["labels"]:
            raise Refused(f"{acc}: {g}: stage {p['stage']!r} is not a label of the {p['end']} end "
                          f"({a[p['end']]['labels']})")
        _check_readers(f"{acc}: {g}", p)
        samples[g] = dict(p)
    for e in ("early", "late"):
        st = {p["stage"] for p in samples.values() if p["end"] == e}
        if len(st) > 1:
            raise Refused(f"{acc}: the {e} end's pins name two stages {sorted(st)}: pin one")
    for p in pb.get("series") or []:
        if not isinstance(p, Mapping) or sorted(p) != sorted(SERIES_PIN_KEYS):
            raise Refused(f"{acc}: series pin {p!r} must have exactly the keys {SERIES_PIN_KEYS}")
        s = p["accession"]
        if s not in a.get("series", []) or s in series:
            raise Refused(f"{acc}: series pin {s!r} is not a series of the atlas, or is pinned twice")
        _check_readers(f"{acc}: {s}", p)
        series[s] = dict(p)
    return {"samples": samples, "series": series}


def pins_cover_both_ends(pins: Mapping) -> bool:
    return {p["end"] for p in pins["samples"].values()} == {"early", "late"}


# ------------------------------------------------------------------------------------------------ processing text
#
# Positive clearance (the reviews' blocker 4): no pattern decides that an atlas's processing was per sample. The plan
# prints, verbatim, every USED sample's processing, extraction and description text and every listed series' summary,
# overall design and supplementary file names, in the canonical form below; a person reads it and writes its sha256
# into the atlas's processing_cleared entry; plan, fetch and load pass only while the live text hashes to that entry.
# Every line of it also lies inside a pinned record (the series pins and the sample pins); it is kept as its own
# clearance because it puts the processing of all USED samples side by side, which no single record shows.
#
# Canonical form (pinned by tests/test_t1_ext_stage.py::test_processing_text_canonical_form):
#   line 1: "t1x-processing-v2 <ACC>"
#   per series, in the atlas's order: "SERIES <GSE>", then for each key of SERIES_TEXT in that order one line
#     "<key>[<i>]: <value>" per value (i from 1, the record's order), or "<key>: (none)" when the key has no value;
#     then "series_supplementary_file[<i>]: <file name>" per supplementary file of the series record (the verifier's
#     lesser point 4), or "series_supplementary_file: (none)";
#   per USED sample, sorted by accession as text: "SAMPLE <GSM> (<GSE>)", then the keys of SAMPLE_TEXT the same way;
#   under any value holding a character outside printable ASCII (32-126; a tab counts), one more line
#     "<key>[<i>] escaped: <the value with Python's unicode_escape, decoded as ASCII>";
#   lines joined by "\n", one "\n" at the end, encoded as UTF-8. The sha256 is of those bytes.
PROCESSING_FORMAT = "t1x-processing-v2"
SERIES_TEXT = ("summary", "overall_design")
SAMPLE_TEXT = ("data_processing", "extract_protocol_ch1", "description", "growth_protocol_ch1",
               "treatment_protocol_ch1")


def _plain(v: str) -> bool:
    return all(32 <= ord(c) < 127 for c in v)


def _text_lines(key: str, values: Iterable[str]) -> list[str]:
    vals = list(values or [])
    if not vals:
        return [f"{key}: (none)"]
    out = []
    for i, v in enumerate(vals, 1):
        out.append(f"{key}[{i}]: {v}")
        if not _plain(v):
            out.append(f"{key}[{i}] escaped: {v.encode('unicode_escape').decode('ascii')}")
    return out


def processing_text(acc: str, series: Iterable[tuple[str, object]], samples: Mapping[str, Mapping],
                    series_text: Mapping[str, Mapping[str, Iterable[str]]] | None = None) -> str:
    """The canonical processing text. series: (GSE, its record) in the atlas's order; samples: GSM -> {"series": GSE,
    "text": {key: values}} for the USED samples; series_text (given by the plan, the fetch and the load since the
    stage-parser exception of ruling row R27): GSE -> its series-wide fields, printed after the series' file names as
    "series_text <key>[<i>]: <value>" lines, keys sorted as text, or "series_text: (none)", so that the
    readers' clearance covers which fields the exception reads as series text."""
    lines = [f"{PROCESSING_FORMAT} {acc}"]
    for s, rec in series:
        d = record_dict(rec)
        lines.append(f"SERIES {s}")
        for k in SERIES_TEXT:
            lines += _text_lines(k, d.get(k, []))
        lines += _text_lines("series_supplementary_file", [file_name(u) for u in supplementary(rec)])
        if series_text is not None:
            st = series_text.get(s) or {}
            if not st:
                lines.append("series_text: (none)")
            for k in sorted(st):
                lines += _text_lines(f"series_text {k}", st[k])
    for gsm in sorted(samples):
        lines.append(f"SAMPLE {gsm} ({samples[gsm]['series']})")
        for k in SAMPLE_TEXT:
            lines += _text_lines(k, samples[gsm]["text"].get(k, []))
    return "\n".join(lines) + "\n"


def sample_text(rec) -> dict[str, list[str]]:
    d = record_dict(rec)
    return {k: list(d.get(k, [])) for k in SAMPLE_TEXT}


def text_sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def clearance(a: Mapping, sha: str) -> tuple[bool, str]:
    """(passed, why) for an atlas entry's processing_cleared against the sha256 of the text as it stands now."""
    c = a.get("processing_cleared") or {}
    if not c.get("sha256"):
        return False, (f"processing text not cleared: a person reads processing_text.txt (sha256 {sha}) and writes "
                       f"processing_cleared: {{sha256: {sha}, read_by: ..., date: ...}} for {a.get('accession')}")
    if c["sha256"] != sha:
        return False, (f"the processing text now hashes to {sha}, not to processing_cleared {c['sha256']}: it changed "
                       f"or other samples are used; read it again")
    return True, f"cleared by {c.get('read_by', '?')} on {c.get('date', '?')}"


def label_day(label: str) -> float:
    """E-day of an atlas.yaml stage label (E8.5 -> 8.5); windows.parse_stage, fail closed."""
    return windows.parse_stage(label)


def bracket_ok(early: float, late: float) -> None:
    """Both bracket ends outside the T1 window (windows.allowed with tasks T1 only), early before late."""
    for d in (early, late):
        if not windows.allowed(d, TASKS):
            raise Refused(f"stage E{d:g} is inside the T1 window")
    if not early < late:
        raise Refused("the early end must come before the late end")


# ------------------------------------------------------------------------------------------------ per-field sha256
#
# Row T72, cell 3: "The plan prints, for every sample record of the series (window samples included), the sha256 of
# each protocol and processing field and nothing else of a window sample. The readers compare those hashes and open no
# record of a never-fetched sample." A field is one key with all its values; its sha256 is that of its canonical
# lines (_text_lines: '<key>[<i>]: <value>' per value, plus the escaped line under a value holding a character outside
# printable ASCII), joined by newlines with one at the end, UTF-8: the same lines the field has in record_text.


def field_sha256(key: str, values: Iterable[str]) -> str:
    return text_sha256("\n".join(_text_lines(key, values)) + "\n")


def protocol_field_sha256(sample_records: Mapping[str, object]) -> dict[str, dict[str, str]]:
    """GSM -> {field: sha256} for every protocol and processing field (PROTOCOL_FIELD, not NOT_STAGE) of every sample
    record of one series, window samples included; a key no sample holds is absent, a key one sample lacks is
    '(absent)' for it, so the readers see every difference."""
    recs = {g: record_dict(r) for g, r in sample_records.items()}
    keys = sorted({k for d in recs.values() for k, v in d.items()
                   if PROTOCOL_FIELD.fullmatch(k) and not NOT_STAGE.fullmatch(k)})
    return {g: {k: field_sha256(k, d[k]) if k in d else "(absent)" for k in keys} for g, d in sorted(recs.items())}
