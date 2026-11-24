"""CHIPQ-C, the partners' part (T3 candidate of ledger row K28; register row X018 as reopened for it): steps 1 to 3 of
the specification's file 02, the structure that the database-only gate G0a tests. Database files only: no fold, no
barred factor, no panel, no score.

Nothing is fetched. The source is the copy chipq_fetch.py pinned under its data root: the analysis list, the
experiment list, the target-gene files and the two screen records, each read through its sha256 in MANIFEST.json. The
archive's co-localisation files are not used (file 01).

  Step 1. A USABLE experiment of an antigen: kept (the committed chipq_kept.yaml for the two queries, the keyword
          filter's placebo_kept.json for every other antigen), not dropped by classes (v) and (vi) of file 03, with a
          column in its antigen's target-gene file and a score above 0 at MIN_BOUND_GENES genes or more. Its class is
          the "Cell type class" field of its row in the pinned experiment list.
  Step 2. L, a query's class: the class with the most usable experiments of the query; a tie goes to the class name
          first in byte order. The query needs QUERY_MIN_USABLE or more there (bar a1).
  Step 3. Candidates: the antigens with CANDIDATE_MIN_USABLE or more usable experiments in L, except the query and
          CANDIDATE_EXCLUDE. At least MIN_CANDIDATES (bar a2).
Classes (v) and (vi) run in code for every antigen and both queries before any count (`screen`): (v) drops a row whose
screened field holds one of CLASS_V_WORDS in any case; (vi) excludes a cell-type class when CLASS_VI_SHARE or more of
its mm10 rows match one of the committed keyword filter's entries CLASS_VI_ENTRIES, and the classes named
CLASS_VI_NAMED.

Every constant of those steps is one named constant below, fixed before any count; none comes from a score and none
changes after a count. The points the specification leaves to a reading are in the block "READINGS".

Only chipq.open_target opens a target-gene file (its every-row mode: every row, kept columns only; never the Average
or the STRING column). This module prints, logs and writes nothing. `g0a` returns counts only. `structure` returns
names (a query's class, the candidate antigens, experiment ids) for the later steps, in memory only: no caller prints,
logs or writes them, and no message here carries an antigen, a class, a gene symbol, a score or an experiment id.

Steps 4 to 13 of file 02 (the strata, the pair value, the picks, the tables, the placebos) are written only after G0a
has passed on the real database. They build on what `structure` returns and change nothing above them.
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import chipq as cq  # noqa: E402  (the kept list and the one reader of a target-gene file)

cf = cq.cf  # chipq_fetch: the manifest, the two lists, the keyword filter and its records
Refused = cq.Refused

TEST_QUERY = cq.TEST_QUERY  # "Ctnnb1"
PRACTICE_QUERY = cq.PRACTICE_QUERY  # "Gata4"
QUERIES = (TEST_QUERY, PRACTICE_QUERY)

# ---------------------------------------------------------------------------------------------------------------------
# The specification's constants (file 02, steps 1 to 3; file 03, classes (v) and (vi); file 04, bars a1 and a2)
MIN_BOUND_GENES = 200  # step 1: a score above 0 at 200 genes or more
QUERY_MIN_USABLE = 2  # step 2, bar a1: the query needs 2 or more usable experiments in its class
CANDIDATE_MIN_USABLE = 3  # step 3: a candidate has 3 or more usable experiments in the class
MIN_CANDIDATES = 11  # step 3, bar a2: at least 11 candidates
CANDIDATE_EXCLUDE = ("Ctnnb1", "Gata4", "Gata6", "Mab21l2")  # step 3: never a candidate; nor the query itself
CLASS_V_WORDS = ("carcinoma", "sarcoma", "lymphoma", "leukemia", "leukaemia", "myeloma", "melanoma", "glioma",
                 "blastoma", "adenoma", "tumor", "tumour", "cancer", "neoplas", "malignan", "metasta", "transform",
                 "immortal")  # class (v): a row drops when a screened field holds one of these, in any case
CLASS_VI_ENTRIES = (15, 17, 18, 19, 20)  # class (vi): these entries, counted from 0, of chipq_fetch.DROP_PATTERNS,
CLASS_VI_OTHER_WORDS = ("stabilised", "activated", "exon-3-deleted", "constitutively active", "stimulated")  # that
#                        is its class (ii) entries other than these five words (class_vi_patterns checks both)
CLASS_VI_SHARE = (5, 100)  # class (vi): a class is excluded when 5% or more of its mm10 rows match
CLASS_VI_NAMED = ("Unclassified", "No description", "Others")  # class (vi): excluded too, if the list has them
CLASS_FIELD = "Cell type class"  # step 1: an experiment's class is this field of its experiment-list row

# READINGS. The literal reading of each point the specification leaves open, fixed here before any count.
#   R1 "its mm10 rows (antigens with a target-gene file)", class (vi): the experiment list's rows of genome mm10 whose
#      antigen the pinned analysis list records with a target-gene file for mm10 (chipq_fetch.analysis_rows, the
#      antigens the fetch's rows step screened; the two queries are among them). Not only the antigens whose file was
#      fetched: the share is "a fact of the database".
#   R2 "match one of entries": as the keyword filter matches, the entry's pattern found in any screened field of the
#      row (chipq_fetch.screened_fields).
#   R3 "5% or more": exact, in integers: 100 * matching rows >= 5 * rows.
#   R4 the number of classes excluded: counted over the classes that have an mm10 row of R1; a class named in
#      CLASS_VI_NAMED counts when it is among them ("if the list has them").
#   R5 "a score above 0 at 200 genes or more": rows of the antigen's file, one row a gene (the reader refuses a file
#      with a gene symbol on two rows, or a row without one).
#   R6 "with a column in its antigen's file": an antigen for which the fetch holds no file has no column, so no usable
#      experiment: one the targets pass recorded as unfetched, and one with fewer than chipq_fetch.PLACEBO_MIN_KEPT
#      kept experiments (the fetch takes no file for it). Any other antigen without a file in the manifest is refused
#      by the reader (an incomplete fetch), as chipq.placebo_pool refuses it.
#   R7 a file is opened only when the antigen's count can matter: a query with a kept experiment left after classes
#      (v) and (vi), and an antigen with CANDIDATE_MIN_USABLE or more such experiments in the query's class. An
#      antigen with fewer cannot reach the bar, whatever its file holds.
#   R8 a query with no usable experiment has no class: both of its bars are then not met ("a bar that cannot be
#      computed is not met") and its two counts are 0.
#   R9 byte order: of the UTF-8 encoding.
_CLASS_COL = cf.EXPERIMENT_COLUMNS.index(CLASS_FIELD)
_V = re.compile("|".join(re.escape(w) for w in CLASS_V_WORDS), re.IGNORECASE)


def _bytes(name: str) -> bytes:
    return name.encode("utf-8")


# ---------------------------------------------------------------------------------------------------------------------
# Classes (v) and (vi) of file 03, in code


def dropped_v(fields) -> bool:
    """Class (v): a screened field holds one of CLASS_V_WORDS, in any case."""
    return any(_V.search(str(f)) for f in fields)


def class_vi_patterns() -> tuple:
    """The compiled patterns of the committed keyword filter's entries CLASS_VI_ENTRIES. Refused unless those entries
    are exactly the filter's class (ii) entries other than CLASS_VI_OTHER_WORDS, the specification's two ways of
    naming them: a changed filter is never read by position."""
    words_ii = [w for c, w, _ in cf.DROP_PATTERNS if c == "ii"]
    by_word = tuple(i for i, (c, w, _) in enumerate(cf.DROP_PATTERNS) if c == "ii" and w not in CLASS_VI_OTHER_WORDS)
    if (by_word != CLASS_VI_ENTRIES or not set(CLASS_VI_OTHER_WORDS) <= set(words_ii)
            or len(words_ii) != len(set(words_ii))
            or [(c, w) for c, w, _ in cf._DROP] != [(c, w) for c, w, _ in cf.DROP_PATTERNS]):
        raise Refused("the committed keyword filter is not the one class (vi) counts its entries in: refused")
    return tuple(cf._DROP[i][2] for i in CLASS_VI_ENTRIES)


def matches_vi(fields, patterns: tuple | None = None) -> bool:
    """A row counts towards its class's share (class (vi)) when one of the entries' patterns is found in any of its
    screened fields."""
    patterns = class_vi_patterns() if patterns is None else patterns
    return any(rx.search(str(f)) for rx in patterns for f in fields)


def share_excludes(matching: int, rows: int) -> bool:
    """CLASS_VI_SHARE or more of the rows match (exact, in integers)."""
    return rows > 0 and CLASS_VI_SHARE[1] * matching >= CLASS_VI_SHARE[0] * rows


class Screen:
    """What the pinned lists and records say, read once and before any count: the classes class (vi) excludes; for
    every antigen with a kept list its kept experiments left after classes (v) and (vi), each with its class; the
    antigens for which the fetch holds no file. It carries the data root, the manifest and the pinned kept list that
    every later read uses, and the usable experiments already counted. Names stay in memory."""

    __slots__ = ("root", "man", "kept", "excluded", "eligible", "no_file", "_usable")

    def __init__(self, root: Path, man: dict, kept, excluded: frozenset, eligible: dict, no_file: frozenset):
        self.root, self.man, self.kept = root, man, kept
        self.excluded, self.eligible, self.no_file = excluded, eligible, no_file
        self._usable: dict[str, dict[str, str]] = {}


def screen(root=None, kept=None, man: dict | None = None) -> Screen:
    """Classes (v) and (vi) for every antigen and both queries, from the two pinned lists, before any count.

    Refused: while the fetch's targets pass is incomplete (as chipq's own steps); when the kept list is not the pinned
    one or has no entry for a query; when a kept experiment has no mm10 row in the experiment list, its row names
    another antigen, an id is kept for two antigens, or an mm10 id appears twice (the pinned files disagree); when the
    keyword filter is not the one class (vi) counts in. No message names an antigen, a class or an experiment."""
    root = cf.root_of(root)
    man = cf.load_manifest(root) if man is None else man
    unfetched = cf.placebo_unfetched(root, man)
    kept = cq.load_kept() if kept is None else kept
    if not isinstance(kept, cq.Kept):
        raise Refused("the kept list must come from load_kept (the pinned chipq_kept.yaml): refused")
    patterns = class_vi_patterns()
    targets, exp_raw = cf._lists(root, man)
    placebo = cf.placebo_kept(root, man)
    ids_of: dict[str, list[str]] = {}
    for q in QUERIES:
        if q not in kept:
            raise Refused("the pinned kept list has no entry for a query: refused")
        ids_of[q] = list(kept[q])
    for a in sorted(placebo, key=_bytes):
        if a in cf.QUERY_ANTIGENS:
            raise Refused("placebo_kept.json holds a query antigen: refused")
        ids_of[a] = list(placebo[a])
    owner: dict[str, str] = {}
    for a, ids in ids_of.items():
        for i in ids:
            if i in owner:
                raise Refused("an experiment is kept twice (for two antigens, or twice for one): refused")
            owner[i] = a
    rows: Counter = Counter()
    hits: Counter = Counter()
    found: dict[str, tuple[str, bool]] = {}
    seen: set[str] = set()
    for f in cf.experiment_rows(exp_raw):
        if f[cf.EXP_GENOME_COL] != cf.GENOME:
            continue
        i = f[cf.EXP_ID_COL]
        if i in seen:
            raise Refused("an mm10 experiment id appears twice in the experiment list: refused")
        seen.add(i)
        a, cls, fields = f[cf.EXP_ANTIGEN_COL], f[_CLASS_COL], cf.screened_fields(f)
        if a in targets:  # reading R1
            rows[cls] += 1
            if matches_vi(fields, patterns):
                hits[cls] += 1
        if i in owner:
            if owner[i] != a:
                raise Refused("a kept experiment's row in the experiment list names another antigen: refused")
            found[i] = (cls, dropped_v(fields))
    missing = len(owner) - len(found)
    if missing:
        raise Refused(f"{missing} kept experiment(s) have no mm10 row in the experiment list: refused")
    excluded = frozenset(c for c in rows if c in CLASS_VI_NAMED or share_excludes(hits[c], rows[c]))
    eligible = {a: [(i, found[i][0]) for i in ids
                    if not found[i][1] and found[i][0] not in excluded and found[i][0] not in CLASS_VI_NAMED]
                for a, ids in ids_of.items()}
    no_file = frozenset(a for a in placebo if a in unfetched or len(placebo[a]) < cf.PLACEBO_MIN_KEPT)  # reading R6
    return Screen(root, man, kept, excluded, eligible, no_file)


# ---------------------------------------------------------------------------------------------------------------------
# Steps 1 to 3


def usable(sc: Screen, antigen: str) -> dict[str, str]:
    """Step 1 for one antigen: {experiment id: its class} of its usable experiments. Its target-gene file is opened
    once, through chipq.open_target's every-row mode, for the columns of its kept experiments left after classes (v)
    and (vi) only; not at all when none is left or the fetch holds no file for it."""
    if antigen in sc._usable:
        return sc._usable[antigen]
    el = sc.eligible.get(antigen)
    if el is None:
        raise Refused("the antigen has no kept list: refused")
    out: dict[str, str] = {}
    if el and antigen not in sc.no_file:
        E, ids, _ = cq.open_target(antigen, None, sc.root, sc.kept if antigen in QUERIES else None,
                                   [i for i, _ in el], sc.man, every_row=True)
        n_bound = (E > 0).sum(axis=0)  # per experiment, the genes with a score above 0
        cls = dict(el)
        out = {i: cls[i] for i, n in zip(ids, n_bound, strict=True) if int(n) >= MIN_BOUND_GENES}
    sc._usable[antigen] = out
    return out


def query_class(us: dict[str, str]) -> tuple[str | None, int]:
    """Step 2: (L, the number of usable experiments in L) from a query's usable experiments {id: class}: the class
    with the most, a tie to the class name first in byte order; (None, 0) when there is none (reading R8)."""
    counts = Counter(us.values())
    if not counts:
        return None, 0
    L = min(counts, key=lambda c: (-counts[c], _bytes(c)))
    return L, counts[L]


def candidates(sc: Screen, q: str, L: str) -> dict[str, list[str]]:
    """Step 3: {antigen: its usable experiments in L, sorted by id} for every antigen with CANDIDATE_MIN_USABLE or
    more of them, except q and CANDIDATE_EXCLUDE; antigens in byte order of their names. An antigen whose kept
    experiments left after classes (v) and (vi) number fewer than that in L is not opened (reading R7)."""
    out: dict[str, list[str]] = {}
    for a in sorted(sc.eligible, key=_bytes):
        if a == q or a in CANDIDATE_EXCLUDE or a in QUERIES:
            continue
        if sum(1 for _, c in sc.eligible[a] if c == L) < CANDIDATE_MIN_USABLE:
            continue
        ids = sorted(i for i, c in usable(sc, a).items() if c == L)
        if len(ids) >= CANDIDATE_MIN_USABLE:
            out[a] = ids
    return out


def structure(sc: Screen, q: str) -> dict:
    """Steps 1 to 3 for one query: {"class": L or None, "query_ids": its usable experiments in L (sorted),
    "candidates": {antigen: its usable experiments in L}}. Names, for the later steps only: never printed, logged or
    written."""
    if q not in QUERIES:
        raise Refused("not a query of this candidate: refused")
    us = usable(sc, q)
    L, n = query_class(us)
    if L is None:
        return {"class": None, "query_ids": [], "candidates": {}}
    ids = sorted(i for i, c in us.items() if c == L)
    if len(ids) != n:
        raise Refused("the query's class count is not its experiments': refused")
    return {"class": L, "query_ids": ids, "candidates": candidates(sc, q, L)}


def bars(st: dict) -> dict:
    """Bars a1 and a2 of gate G0a on one query's structure, as counts and booleans only. A query without a class meets
    neither bar: a bar that cannot be computed is not met (reading R8)."""
    n, m = len(st["query_ids"]), len(st["candidates"])
    has_class = st["class"] is not None
    return {"usable_in_class": n, "candidates": m, "a1": has_class and n >= QUERY_MIN_USABLE,
            "a2": has_class and m >= MIN_CANDIDATES}


def g0a(root=None, kept=None, man: dict | None = None) -> dict:
    """Gate G0a's counts, and nothing else: the number of classes that class (vi) excludes and, per query, its usable
    experiments in its class, its candidates and whether bars a1 and a2 hold. No name of any kind."""
    sc = screen(root, kept, man)
    per = {q: bars(structure(sc, q)) for q in QUERIES}
    return {"classes_excluded": len(sc.excluded), "queries": per,
            "passes": all(v["a1"] and v["a2"] for v in per.values())}


def readings() -> dict:
    """The constants of steps 1 to 3, for a record."""
    return {"MIN_BOUND_GENES": MIN_BOUND_GENES, "QUERY_MIN_USABLE": QUERY_MIN_USABLE,
            "CANDIDATE_MIN_USABLE": CANDIDATE_MIN_USABLE, "MIN_CANDIDATES": MIN_CANDIDATES,
            "CANDIDATE_EXCLUDE": list(CANDIDATE_EXCLUDE), "CLASS_V_WORDS": list(CLASS_V_WORDS),
            "CLASS_VI_ENTRIES": list(CLASS_VI_ENTRIES), "CLASS_VI_SHARE": list(CLASS_VI_SHARE),
            "CLASS_VI_NAMED": list(CLASS_VI_NAMED), "CLASS_FIELD": CLASS_FIELD}
