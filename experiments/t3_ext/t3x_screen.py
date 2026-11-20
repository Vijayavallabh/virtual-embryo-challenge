"""The held-out-gene screen of the outside-knockout tables. It cannot be switched off: the patterns are module
constants, no function takes a flag that skips them, and datasets.yaml can only ADD refusal terms or clear a
review term for one named sample (never a held-out term).

Held out (docs/DATA_POLICY.md 'Genotype'; the organisers' ruling): beta-catenin (Ctnnb1) and Gata4, any allele, and
any perturbation that phenocopies them. A sample or cell whose genotype, condition, title, source, characteristics,
description, protocol, file name or any other text field matches HELD_OUT is removed before any count is read; a
dataset whose labels cannot be screened is refused. REVIEW terms (Gata6, which the Gata4 validation genotype carries;
Wnt-pathway gain or loss lines and drugs that are not named in HELD_OUT) refuse the dataset until the main session
clears that term for that sample in datasets.yaml with a reason.

Spellings caught (lower-cased text; 'ß' is kept, not folded to 'ss'): ctnnb with or without the 1 (sgCtnnb1, Ctnnb1
fl/fl, Ctnnb1-dTAG, shCtnnb1 ...), catnb, beta / b / β / ß + optional separator + 'cat' with optional 'enin' when
not preceded by a letter (beta-catenin, β-cat, bcat, b-cat fl, ß-catenin), gata4 with an optional separator, and
the Wnt / beta-catenin drugs CHIR99021, XAV939, IWR-1, ICG-001, PNU-74654, iCRT*, IWP-*, LGK974, Wnt-C59. The
pattern over-refuses some unrelated names on purpose (Bcat1, Ctnnbip1, Ctnnbl1): fail closed.
"""

from __future__ import annotations

import re

import pandas as pd
from t3x_common import Refused

_SEP = r"[\s_\-./:]*"
HELD_OUT = re.compile(
    "|".join([
        r"ctnnb",
        r"catnb",
        rf"(?<![a-z])(?:beta|β|ß|b){_SEP}cat(?:enin)?(?![a-z])",
        rf"gata{_SEP}4(?!\d)",
        rf"chir{_SEP}99021", rf"xav{_SEP}939", rf"iwr{_SEP}1(?!\d)", rf"icg{_SEP}001", rf"pnu{_SEP}74654",
        r"(?<![a-z])icrt", rf"(?<![a-z])iwp{_SEP}\d", rf"lgk{_SEP}974", rf"wnt{_SEP}c{_SEP}59",
    ]),
    re.IGNORECASE,
)
REVIEW = re.compile(
    "|".join([
        rf"gata{_SEP}6(?!\d)",
        r"(?<![a-z])apc(?![a-z])",
        rf"(?<![a-z])axin{_SEP}[12]?(?![a-z0-9])",
        rf"(?<![a-z])gsk{_SEP}3",
        r"(?<![a-z])chir(?![a-z])",
        r"catenin",
        r"porcn|porcupine",
        r"wnt.{0,30}(?:agonist|inhibitor|activator|antagonist)",
        r"(?:agonist|inhibitor|activator|antagonist).{0,30}wnt",
    ]),
    re.IGNORECASE,
)
# a cell-hashing or demultiplexing value that is not one sample: dropped, never defaulted to control
UNASSIGNED = re.compile(r"^\s*$|negative|doublet|multiplet|unassigned|ambiguous|unknown|^nan$|^none$", re.IGNORECASE)
GUIDE_COLUMN = re.compile(r"sgrna|guide|grna|protospacer|crispr_?id", re.IGNORECASE)


def _text(v) -> str:
    return str(v).lower()


def match(text) -> tuple[str, str] | None:
    """('held_out', term) or ('review', term) for the first term found in one text, else None."""
    t = _text(text)
    m = HELD_OUT.search(t)
    if m:
        return "held_out", m.group(0)
    m = REVIEW.search(t)
    if m:
        return "review", m.group(0)
    return None


def check_review_cleared(entries: list) -> None:
    """datasets.yaml review_cleared entries: each names one sample, one REVIEW term and a reason, and may never name
    a held-out term."""
    for e in entries:
        if not isinstance(e, dict) or set(e) != {"sample", "term", "reason"}:
            raise Refused(f"review_cleared entry {e!r} must have exactly sample, term and reason")
        if HELD_OUT.search(_text(e["term"])) or not REVIEW.fullmatch(_text(e["term"])):
            raise Refused(f"review_cleared term {e['term']!r} is not a review-only term: a held-out term is never "
                          "cleared")


def screen_record(sample_id: str, fields: dict, extra_refuse: list[str] | None = None,
                  review_cleared: list | None = None) -> dict:
    """Screen every text field of one sample record (title, characteristics, source, description, protocols, file
    names ...). Returns {'sample', 'verdict': 'pass'|'held_out'|'review', 'hits': [(field, term)]}. A record with no
    text field at all cannot be screened and is refused."""
    flat = [x for v in (fields or {}).values() for x in (v if isinstance(v, list | tuple) else [v])]
    if not any(str(x).strip() for x in flat):
        raise Refused(f"{sample_id}: no text field to screen: refused (unscreenable)")
    cleared = {(e["sample"], _text(e["term"])) for e in (review_cleared or [])}
    extra = [re.compile(p, re.IGNORECASE) for p in (extra_refuse or [])]
    hits, verdict = [], "pass"
    for k, v in fields.items():
        vals = v if isinstance(v, list | tuple) else [v]
        for one in vals:
            got = match(one)
            if got is None:
                for p in extra:
                    m = p.search(_text(one))
                    if m:
                        got = ("held_out", m.group(0))
                        break
            if got is None:
                continue
            kind, term = got
            if kind == "review" and (sample_id, _text(term)) in cleared:
                continue
            hits.append((k, kind, term))
            if kind == "held_out":
                verdict = "held_out"
            elif verdict == "pass":
                verdict = "review"
    return {"sample": sample_id, "verdict": verdict, "hits": hits}


def screen_cells(obs: pd.DataFrame, genotype_columns: list[str], hash_columns: list[str] | None = None,
                 guide_map: dict[str, str] | None = None, review_cleared: list | None = None) -> tuple[pd.Series, dict]:
    """Per-cell screen, applied to a cell-label table before any count row is read.

    Refused: a named genotype column missing (no per-cell genotype); a review term in any string column; a guide or
    sgRNA column without a guide-to-gene map, or a guide the map does not hold. Removed (kept False): every cell with
    a held-out term in ANY string column (after mapping guides to genes), a missing or empty genotype, or a hashing
    value that is negative, a doublet or unassigned. Returns (keep mask, record of what was removed and why)."""
    if not genotype_columns:
        raise Refused("no per-cell genotype column named: the cells cannot be screened")
    missing = [c for c in genotype_columns if c not in obs.columns]
    if missing:
        raise Refused(f"per-cell genotype column(s) {missing} absent: the cells cannot be screened")
    df = obs.copy()
    guide_cols = [c for c in df.columns if GUIDE_COLUMN.search(str(c))]
    if guide_cols:
        if not guide_map:
            raise Refused(f"guide column(s) {guide_cols} without a guide-to-gene map: refused")
        for c in guide_cols:
            vals = df[c].astype(str)
            unknown = sorted(set(vals) - set(guide_map))
            if unknown:
                raise Refused(f"guide column {c}: {len(unknown)} guide(s) not in the map: refused")
            df[c + "__gene"] = vals.map(guide_map)
    keep = pd.Series(True, index=df.index)
    why: dict[str, int] = {}
    cleared = {_text(e["term"]) for e in (review_cleared or [])}
    for c in df.columns:
        col = df[c]
        if not (col.dtype == object or isinstance(col.dtype, pd.CategoricalDtype) or pd.api.types.is_string_dtype(col)):
            continue
        for v in pd.unique(col.astype(str)):
            got = match(v)
            if got is None:
                continue
            kind, term = got
            if kind == "review" and _text(term) not in cleared:
                raise Refused(f"cell column {c!r} carries review term {term!r}: refused until cleared")
            if kind == "held_out":
                hit = (col.astype(str) == v) & keep
                why[f"held_out:{c}"] = why.get(f"held_out:{c}", 0) + int(hit.sum())
                keep &= ~hit
    for c in genotype_columns:
        bad = df[c].isna() | df[c].astype(str).str.fullmatch(r"\s*|nan|none|NA", case=False)
        why[f"no_genotype:{c}"] = int((bad & keep).sum())
        keep &= ~bad
    for c in hash_columns or []:
        if c not in df.columns:
            raise Refused(f"hashing column {c!r} absent")
        bad = df[c].isna() | df[c].astype(str).str.contains(UNASSIGNED)
        why[f"unassigned:{c}"] = int((bad & keep).sum())
        keep &= ~bad
    return keep, {"columns_screened": [str(c) for c in df.columns], "removed": why,
                  "kept": int(keep.sum()), "of": int(len(keep))}
