"""The held-out-gene screen of the outside-knockout tables. It cannot be switched off: the patterns are module
constants, no function takes a flag that skips them, and datasets.yaml can only ADD refusal terms or clear a
review term for one named sample or series (never a held-out term).

Held out (docs/DATA_POLICY.md 'Genotype'; the organisers' ruling): beta-catenin (Ctnnb1) and Gata4, any allele, and
any perturbation that phenocopies them. A sample or cell whose genotype, condition, title, source, characteristics,
description, protocol, file name or any other text field matches HELD_OUT is removed before any count is read; a
dataset whose labels cannot be screened is refused. REVIEW terms (pathway-level stand-ins and neighbours that are not
the held-out genes themselves) refuse the dataset until the main session clears that term for that sample (or that
series) in datasets.yaml with a reason. Over-refusal is by design: fail closed.

Every text is normalised before matching, and every variant is screened: the text as given, a latin-1 mojibake
repair ('Î²-catenin' -> 'β-catenin'), URL and HTML unescaping ('%CE%B2', '&beta;'), Unicode NFKC (fullwidth letters,
superscripts, 'ϐ' -> 'β'), zero-width characters and soft hyphens removed, every Unicode dash mapped to '-', then
lower case ('ß' is kept, not folded to 'ss').

Spelling families caught as HELD-OUT (each has a test string in tests/test_t3_ext_screen.py):
  Ctnnb1  symbol in any case with or without separators (ctnnb1, CTNNB1, Ctnn-b1, ctnnb, zebrafish ctnnb2), old symbol
          Catnb, the ENU allele name Bfc (batface); full names 'catenin beta 1', 'catenin (cadherin-associated
          protein), beta 1', 'cadherin associated protein beta'; beta / b / β / ß / ꞵ + optional separator + 'cat' or
          'catenin' (beta-catenin, β-cat, bcat, b-cat, ß-catenin, B_catenin, betacatenin, 'ss-catenin'), not preceded
          by a letter unless the letters are a perturbation prefix (sg, g, sh, si, cr, crispr, crispri, crispra, mo,
          anti, ko, dn: sgbcat, shBcat, sgβcat) and not followed by a letter unless it is an allele suffix (bcatKO,
          bcatfl, bcatcKO, βcatMO); 'armadillo' (the fly orthologue) unless it names a repeat or domain; the
          phospho-site alleles S33Y, S37A, T41A, S45F, S45del (stabilised forms); allele and line names carry the
          symbol (Ctnnb1 flox, fl/fl, lox(ex3), Ctnnb1tm1Mmt, Ctnnb1-dTAG) and are caught by it.
  Gata4   gata + optional separator + 4 (Gata4, GATA-4, gata 4, Gata4;Gata6, Gata4/6), a list of GATA numbers that
          includes 4 (GATA5/4, Gata6, 4), Roman 'GATA IV', full names 'GATA binding protein 4' and 'GATA-binding factor
          4'; G4 glued to an allele word (G4 fl/fl, G4 KO, G4 -/-). A Gata4 reporter knock-in (Gata4-H2B-GFP and the
          like) is HELD-OUT too: many reporter knock-ins disrupt one copy, and the label alone cannot tell, so it fails
          closed.
  Both    the mouse and human Ensembl gene identifiers (ENSMUSG00000006932 Ctnnb1, ENSG00000168036 CTNNB1,
          ENSMUSG00000021944 Gata4, ENSG00000136574 GATA4) and the Entrez Gene IDs (12387, 1499, 14463, 2626) when
          labelled as gene IDs; written from memory, not checked against Ensembl or NCBI in this session. HGNC, MGI
          and UniProt identifiers are left out because they were not certain.
  Drugs   compounds that act on beta-catenin itself or on the TCF complex, plus those stage 1 already held:
          CHIR99021 and CHIR98014, XAV939, IWR-1, ICG-001, PRI-724, PNU-74654, iCRT*, FH535, PKF115-584, CWP232291,
          tegavivint (BC2059), KYA1797K, MSAB, IWP-*, LGK974 / WNT974, Wnt-C59.
  Guides, hairpins, siRNA, morpholinos, degrons and CRISPRi/a carry the symbol or the short name and are caught by it
  (sgCtnnb1, shGata4, Ctnnb1 siRNA, MO-ctnnb1, Gata4-AID, sgbcat).

REVIEW terms (refuse until a named human clearance): Gata5, Gata6, Roman GATA V/VI, a bare 'G4'; Apc (and Apc Min,
'adenomatous polyposis coli'); Axin1/2 (conductin); Gsk3 and 'glycogen synthase kinase'; GSK3 inhibitors (LiCl,
lithium, kenpaullone, SB216763, SB415286, TWS119, AZD1080, LY2090314, BIO / 6-bromoindirubin); bare CHIR, IWP, IWR,
XAV, LGK, ICG, iCRT and C59; tankyrase, Tnks1/2 and named tankyrase inhibitors (G007-LK, G244-LM, NVP-TNKS656, JW55);
porcupine / Porcn and ETC-159, GNF-6231; Wnt with inhibit / antagonist / agonist / activat / block / suppress /
knock / loss / ablat / delet / null / mutant within 30 characters either way; Wnt ligand symbols (Wnt3, Wnt3a ...);
R-spondin, Rspo1-4, Dkk1-4 / Dickkopf, Sfrp1-5; Lrp5 / Lrp6; Tcf7, Tcf7l1, Tcf7l2, Tcf1/3/4, Lef1, TCF/LEF; Wls,
Gpr177, Wntless; the word phenocopy; exon-3 alleles without a symbol (lox(ex3), delta ex3, ex3 flox); a bare
deltaN / ΔN; pannier (the fly Gata4 orthologue); any 'catenin' not caught above (alpha-catenin).

Known over-refusals (fail closed): Bcat1, Ctnnbip1, Ctnnbl1, 'armadillo protein', 'alpha-catenin, beta-actin', a
G4 ES cell line (review), Tcf3 / Tcf4 meaning E2A / E2-2 (review). Known limits: a label that names neither the gene
nor the drug ('Mesp1-Cre; F/F' alone) cannot be caught by any text screen; such a dataset is caught only if a
reviewer reads its plan (README, 'Running on the DGX').
"""

from __future__ import annotations

import html
import re
import unicodedata
import urllib.parse

import pandas as pd
from t3x_common import Refused

_SEP = r"[\s_\-./:;,()\[\]{}'\"*+&#|~^]*"
_B = r"(?:beta|β|ß|ꞵ|b)"
# a beta/b short name may follow a non-letter or one of these perturbation prefixes (fixed-width lookbehinds)
_PRE = ("(?:(?<![a-z])|(?<=g)|(?<=sh)|(?<=si)|(?<=cr)|(?<=crispr)|(?<=crispri)|(?<=crispra)|(?<=mo)|(?<=anti)"
        "|(?<=ko)|(?<=dn))")
# ... and may be followed by a non-letter or an allele / perturbation suffix
_SUF = r"(?:(?![a-z])|(?=(?:c?ko|fl|flox|lox|kd|null|del|mut|mo|sg|sh|si|aid|dtag|deg|het|hom|ex|oe|s)(?![a-z])))"
_GLIST = r"(?:[0-9]\s*(?:/|,|&|\+|;|and|or)\s*(?:gata" + _SEP + r")?)*"

HELD_OUT_PATTERNS = [
    # Ctnnb1 / beta-catenin
    rf"ctnn{_SEP}b",
    r"catnb",
    rf"{_PRE}{_B}{_SEP}cat(?:enin)?{_SUF}",
    rf"(?<![a-z])ss{_SEP}catenin",
    rf"î[²2]{_SEP}cat",
    rf"catenin[^a-z0-9]*(?:cadherin{_SEP}associated{_SEP}protein[^a-z0-9]*)?{_B}(?![a-z])",
    rf"cadherin{_SEP}associated{_SEP}protein[^a-z0-9]*{_B}(?![a-z])",
    rf"(?<![a-z])armadillo(?!{_SEP}(?:repeat|domain|like|family|motif|fold))",
    r"(?<![a-z0-9])bfc(?![a-z0-9])",
    r"(?<![a-z0-9])(?:s33y|s37a|t41a|s45f|s45del)(?![a-z0-9])",
    # Gata4
    rf"gata{_SEP}(?:binding{_SEP}(?:protein|factor){_SEP})?{_GLIST}4(?![0-9])",
    rf"gata{_SEP}iv(?![a-z])",
    rf"(?<![a-z0-9])g4(?:{_SEP}(?:fl|flox|lox|c?ko|null|het|del|mut|knock[a-z]*)(?![a-z])|\s*[-+]\s*/\s*[-+])",
    # identifiers (from memory; see the docstring)
    r"ensmusg00000006932|ensg00000168036|ensmusg00000021944|ensg00000136574",
    rf"(?:entrez|gene{_SEP}id|ncbi{_SEP}gene)[^0-9]{{0,12}}(?:12387|1499|14463|2626)(?![0-9])",
    # drugs on beta-catenin or the TCF complex, and those stage 1 held
    rf"chir{_SEP}9[89]0[12][14]", rf"xav{_SEP}939", rf"iwr{_SEP}1(?!\d)", rf"icg{_SEP}001", rf"pri{_SEP}724",
    rf"pnu{_SEP}74654", r"(?<![a-z])icrt", rf"(?<![a-z])fh{_SEP}535", rf"pkf{_SEP}115", rf"cwp{_SEP}232291",
    r"tegavivint", rf"(?<![a-z])bc{_SEP}2059", rf"kya{_SEP}1797", r"(?<![a-z])msab(?![a-z])",
    rf"(?<![a-z])iwp{_SEP}\d", rf"lgk{_SEP}974", rf"wnt{_SEP}974", rf"wnt{_SEP}c{_SEP}59",
]
REVIEW_PATTERNS = [
    rf"gata{_SEP}[56](?![0-9])", rf"gata{_SEP}vi?(?![a-z])", r"(?<![a-z0-9])g4(?![a-z0-9])",
    rf"(?<![a-z])apc(?:{_SEP})?(?:min)?(?![a-z])", rf"adenomatous{_SEP}polyposis",
    rf"(?<![a-z])axin{_SEP}[12]?(?![a-z0-9])", r"conductin(?![a-z])",
    rf"(?<![a-z])gsk{_SEP}3", rf"glycogen{_SEP}synthase{_SEP}kinase",
    r"(?<![a-z])licl(?![a-z])", r"lithium", r"kenpaullone", rf"sb{_SEP}216763", rf"sb{_SEP}415286",
    rf"tws{_SEP}119", rf"azd{_SEP}1080", rf"ly{_SEP}2090314", rf"6{_SEP}bromoindirubin", rf"bio{_SEP}acetoxime",
    r"(?<![a-z])(?:chir|iwp|iwr|xav|lgk|icg|icrt)(?![a-z])", rf"(?<![a-z0-9])c{_SEP}59(?![0-9])",
    r"tankyrase", rf"(?<![a-z])tnks{_SEP}[12]?(?![a-z0-9])", rf"g007{_SEP}lk", rf"g244{_SEP}lm", r"tnks656",
    rf"(?<![a-z])jw{_SEP}55(?![0-9])",
    r"porcn|porcupine", rf"etc{_SEP}159", rf"gnf{_SEP}6231",
    r"wnt.{0,30}(?:inhibit|antagon|agonist|activat|block|suppress|knock|loss|ablat|delet|null|mutant)",
    r"(?:inhibit|antagon|agonist|activat|block|suppress|knock|loss|ablat|delet).{0,30}wnt",
    rf"(?<![a-z])wnt{_SEP}\d{{1,2}}[ab]?(?![0-9])",
    rf"r{_SEP}spondin", rf"(?<![a-z])rspo{_SEP}[1-4]", rf"(?<![a-z])dkk{_SEP}[1-4]", r"dickkopf",
    rf"(?<![a-z])sfrp{_SEP}[1-5]",
    rf"(?<![a-z])lrp{_SEP}[56](?![0-9])",
    rf"(?<![a-z])tcf{_SEP}7(?:{_SEP}l{_SEP}[12])?", rf"(?<![a-z])tcf{_SEP}[134](?![0-9])",
    rf"(?<![a-z])lef{_SEP}1(?![0-9])", rf"(?<![a-z])(?:tcf{_SEP}lef|lef{_SEP}tcf)(?![a-z])",
    r"(?<![a-z])wls(?![a-z])", rf"gpr{_SEP}177", r"wntless",
    r"phenocop",
    rf"(?<![a-z])(?:δ|delta|del|d){_SEP}ex(?:on)?{_SEP}3(?![0-9])", rf"lox{_SEP}ex(?:on)?{_SEP}3(?![0-9])",
    rf"(?<![a-z])ex(?:on)?{_SEP}3{_SEP}(?:fl|flox|lox|del)",
    rf"(?<![a-z])(?:δ|delta){_SEP}n(?![a-z])",
    r"pannier",
    r"catenin",
]
HELD_OUT = re.compile("|".join(HELD_OUT_PATTERNS), re.IGNORECASE)
REVIEW = re.compile("|".join(REVIEW_PATTERNS), re.IGNORECASE)
# BIO (the GSK3 inhibitor) only in capitals as a word: 'bio' in lower case is too common to be a review term
REVIEW_CASE = re.compile(r"(?<![A-Za-z0-9])BIO(?![A-Za-z0-9])")
# a cell-hashing or demultiplexing value that is not one sample: dropped, never defaulted to control
UNASSIGNED = re.compile(r"^\s*$|negative|doublet|multiplet|unassigned|ambiguous|unknown|^nan$|^none$", re.IGNORECASE)
GUIDE_COLUMN = re.compile(r"sgrna|guide|grna|protospacer|crispr_?id", re.IGNORECASE)

_ZERO_WIDTH = dict.fromkeys(map(ord, "­​‌‍⁠﻿"), None)
_DASHES = dict.fromkeys(map(ord, "‐‑‒–—―−﹘﹣－"), "-")


def _variants(v) -> list[str]:
    """Every normalised form of one text, not lower-cased (see the module docstring)."""
    s = str(v)
    raw = [s]
    try:
        fixed = s.encode("latin-1").decode("utf-8")
        if fixed != s:
            raw.append(fixed)
    except (UnicodeEncodeError, UnicodeDecodeError):
        pass
    out = []
    for r in raw:
        for t in (r, html.unescape(urllib.parse.unquote(r))):
            t = unicodedata.normalize("NFKC", t).translate(_ZERO_WIDTH).translate(_DASHES)
            if t not in out:
                out.append(t)
    return out


def _text(v) -> str:
    return unicodedata.normalize("NFKC", str(v)).translate(_ZERO_WIDTH).translate(_DASHES).lower()


def hits(text) -> list[tuple[str, str]]:
    """Every (kind, term) found in one text over all its normalised variants: held-out terms first, then review
    terms. Terms are lower case."""
    out: list[tuple[str, str]] = []
    vs = _variants(text)
    for kind, rx in (("held_out", HELD_OUT), ("review", REVIEW)):
        for v in vs:
            for m in rx.finditer(v.lower()):
                if m.group(0) and (kind, m.group(0)) not in out:
                    out.append((kind, m.group(0)))
    for v in vs:
        for m in REVIEW_CASE.finditer(v):
            if ("review", m.group(0).lower()) not in out:
                out.append(("review", m.group(0).lower()))
    return out


def match(text) -> tuple[str, str] | None:
    """('held_out', term) if any held-out term is found, else ('review', term) for the first review term, else
    None."""
    h = hits(text)
    return next((x for x in h if x[0] == "held_out"), h[0] if h else None)


def is_review_term(term: str) -> bool:
    t = _text(term)
    return bool(REVIEW.fullmatch(t) or REVIEW_CASE.fullmatch(t.upper()))


def check_review_cleared(entries: list) -> None:
    """datasets.yaml review_cleared entries: each names one sample (a GSM, or a series accession for a series-level
    field), one REVIEW term and a reason, and may never name a held-out term."""
    for e in entries:
        if not isinstance(e, dict) or set(e) != {"sample", "term", "reason"}:
            raise Refused(f"review_cleared entry {e!r} must have exactly sample, term and reason")
        if not str(e["sample"]).strip() or not str(e["reason"]).strip():
            raise Refused(f"review_cleared entry {e!r}: the sample and the reason must be named")
        if HELD_OUT.search(_text(e["term"])) or not is_review_term(e["term"]):
            raise Refused(f"review_cleared term {e['term']!r} is not a review-only term: a held-out term is never "
                          "cleared")


def _cleared(review_cleared, ident: str) -> set[str]:
    check_review_cleared(review_cleared or [])
    return {_text(e["term"]) for e in (review_cleared or []) if str(e["sample"]) == ident}


def screen_text(ident: str, text, extra_refuse: list[str] | None = None, review_cleared: list | None = None
                ) -> list[tuple[str, str]]:
    """The uncleared hits of one text that belongs to `ident` (a sample or series accession)."""
    cleared = _cleared(review_cleared, ident)
    got = [(k, t) for k, t in hits(text) if not (k == "review" and t in cleared)]
    for p in extra_refuse or []:
        rx = re.compile(p, re.IGNORECASE)
        for v in _variants(text):
            m = rx.search(v.lower())
            if m and ("held_out", m.group(0)) not in got:
                got.append(("held_out", m.group(0)))
    return got


def screen_record(sample_id: str, fields: dict, extra_refuse: list[str] | None = None,
                  review_cleared: list | None = None) -> dict:
    """Screen every text field of one sample record (title, characteristics, source, description, protocols, file
    names ...). Returns {'sample', 'verdict': 'pass'|'held_out'|'review', 'hits': [(field, kind, term)]}. A record
    with no text field at all cannot be screened and is refused."""
    flat = [x for v in (fields or {}).values() for x in (v if isinstance(v, list | tuple) else [v])]
    if not any(str(x).strip() for x in flat):
        raise Refused(f"{sample_id}: no text field to screen: refused (unscreenable)")
    out, verdict = [], "pass"
    for k, v in fields.items():
        for one in (v if isinstance(v, list | tuple) else [v]):
            for kind, term in screen_text(sample_id, one, extra_refuse, review_cleared):
                out.append((k, kind, term))
                if kind == "held_out":
                    verdict = "held_out"
                elif verdict == "pass":
                    verdict = "review"
    return {"sample": sample_id, "verdict": verdict, "hits": out}


def screen_cells(obs: pd.DataFrame, genotype_columns: list[str], hash_columns: list[str] | None = None,
                 guide_map: dict[str, str] | None = None, review_cleared: list | None = None,
                 sample_id: str = "") -> tuple[pd.Series, dict]:
    """Per-cell screen, applied to a cell-label table of one sample before any count row is read.

    Refused: a named genotype column missing (no per-cell genotype); a held-out or review term in a column NAME; a
    review term in any string column not cleared for `sample_id`; a guide or sgRNA column without a guide-to-gene
    map, or a guide the map does not hold. Removed (kept False): every cell with a held-out term in ANY string column
    (after mapping guides to genes), a missing or empty genotype, or a hashing value that is negative, a doublet or
    unassigned. Returns (keep mask, record of what was removed and why)."""
    if not genotype_columns:
        raise Refused("no per-cell genotype column named: the cells cannot be screened")
    missing = [c for c in genotype_columns if c not in obs.columns]
    if missing:
        raise Refused(f"per-cell genotype column(s) {missing} absent: the cells cannot be screened")
    for c in obs.columns:
        got = screen_text(sample_id, c, review_cleared=review_cleared)
        if got:
            raise Refused(f"cell column name {c!r} carries {got[0][0]} term {got[0][1]!r}: refused")
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
    for c in df.columns:
        col = df[c]
        if not (col.dtype == object or isinstance(col.dtype, pd.CategoricalDtype) or pd.api.types.is_string_dtype(col)):
            continue
        for v in pd.unique(col.astype(str)):
            got = screen_text(sample_id, v, review_cleared=review_cleared)
            if not got:
                continue
            if any(k == "held_out" for k, _ in got):
                hit = (col.astype(str) == v) & keep
                why[f"held_out:{c}"] = why.get(f"held_out:{c}", 0) + int(hit.sum())
                keep &= ~hit
            else:
                raise Refused(f"cell column {c!r} carries review term {got[0][1]!r}: refused until cleared")
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
