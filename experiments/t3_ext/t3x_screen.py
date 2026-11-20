"""The held-out-gene screen of the outside-knockout tables: the SECOND line of defence (stage 3).

The first line is the pinned sample list (datasets.yaml 'pinned', README "The pinned sample list"): two readers read
every sample's verbatim metadata and pin, by accession and by the sha256 of that text, the exact samples that may be
used. This screen still runs on every sample of every screened series at plan time and again on the pinned samples at
fetch and process time. A held-out hit anywhere refuses the whole series; a review term refuses it until it is cleared
for one named sample. It is built to over-refuse: a text blacklist can never be shown complete, so it errs towards
refusing.

It cannot be switched off: the patterns are module constants, no function takes a flag that skips them, and
datasets.yaml can only ADD refusal terms or clear a review term for one named sample or series (never a held-out
term, and never on a sample whose text matches an EXTENDED form).

Held out (docs/DATA_POLICY.md 'Genotype'; the organisers' ruling): beta-catenin (Ctnnb1) and Gata4, any allele, and
any perturbation that phenocopies them.

Every text is normalised before matching, and every variant is screened: the text as given, a latin-1 mojibake
repair ('Î²-catenin' -> 'β-catenin'), URL and HTML unescaping ('%CE%B2', '&beta;'), Unicode NFKC (fullwidth letters,
superscripts, 'ϐ' -> 'β'), zero-width characters and soft hyphens removed, every Unicode dash mapped to '-', a fold
of Cyrillic look-alike letters to Latin ('Gаta4' with a Cyrillic 'а' -> 'Gata4'), then lower case ('ß' is kept).
Stage 4 adds, each as a further variant: every Unicode control (Cc) and format (Cf) character and U+034F deleted
('Ctnn<RLO>b1', 'Ga<LRM>ta4', 'Gata<BS>4', 'Ctn<CGJ>nb1'); the fold extended to Greek and Latin look-alikes (Α Β Τ Ε Ο
Ν Κ Ι, α ο τ β ν ι κ, ɑ and small capitals; so 'Gαta4', 'GΑΤΑ4', 'Gɑta4' and 'CTNNβ1' fold to the gene); every
non-ASCII decimal digit mapped to its ASCII digit ('Gata٤' -> 'Gata4'), and, in one more variant, deleted ('Gata١4'
-> 'Gata4'). A text holding any Cc or Cf character (TAB excepted: it hides nothing) also gets the review term
'control-or-format character', clearable for one named sample like any review term.

HELD-OUT spelling families (each has test strings in tests/test_t3_ext_screen.py):
  Ctnnb1  ctnn + optional separators + b anywhere (ctnnb1, Ctnn-b1, ctnnb2, Ctnnb1tm1Mmt, Ctnnbip1); old symbol Catnb;
          the ENU allele Bfc and the word batface; the Greek short forms β / ß / ꞵ + optional separator + 'cat' with
          NO boundary on either side (CAβcat, βcatf/f, β-cateninfloxed); 'beta' + separator + 'cat' with no boundary
          either side (ibeta-catenin, beta-catGOF), except before catalogue wording; 'bet' + 'catenin'; ASCII 'b' +
          optional separator + 'cat' with no boundary on the right and none on the left, except (a) after 'su' or
          'bo' (subcategory, bobcat) and (b) before catalogue wording ('Cat#', 'Cat. No.', 'cat no', 'catalogue',
          'category'): those two cases are REVIEW, not pass. So bcatGOF, bcatf/f, deltaNbcat, iBcat, shRNAbcat, CAbcat,
          lentiBcat and 'plate B cat' are held out. 'ss-catenin' (ß folded by hand); full names 'catenin beta 1',
          'catenin (cadherin-associated protein), beta 1', 'cadherin associated protein beta'; 'armadillo' (the fly
          orthologue) unless it names a repeat or domain; the phospho-site alleles S33Y, S37A, T41A, S45F, S45del.
  Gata4   gata + optional separator + 4 (Gata4, GATA-4, Gata4;Gata6, Gata4/6), a list of GATA numbers that includes 4
          (GATA5/4, Gata6, 4), a range that includes 4 (Gata1-4, Gata 2-6, Gata 1 to 4), 'four' (Gata-four), Roman
          'GATA IV', the full and shortened names 'GATA binding protein 4', 'GATA-binding factor 4', 'GATA binding 4';
          G4 + any allele or perturbation word, glued or separated (G4 fl/fl, G4f/f, G4fx, G4 KO, G4 -/-, G4 KD, G4
          siRNA, G4 MO, G4 deficient, G4 mutant, G4 morphant); sgG4, shG4, siG4, crG4, moG4. A Gata4 reporter
          knock-in (Gata4-H2B-GFP) is held out too (fails closed).
  Both    the mouse and human Ensembl gene identifiers (ENSMUSG00000006932 Ctnnb1, ENSG00000168036 CTNNB1,
          ENSMUSG00000021944 Gata4, ENSG00000136574 GATA4) and the Entrez Gene IDs (12387, 1499, 14463, 2626) when
          labelled as gene IDs; written from memory, not checked against Ensembl or NCBI. HGNC, MGI and UniProt
          identifiers are left out because they were not certain.
  Drugs   compounds that act on beta-catenin itself or on the TCF complex, plus those stage 1 held: CHIR99021 and
          CHIR98014, CT99021 / CT-99021 (CHIR99021's other name), XAV939, IWR-1, ICG-001, PRI-724, PNU-74654, iCRT*,
          FH535, PKF115-584, CWP232291, tegavivint (BC2059), KYA1797K, MSAB, IWP-*, LGK974 / WNT974, Wnt-C59.
  Guides, hairpins, siRNA, morpholinos, degrons and CRISPRi/a carry the symbol or a short name and are caught by it.

REVIEW terms (refuse until a named human clearance): Gata5, Gata6, Roman GATA V/VI, a bare 'G4' (a G4 ES cell line),
G4 + allele word glued to a letter before it (Atg4 KO: a gene that ends in g4) or with letters glued after a short
allele word ('G4 mouse', 'G4 KOmutant'), a typo-tolerant Gata4
(g + 2 to 4 of a/t + 4: Gat4, Gtaa4) and Ctnnb1 (c + 2 to 4 of t/n + 1 or 2 b: Ctnb1, Cttnb1, Ctnbb1, Cnntb1,
Ctnnnb1), 'GATA binding' with no number; every b / beta + 'cat' that the held-out side leaves out (subcategory,
bobcat, 'Chip B, Cat#', 'Buffer B (Cat. No. 123)', 'LB catalogue', 'TGFbeta (Cat# ...)'); any 'catenin' not caught
above (alpha-catenin, delta-catenin, 'catenin alpha'); Apc (and Apc Min, 'adenomatous polyposis coli'); Axin1/2
(conductin); Gsk3 and 'glycogen synthase kinase'; GSK3 inhibitors (LiCl, lithium, kenpaullone, SB216763, SB415286,
TWS119, AZD1080, LY2090314, BIO / 6-bromoindirubin, CP21R7, IM-12); bare CHIR, Chiron, a bare 'Chi' ('Chi pulse',
'3 uM Chi'), 2i / 3i media ('2i/LIF', '2i medium'); bare IWP, IWR, XAV, LGK, ICG, iCRT and C59; tankyrase, Tnks1/2
and named tankyrase inhibitors (G007-LK, G244-LM, NVP-TNKS656, JW55); porcupine / Porcn, ETC-159, GNF-6231; Wnt
with inhibit / antagonist / agonist / activat / block / suppress / knock / loss / ablat / delet / null / mutant
within 30 characters either way, and Wnt-off / -low / -neg / -deficient; Wnt ligand symbols (Wnt3, Wnt3a ...);
R-spondin, Rspo1-4, Dkk1-4 / Dickkopf, Sfrp with or without a digit; Lrp5 / Lrp6; Tcf7, Tcf7l1, Tcf7l2, Tcf1/3/4,
Lef1, TCF/LEF, also after a dominant-negative or truncation prefix (dnTCF, dnTcf4, dnLef1, dnTCF7L2, ΔNTcf3,
caTcf), 'dominant negative ... Tcf/Lef', the fly dTCF and pangolin; pyrvinium, NSC668036, niclosamide; Wls, Gpr177,
Wntless; the word phenocopy; exon-3 alleles without a symbol (lox(ex3), delta ex3, ex3 flox); a bare deltaN / ΔN;
pannier (the fly Gata4 orthologue); 'arm' + RNAi / dsRNA / KD / mutant / null / KO / -/- (the fly beta-catenin).
Stage 4 (the spelling gaps the reviewers of 9c981e1 listed): spaced or split symbols ('C-t-n-n-b-1', 'Ctn-nb1',
'Cat-nb1', 'G a t a 4', 'Ga-ta-4', 'GATA bp 4', 'G-4 KO', 'CTNNβ1'), short forms ('bkat', 'βc-KO', 'bc-cKO'), catenin
after a digit or LaTeX ('$\\beta$-catenin', 'beta$catenin', 'β1-catenin', 'beta1-catenin'; also EXTENDED, so a bare
'catenin' clearance never passes them), Fzd / Frizzled, Dvl / Dishevelled, Csnk1a1 / CK1alpha, Amer1 / WTX, Bcl9,
Pygo1/2, Rnf43 / Znrf3, BML-284, SKL2001, WAY-262611, QS11, indirubin, Foxy-5, KY02111, LF3, and the README's
phenocopy list Hand2, Mef2c, Zfpm2 / Fog2.
Stage 5 (the reviewers of bf78b1a): the aliases Evi (Wls), Notum, ICAT, Chibby / Cby1, dHand / eHand, Jup /
plakoglobin / gamma-catenin, and the Wnt / beta-catenin targets Tbx6, Sp5 / Sp8, Tbxt / Brachyury / 'T/T', Mixl1 and
Cdx2. Field NAMES are screened like values (screen_record, fetch.screen_all), and a series' contact fields are screened.

EXTENDED forms: a text that matches one (a Ctnnb1 or Gata4 typo form, beta / bet + 'cat' in any context but
catalogue wording, 'catenine') voids every review clearance of that sample (or series): its review hits refuse
whatever datasets.yaml says. A true miss there needs a reviewed code change, not a clearance.

Known over-refusals (fail closed): Bcat1, Bcat2, Ctnnbip1, Ctnnbl1, 'plate B cat', any word ending in b before 'cat'
outside catalogue wording ('lab cat'), 'alpha-catenin, beta-actin' (held: 'catenin, b'), a G4 ES cell line (review),
Atg4 or Gng4 knockouts (review), Tcf3 / Tcf4 meaning E2A / E2-2 (review), 'Chi' as a surname (review), APC as a
fluorochrome (review), Wnt1-Cre (review), ICG staining (review). Known limits: a label that names neither the gene
nor the drug ('Mesp1-Cre; F/F' alone) cannot be caught by any text screen; that is what the pinned list's two readers
are for (README, "What the two readers check").
"""

from __future__ import annotations

import html
import re
import unicodedata
import urllib.parse

import pandas as pd
from t3x_common import Refused

_SEP = r"[\s_\-./:;,()\[\]{}'\"*+&#|~^]*"
# the separators plus LaTeX's '$' and '\' (stage 4: '$\beta$-catenin'), and a beta written in any of its forms
_SEPX = r"[\s_\-./:;,()\[\]{}'\"*+&#|~^$\\]*"
_LATEX_B = r"(?:beta|β|ß|ꞵ)"
_B = r"(?:beta|β|ß|ꞵ|b)"
# catalogue wording after 'cat': 'Cat#', 'Cat. No.', 'cat no', 'Cat nr', 'Cat number', 'catalogue', 'category'
_CATALOGUE = r"\s*(?:\.|#|no(?![a-z])|nr(?![a-z])|number|alog|egor)"
_GLIST = r"(?:[0-9]\s*(?:/|,|&|\+|;|and|or)\s*(?:gata" + _SEP + r")?)*"
_GATA_NAME = rf"gata{_SEP}(?:binding{_SEP}(?:(?:protein|factor){_SEP})?)?"
# allele and perturbation words after G4
_G4W = (rf"(?:fl|flox[a-z]*|lox|fx|f{_SEP}f|c?kos?|knock[a-z]*|ki|null|het|hom|del[a-z]*|mut[a-z]*|kd|sirna|si"
        r"|shrna|sh|sgrna|sg|mo|morph(?:ant|olino)s?|deficien[a-z]*|hypomorph[a-z]*|lof|gof|cre|oe)")
# a dominant-negative or truncation prefix may stand before a Tcf / Lef symbol
_TCF_PRE = r"(?:(?<![a-z])|(?<=dn)|(?<=δn)|(?<=deltan)|(?<=ca))"

HELD_OUT_PATTERNS = [
    # Ctnnb1 / beta-catenin
    rf"ctnn{_SEP}b",
    r"catnb",
    rf"(?:β|ß|ꞵ){_SEP}cat",
    rf"beta{_SEP}cat(?!{_CATALOGUE})",
    rf"bet{_SEP}catenin",
    rf"(?<!su)(?<!bo)b{_SEP}cat(?!{_CATALOGUE})",
    rf"ss{_SEP}catenin",
    rf"î[²2]{_SEP}cat",
    rf"catenin[^a-z0-9]*(?:cadherin{_SEP}associated{_SEP}protein[^a-z0-9]*)?(?:beta|β|ß|ꞵ|b(?![a-z]))",
    rf"cadherin{_SEP}associated{_SEP}protein[^a-z0-9]*(?:beta|β|ß|ꞵ|b(?![a-z]))",
    rf"armadillo(?!{_SEP}(?:repeat|domain|like|family|motif|fold))",
    r"(?<![a-z0-9])bfc(?![a-z0-9])", r"batface",
    r"s33y|s37a|t41a|s45f|s45del",
    # Gata4
    rf"{_GATA_NAME}{_GLIST}4(?![0-9])",
    rf"gata{_SEP}[1-4]\s*(?:-|to)\s*(?:gata{_SEP})?[4-6](?![0-9])",
    rf"{_GATA_NAME}four",
    rf"gata{_SEP}iv",
    rf"(?<![a-z0-9])g4{_SEP}{_G4W}(?![a-z])", r"(?<![a-z0-9])g4\s*[-+]\s*/\s*[-+]",
    r"(?:sg|sh|si|cr|mo)g4(?![0-9])",
    # identifiers (from memory; see the docstring)
    r"ensmusg00000006932|ensg00000168036|ensmusg00000021944|ensg00000136574",
    rf"(?:entrez|gene{_SEP}id|ncbi{_SEP}gene)[^0-9]{{0,12}}(?:12387|1499|14463|2626)(?![0-9])",
    # drugs on beta-catenin or the TCF complex, and those stage 1 held
    rf"chir{_SEP}9[89]0[12][14]", rf"ct{_SEP}99021", rf"xav{_SEP}939", rf"iwr{_SEP}1(?!\d)", rf"icg{_SEP}001",
    rf"pri{_SEP}724", rf"pnu{_SEP}74654", r"icrt", rf"fh{_SEP}535", rf"pkf{_SEP}115", rf"cwp{_SEP}232291",
    r"tegavivint", rf"bc{_SEP}2059", rf"kya{_SEP}1797", r"msab", rf"iwp{_SEP}\d", rf"lgk{_SEP}974",
    rf"wnt{_SEP}974", rf"wnt{_SEP}c{_SEP}59",
]
REVIEW_PATTERNS = [
    # the b / beta + cat forms the held-out side leaves out (catalogue wording, subcategory, bobcat): REVIEW, never
    # pass. The specific forms come first so that the term a clearance names is the whole phrase.
    rf"b{_SEP}cat{_CATALOGUE}", rf"beta{_SEP}cat{_CATALOGUE}", rf"(?:su|bo)b{_SEP}cat",
    rf"bet{_SEP}cat", rf"b{_SEP}cat",
    # typo-tolerant symbols
    r"(?<![a-z])c[tn]{2,4}b{1,2}(?:[0-9]|(?![a-z]))", rf"(?<![a-z])g[at]{{2,4}}{_SEP}4(?![0-9])",
    # Gata neighbours and G4
    rf"gata{_SEP}[56](?![0-9])", rf"gata{_SEP}vi?(?![a-z])", rf"gata{_SEP}binding",
    r"(?<![a-z0-9])g4(?![a-z0-9])", rf"g4{_SEP}{_G4W}(?![a-z])", r"g4\s*[-+]\s*/\s*[-+]",
    rf"g4{_SEP}(?:fl|fx|f{_SEP}f|kd|ko|mo|si|sh|sg|mut|del)",
    # the pathway
    rf"(?<![a-z])apc(?:{_SEP})?(?:min)?(?![a-z])", rf"adenomatous{_SEP}polyposis",
    rf"(?<![a-z])axin{_SEP}[12]?(?![a-z0-9])", r"conductin(?![a-z])",
    rf"(?<![a-z])gsk{_SEP}3", rf"glycogen{_SEP}synthase{_SEP}kinase",
    r"(?<![a-z])licl(?![a-z])", r"lithium", r"kenpaullone", rf"sb{_SEP}216763", rf"sb{_SEP}415286",
    rf"tws{_SEP}119", rf"azd{_SEP}1080", rf"ly{_SEP}2090314", rf"6{_SEP}bromoindirubin", rf"bio{_SEP}acetoxime",
    rf"cp{_SEP}21r7", rf"(?<![a-z])im{_SEP}12(?![0-9])",
    r"(?<![a-z])(?:chir|iwp|iwr|xav|lgk|icg|icrt)(?![a-z])", rf"(?<![a-z0-9])c{_SEP}59(?![0-9])",
    r"(?<![a-z])chiron", r"(?<![a-z])chi(?![a-z])",
    rf"(?<![a-z0-9])[23]i(?:{_SEP}(?:lif|medium|media|condition|culture))?(?![a-z0-9])",
    r"tankyrase", rf"(?<![a-z])tnks{_SEP}[12]?(?![a-z0-9])", rf"g007{_SEP}lk", rf"g244{_SEP}lm", r"tnks656",
    rf"(?<![a-z])jw{_SEP}55(?![0-9])",
    r"porcn|porcupine", rf"etc{_SEP}159", rf"gnf{_SEP}6231",
    r"wnt.{0,30}(?:inhibit|antagon|agonist|activat|block|suppress|knock|loss|ablat|delet|null|mutant)",
    r"(?:inhibit|antagon|agonist|activat|block|suppress|knock|loss|ablat|delet).{0,30}wnt",
    rf"wnt{_SEP}(?:off|low|lo|neg|negative|deficien)",
    rf"(?<![a-z])wnt{_SEP}\d{{1,2}}[ab]?(?![0-9])",
    rf"r{_SEP}spondin", rf"(?<![a-z])rspo{_SEP}[1-4]", rf"(?<![a-z])dkk{_SEP}[1-4]", r"dickkopf",
    r"(?<![a-z])sfrp",
    rf"(?<![a-z])lrp{_SEP}[56](?![0-9])",
    rf"{_TCF_PRE}tcf{_SEP}7(?:{_SEP}l{_SEP}[12])?", rf"{_TCF_PRE}tcf{_SEP}[134](?![0-9])",
    rf"{_TCF_PRE}lef{_SEP}1(?![0-9])", rf"{_TCF_PRE}(?:tcf{_SEP}lef|lef{_SEP}tcf)(?![a-z])",
    rf"(?<![a-z])(?:dn|δn|deltan|ca){_SEP}(?:tcf|lef)",rf"dominant{_SEP}negative.{{0,20}}(?:tcf|lef)",
    r"(?<![a-z])dtcf(?![a-z])", r"pangolin",
    r"pyrvinium", rf"nsc{_SEP}668036", r"niclosamide",
    r"(?<![a-z])wls(?![a-z])", rf"gpr{_SEP}177", r"wntless",
    r"phenocop",
    rf"(?<![a-z])(?:δ|delta|del|d){_SEP}ex(?:on)?{_SEP}3(?![0-9])", rf"lox{_SEP}ex(?:on)?{_SEP}3(?![0-9])",
    rf"(?<![a-z])ex(?:on)?{_SEP}3{_SEP}(?:fl|flox|lox|del)",
    rf"(?<![a-z])(?:δ|delta){_SEP}n(?![a-z])",
    r"pannier",
    rf"(?<![a-z])arm{_SEP}(?:rnai|ds?rna|kd|mutant|null|ko|-/-)",
    # stage 4: the spelling gaps both reviewers of 9c981e1 listed (all REVIEW: over-refuse). Spaced and split symbols:
    rf"(?<![a-z])c{_SEP}t{_SEP}n{_SEP}n{_SEP}b", rf"(?<![a-z])cat{_SEP}nb", rf"ctnn{_SEP}(?:β|ß|ꞵ)",
    rf"(?<![a-z])g{_SEP}a{_SEP}t{_SEP}a{_SEP}4(?![0-9])", rf"gata{_SEP}bp{_SEP}4",
    rf"(?<![a-z0-9])g{_SEP}-{_SEP}4{_SEP}{_G4W}(?![a-z])",
    rf"(?<![a-z])b{_SEP}kat", rf"(?:β|ß|ꞵ){_SEP}c(?![a-z])", rf"(?<![a-z])bc{_SEP}c?kos?(?![a-z])",
    # catenin after a digit or LaTeX ('$\beta$-catenin', 'beta$catenin', 'β1-catenin', 'beta1-catenin'); also EXTENDED
    rf"{_LATEX_B}{_SEPX}1?{_SEPX}catenin", rf"(?:β|ß|ꞵ)1{_SEP}cat",
    # pathway members and compounds
    rf"(?<![a-z])fzd{_SEP}\d*", r"frizzled", rf"(?<![a-z])dvl{_SEP}[1-3]?(?![a-z])", r"dishevelled",
    rf"csnk{_SEP}1{_SEP}a{_SEP}1", rf"(?<![a-z])ck{_SEP}1{_SEP}(?:alpha|α|a)(?![a-z])",
    r"(?<![a-z])(?:amer1|wtx)(?![a-z])", rf"(?<![a-z])bcl{_SEP}9",
    rf"(?<![a-z])pygo{_SEP}[12]?(?![a-z])", rf"(?<![a-z])(?:rnf{_SEP}43|znrf{_SEP}3)",
    rf"bml{_SEP}284", rf"skl{_SEP}2001", rf"way{_SEP}262611", rf"(?<![a-z])qs{_SEP}11(?![0-9])", r"indirubin",
    rf"foxy{_SEP}5(?![0-9])", rf"ky{_SEP}02111", rf"(?<![a-z])lf{_SEP}3(?![0-9])",
    # the README's phenocopy list (Gata4 partners and heart regulators held under the phenocopy clause)
    rf"(?<![a-z])hand{_SEP}2(?![0-9])", rf"(?<![a-z])mef{_SEP}2{_SEP}c", rf"(?<![a-z])zfpm{_SEP}2(?![0-9])",
    rf"(?<![a-z])fog{_SEP}2(?![0-9])",
    # stage 5 (the second reviewer of bf78b1a): aliases of listed terms (Evi = Wls; Notum; ICAT and Chibby, the
    # beta-catenin inhibitors; dHand / eHand, the Hand2 / Hand1 aliases; Jup, plakoglobin, gamma-catenin) ...
    r"(?<![a-z])evi(?![a-z])", r"notum", r"(?<![a-z])icat(?![a-z])", r"chibby|(?<![a-z])cby1?(?![a-z])",
    rf"(?<![a-z])[de]{_SEP}hand(?![a-z0-9])", rf"(?<![a-z])jup(?![a-z])|plakoglobin|(?:gamma|γ){_SEP}catenin",
    # ... and Wnt / beta-catenin targets whose loss phenocopies Ctnnb1 loss in the posterior mesoderm
    rf"(?<![a-z])tbx{_SEP}6(?![0-9])", rf"(?<![a-z])sp{_SEP}[58](?![0-9])",
    rf"(?<![a-z])tbxt(?![a-z])|brachyury|(?<![a-z])t{_SEP}/{_SEP}t(?![a-z])", rf"(?<![a-z])mixl{_SEP}1",
    rf"(?<![a-z])cdx{_SEP}2(?![0-9])",
    r"catenin",
]
# a text matching one of these voids every review clearance of its sample (or series)
EXTENDED_PATTERNS = [
    r"(?<![a-z])c[tn]{2,4}b{1,2}(?:[0-9]|(?![a-z]))", rf"(?<![a-z])g[at]{{2,4}}{_SEP}4(?![0-9])",
    rf"beta?{_SEP}cat(?!{_CATALOGUE})", r"catenine",
    rf"{_LATEX_B}{_SEPX}1?{_SEPX}catenin", rf"(?:β|ß|ꞵ)1{_SEP}cat",
]
HELD_OUT = re.compile("|".join(HELD_OUT_PATTERNS), re.IGNORECASE)
REVIEW = re.compile("|".join(REVIEW_PATTERNS), re.IGNORECASE)
EXTENDED = re.compile("|".join(EXTENDED_PATTERNS), re.IGNORECASE)
# BIO (the GSK3 inhibitor) only in capitals as a word: 'bio' in lower case is too common to be a review term
REVIEW_CASE = re.compile(r"(?<![A-Za-z0-9])BIO(?![A-Za-z0-9])")
# a cell-hashing or demultiplexing value that is not one sample: dropped, never defaulted to control
UNASSIGNED = re.compile(r"^\s*$|negative|doublet|multiplet|unassigned|ambiguous|unknown|^nan$|^none$", re.IGNORECASE)
GUIDE_COLUMN = re.compile(r"sgrna|guide|grna|protospacer|crispr_?id", re.IGNORECASE)

_ZERO_WIDTH = dict.fromkeys(map(ord, "­​‌‍⁠﻿"), None)
_DASHES = dict.fromkeys(map(ord, "‐‑‒–—―−﹘﹣－"), "-")
# Cyrillic, Greek and Latin letters that look like the Latin letters of the two gene names and their short forms
# (stage 4 adds the Greek and Latin ones, and β -> b): folded in one extra variant
_CONFUSABLE = str.maketrans({
    "а": "a", "е": "e", "о": "o", "р": "p", "с": "c", "у": "y", "х": "x", "і": "i", "ј": "j", "ѕ": "s", "ԁ": "d",
    "к": "k", "м": "m", "т": "t", "в": "b", "н": "h", "г": "r", "ɡ": "g",
    "А": "A", "В": "B", "Е": "E", "К": "K", "М": "M", "Н": "H", "О": "O", "Р": "P", "С": "C", "Т": "T", "Х": "X",
    "У": "Y", "І": "I", "Ј": "J", "Ѕ": "S",
    "Α": "A", "Β": "B", "Τ": "T", "Ε": "E", "Ο": "O", "Ν": "N", "Κ": "K", "Ι": "I",
    "α": "a", "ο": "o", "τ": "t", "β": "b", "ν": "v", "ι": "i", "κ": "k",
    "ɑ": "a", "ᴀ": "A", "ʙ": "B", "ɢ": "G", "ᴛ": "T", "ɴ": "N",
})
# a character that is a Unicode control (Cc) or format (Cf) character, or the combining grapheme joiner: deleted in
# an extra variant, and (TAB excepted, which hides nothing) a review term of its own
FORMAT_TERM = "control-or-format character"


def _is_format(ch: str) -> bool:
    return ch == "͏" or unicodedata.category(ch) in ("Cc", "Cf")


def _strip_format(t: str) -> str:
    return "".join(c for c in t if not _is_format(c))


def has_format_char(text) -> bool:
    """True if the text holds a control or format character other than TAB (a review term: FORMAT_TERM)."""
    return any(_is_format(c) and c != "\t" for c in str(text))


def _ascii_digits(t: str, drop: bool = False) -> str:
    """Every non-ASCII character with a decimal value mapped to its ASCII digit (or, with drop, deleted)."""
    out = []
    for c in t:
        if not c.isascii() and unicodedata.decimal(c, None) is not None:
            if not drop:
                out.append(str(unicodedata.decimal(c)))
        else:
            out.append(c)
    return "".join(out)


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
            st = _strip_format(t)
            fold = st.translate(_CONFUSABLE)
            for u in (t, t.translate(_CONFUSABLE), st, fold, _ascii_digits(fold), _ascii_digits(fold, drop=True)):
                if u not in out:
                    out.append(u)
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
    if has_format_char(text):
        out.append(("review", FORMAT_TERM))
    return out


def extended(text) -> list[str]:
    """The EXTENDED forms found in one text (any variant): a non-empty list voids review clearances."""
    out: list[str] = []
    for v in _variants(text):
        for m in EXTENDED.finditer(v.lower()):
            if m.group(0) and m.group(0) not in out:
                out.append(m.group(0))
    return out


def match(text) -> tuple[str, str] | None:
    """('held_out', term) if any held-out term is found, else ('review', term) for the first review term, else
    None."""
    h = hits(text)
    return next((x for x in h if x[0] == "held_out"), h[0] if h else None)


def is_review_term(term: str) -> bool:
    t = _text(term)
    return bool(t == FORMAT_TERM or REVIEW.fullmatch(t) or REVIEW_CASE.fullmatch(t.upper()))


def check_review_cleared(entries: list) -> None:
    """datasets.yaml review_cleared entries: each names one sample (a GSM, or a series accession for a series-level
    field), one REVIEW term and a reason, and may never name a held-out term or an EXTENDED form."""
    if not isinstance(entries, list):
        raise Refused("review_cleared must be a list")
    for e in entries:
        if not isinstance(e, dict) or set(e) != {"sample", "term", "reason"}:
            raise Refused(f"review_cleared entry {e!r} must have exactly sample, term and reason")
        if not str(e["sample"]).strip() or not str(e["reason"]).strip():
            raise Refused(f"review_cleared entry {e!r}: the sample and the reason must be named")
        if HELD_OUT.search(_text(e["term"])) or not is_review_term(e["term"]) or extended(e["term"]):
            raise Refused(f"review_cleared term {e['term']!r} is not a review-only term: a held-out term or an "
                          "extended form is never cleared")


def _cleared(review_cleared, ident: str) -> set[str]:
    check_review_cleared(review_cleared or [])
    return {_text(e["term"]) for e in (review_cleared or []) if str(e["sample"]) == ident}


def screen_text(ident: str, text, extra_refuse: list[str] | None = None, review_cleared: list | None = None,
                void: bool = False) -> list[tuple[str, str]]:
    """The uncleared hits of one text that belongs to `ident` (a sample or series accession). With void (the
    sample's or series' text matches an EXTENDED form) no review term is cleared."""
    cleared = set() if (void or extended(text)) else _cleared(review_cleared, ident)
    if void or extended(text):
        check_review_cleared(review_cleared or [])
    got = [(k, t) for k, t in hits(text) if not (k == "review" and t in cleared)]
    for p in extra_refuse or []:
        rx = re.compile(p, re.IGNORECASE)
        for v in _variants(text):
            m = rx.search(v.lower())
            if m and ("held_out", m.group(0)) not in got:
                got.append(("held_out", m.group(0)))
    return got


def _flat(fields: dict) -> list:
    return [x for v in (fields or {}).values() for x in (v if isinstance(v, list | tuple) else [v])]


def screen_record(sample_id: str, fields: dict, extra_refuse: list[str] | None = None,
                  review_cleared: list | None = None) -> dict:
    """Screen every text field of one sample record (title, characteristics, source, description, protocols, file
    names ...). Returns {'sample', 'verdict': 'pass'|'held_out'|'review', 'hits': [(field, kind, term)],
    'clearance_void': [extended forms]}. A record with no text field at all cannot be screened and is refused. If any
    field matches an EXTENDED form, no review clearance applies to the record. Since stage 5 every field NAME is
    screened like a value too (its hits are reported under '<name> (field name)')."""
    flat = _flat(fields)
    if not any(str(x).strip() for x in flat):
        raise Refused(f"{sample_id}: no text field to screen: refused (unscreenable)")
    void = sorted({t for x in [*flat, *(fields or {})] for t in extended(x)})
    out, verdict = [], "pass"
    for k, v in fields.items():
        vals = v if isinstance(v, list | tuple) else [v]
        for one, label in [(k, f"{k} (field name)"), *((x, k) for x in vals)]:
            for kind, term in screen_text(sample_id, one, extra_refuse, review_cleared, void=bool(void)):
                out.append((label, kind, term))
                if kind == "held_out":
                    verdict = "held_out"
                elif verdict == "pass":
                    verdict = "review"
    return {"sample": sample_id, "verdict": verdict, "hits": out, "clearance_void": void}


def screen_cells(obs: pd.DataFrame, genotype_columns: list[str], hash_columns: list[str] | None = None,
                 guide_map: dict[str, str] | None = None, review_cleared: list | None = None,
                 sample_id: str = "") -> tuple[pd.Series, dict]:
    """Per-cell screen, applied to a cell-label table of one sample before any count row is read. No loader of this
    stage uses it: per-cell genotype designs are excluded by design until a later stage pins them.

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
