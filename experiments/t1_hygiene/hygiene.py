"""R2's gene modules (experiments/harness/TRIES.md R2, registered at 7aa0a2d; docs/PLAN.md 4.1, R2's provenance rule).

Four technical modules are held flat (f = 1) in Feng's direction for the cells M1 moves along it:
- D, dissociation stress: van den Brink 2017 Supplementary Table 1 (EXTERNAL_DATA.md X007);
- A1, haemoglobin (ambient RNA): Feng's gene symbols matching ^Hb[ab]- (X002; its basis is X009);
- N1, mitochondrially encoded: Feng's gene symbols matching ^mt- (X002);
- N2, strongly nucleus-retained ncRNAs: the Bakken 2018 S6 rows that pass R2's filter (X008).

D and N2 are read from the committed extracts (experiments/external/hygiene_sources.py), each refused unless its sha256
is the one configs/t1_hygiene.yaml pins; A1 and N1 are every Feng gene symbol that matches the config's pattern, each
once, as R2 defines them. (A symbol Feng uses more than once is kept, as R2's classes say. dirlib.sym_index gives it no
Feng value on the T1 genes, so its factor is already 1, but holding it still keeps its counts exactly, outside
apply_factors' keep-depth rescaling. N12 lists such genes among each module's genes without a Feng value.) No symbol is
typed by hand.
Each module is mapped to a given T1 gene-name list by exact, case-sensitive symbol, and the symbols with no T1 name are
listed, never guessed. Van den Brink's Supplementary Table 5 (`D_table5`) is read by N12 in place of D and decides
nothing.

R2's union is the mapped genes of the modules whose line is on (D and N2 always; A1 and N1 are lines the captain can
drop). N12 counts it for R2's gate and writes it to experiments/t1_hygiene/results/n12_hygiene_size/hygiene_genes.csv;
`load_union` gives it to the build only when the file's sha256 equals the config's union_sha256, which the build commit
after N12 pins (null until then, so the build refuses --feng-hygiene).

Importing this module reads no file and changes nothing.
"""

from __future__ import annotations

import hashlib
import re
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "configs" / "t1_hygiene.yaml"
MODULES = ("D", "A1", "N1", "N2")  # R2's modules, in this order everywhere
FROM_EXTRACTS = ("D", "N2")
FROM_PATTERNS = ("A1", "N1")  # on Feng's gene symbols; also the lines the captain can drop
TABLE5 = "D_table5"  # van den Brink's Table 5, read by N12 in place of D (deciding nothing)
UNION_COLUMNS = ("gene", "modules")


class Refused(RuntimeError):
    """An input that is not what configs/t1_hygiene.yaml pins: nothing is built or counted from it."""


def load_config(path: str | Path = CONFIG) -> dict:
    return yaml.safe_load(Path(path).read_text())


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(1 << 22):
            h.update(block)
    return h.hexdigest()


def list_sha256(symbols: Iterable[str]) -> str:
    """The sha256 of a gene list as text, one symbol per line (how N12 records each module's list)."""
    return hashlib.sha256("".join(f"{s}\n" for s in symbols).encode()).hexdigest()


def _rel(path: Path, root: Path) -> str:
    try:
        return Path(path).resolve().relative_to(Path(root).resolve()).as_posix()
    except ValueError:
        return str(path)


# ---------------------------------------------------------------------------------------------------------------------
# Modules


def read_extract(name: str, cfg: Mapping | None = None, root: str | Path = REPO) -> list[str]:
    """The symbols of the committed extract `name` (the config's `extracts`), in the source's order. Refused when the
    file's sha256 is not the pinned one, or when its rows are not the registered number of distinct symbols."""
    cfg = load_config() if cfg is None else cfg
    ent = cfg["extracts"][name]
    path = Path(root) / ent["path"]
    if not path.exists():
        raise Refused(f"{ent['path']}: no such extract (experiments/external/hygiene_sources.py writes it)")
    sha = sha256_file(path)
    if not ent.get("sha256") or sha != ent["sha256"]:
        raise Refused(f"{ent['path']}: sha256 {sha} is not the one configs/t1_hygiene.yaml pins ({ent.get('sha256')})")
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
    symbols = [ln.split("\t")[0] for ln in lines[1:]]  # lines[0] names the columns
    if len(symbols) != int(ent["n"]) or len(set(symbols)) != len(symbols):
        raise Refused(f"{ent['path']}: {len(symbols)} rows ({len(set(symbols))} distinct); the config registers "
                      f"{ent['n']}")
    return symbols


def pattern_module(feng_symbols: Iterable[str], pattern: str) -> list[str]:
    """Every Feng symbol that matches `pattern` (re.search, case-sensitive), each once, in Feng's order (R2: "Feng
    symbols matching ^Hb[ab]-", "Feng symbols beginning mt-"; a repeated symbol is listed once and not dropped)."""
    rx = re.compile(pattern)
    return list(dict.fromkeys(s for s in map(str, feng_symbols) if rx.search(s)))


def build_modules(feng_symbols: Sequence[str], cfg: Mapping | None = None, root: str | Path = REPO,
                  table5: bool = False) -> dict[str, list[str]]:
    """R2's modules as source symbols: D and N2 from their extracts, A1 and N1 from the patterns on Feng's gene
    symbols; with table5, also van den Brink's Table 5 (`D_table5`)."""
    cfg = load_config() if cfg is None else cfg
    mods = {"D": read_extract("D", cfg, root),
            "A1": pattern_module(feng_symbols, cfg["patterns"]["A1"]),
            "N1": pattern_module(feng_symbols, cfg["patterns"]["N1"]),
            "N2": read_extract("N2", cfg, root)}
    if table5:
        mods[TABLE5] = read_extract(TABLE5, cfg, root)
    return mods


def map_symbols(symbols: Iterable[str], t1_genes: Iterable[str]) -> dict[str, list[str]]:
    """Exact, case-sensitive mapping to the T1 gene names: the mapped and the unmapped symbols, each in source order."""
    names = {str(g) for g in t1_genes}
    syms = [str(s) for s in symbols]
    return {"mapped": [s for s in syms if s in names], "unmapped": [s for s in syms if s not in names]}


def map_modules(modules: Mapping[str, Sequence[str]], t1_genes: Iterable[str]) -> dict[str, dict[str, list[str]]]:
    """Each module's source symbols, and those mapped and unmapped to the T1 gene names (`map_symbols`)."""
    names = [str(g) for g in t1_genes]
    return {m: {"source": [str(s) for s in syms], **map_symbols(syms, names)} for m, syms in modules.items()}


def union_modules(cfg: Mapping) -> list[str]:
    """The modules that make R2's union: those whose line is on in the config's `lines`, which must name all four with
    D on. A1 and N1 were registered as lines the captain can drop; his sign-off before N12 ran (docs/PLAN.md decision
    log 2026-11-06 (j)) also turned off N2, since N13 could not confirm that T1 is nuclei."""
    lines = cfg["lines"]
    if set(lines) != set(MODULES) or not lines["D"]:
        raise Refused("configs/t1_hygiene.yaml `lines` must name D, A1, N1 and N2, with D on (TRIES R2)")
    return [m for m in MODULES if bool(lines[m])]


def union(mapped: Mapping[str, Mapping[str, Sequence[str]]], modules: Sequence[str]) -> dict[str, list[str]]:
    """gene -> the modules among `modules` whose mapped genes hold it, sorted by gene."""
    out: dict[str, list[str]] = {}
    for m in modules:
        for g in mapped[m]["mapped"]:
            out.setdefault(g, [])
            if m not in out[g]:
                out[g].append(m)
    return {g: out[g] for g in sorted(out)}


# ---------------------------------------------------------------------------------------------------------------------
# The union file: written by N12, read by the build


def union_text(members: Mapping[str, Sequence[str]], header: Sequence[str]) -> str:
    """hygiene_genes.csv: '#' lines (what it is and its provenance), then gene,modules (modules joined by ';')."""
    bad = [g for g in members if not g or any(ch in g for ch in ",;#\r\n") or g != g.strip()]
    if bad:
        raise Refused(f"gene names that the union file cannot hold: {bad[:5]}")
    lines = [f"# {h}" for h in header] + [",".join(UNION_COLUMNS)]
    lines += [f"{g},{';'.join(ms)}" for g, ms in members.items()]
    return "\n".join(lines) + "\n"


def load_union(path: str | Path | None = None, cfg: Mapping | None = None, root: str | Path = REPO) -> dict:
    """The union N12 wrote, for the build: refused unless the file's sha256 equals the config's union_sha256 (null
    until the build commit after N12 pins it). Returns the file (relative to the repository), its sha256, the genes in
    the file's order, the number of them in each module and the sources' sha256 (docs/PLAN.md 4.1, R2's provenance
    rule: each T1 build records the sources' and the union's sha256)."""
    cfg = load_config() if cfg is None else cfg
    path = Path(root) / cfg["union"]["path"] if path is None else Path(path)
    rel = _rel(path, Path(root))
    pinned = cfg.get("union_sha256")
    if not pinned:
        raise Refused("configs/t1_hygiene.yaml pins no union_sha256 yet: it is pinned at the build commit after N12 "
                      "(R2's commit order)")
    if not path.exists():
        raise Refused(f"{rel}: no union file (N12 writes it)")
    sha = sha256_file(path)
    if sha != pinned:
        raise Refused(f"{rel}: sha256 {sha} is not the union_sha256 configs/t1_hygiene.yaml pins ({pinned})")
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
    if not lines or lines[0] != ",".join(UNION_COLUMNS):
        raise Refused(f"{rel}: the first row is not {','.join(UNION_COLUMNS)}")
    genes, counts = [], Counter()
    for ln in lines[1:]:
        gene, mods = ln.split(",")
        genes.append(gene)
        counts.update(mods.split(";"))
    if len(set(genes)) != len(genes) or set(counts) - set(MODULES):
        raise Refused(f"{rel}: a repeated gene or an unknown module")
    return {"file": rel, "sha256": sha, "genes": genes, "modules": {m: counts[m] for m in MODULES if counts[m]},
            "sources": {k: v["sha256"] for k, v in cfg["sources"].items()}}
