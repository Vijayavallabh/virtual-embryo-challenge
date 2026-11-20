"""Shared pieces of the outside-knockout factor tables (experiments/t3_ext): the fixed constants, the dataset list's
loader and checks, the stage cap, sha256, and the guard that keeps every output outside the repository.

Every constant below was fixed on 2026-11-20, before any outside count was read, from design numbers only (the plan
of the T3 outside-knockout campaign and the challenger's required changes). None comes from a score.
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
MIN_UMI = 20             # a gene with fewer summed UMIs in any arm of (group, stage) gets log ratio 0
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


def outside_repo(path, repo: Path | None = None) -> Path:
    """The resolved path, refused if it lies inside the repository (this checkout or the main checkout)."""
    p = Path(path).expanduser().resolve()
    roots = {(repo or REPO).resolve()}
    main = (repo or REPO).resolve()
    if ".claude" in main.parts:  # a worktree under <main>/.claude/worktrees/<name>
        roots.add(Path(*main.parts[:main.parts.index(".claude")]))
    for r in roots:
        if p == r or r in p.parents:
            raise Refused(f"{p} lies inside the repository ({r}): outputs and downloads go outside it")
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

_TOP = {"version", "fixed", "notes", "datasets"}
_KEYS = {"accession", "status", "status_reason", "tier", "study", "kind", "design", "control_relation", "chimera",
         "repository", "screen_series", "stages_used", "control_definition", "arms", "exclude_match",
         "design_refuse", "chemistry", "expect", "sample_file_patterns", "files_note", "licence", "why", "risks",
         "verified", "screen_extra_refuse", "review_cleared", "cell_labels", "citation"}


def design_factor(ds: dict) -> float:
    """d of the shrinkage: 1 for a germline or conditional knockout against littermate or same-capture controls, and
    for a chimera with a wild-type-chimera second difference; 0.5 for everything else (F0 CRISPR, heterozygous or
    same-series controls of unverified relation)."""
    if ds["design"] == "chimera":
        return 1.0 if ds.get("chimera", {}).get("second_difference") else 0.5
    if ds["design"] in ("germline", "conditional") and ds["control_relation"] in ("littermate", "same_capture"):
        return 1.0
    return 0.5


def load_datasets(path: Path | None = None) -> dict[str, dict]:
    """datasets.yaml, checked: known keys only (nothing can switch the screen off), known statuses, tiers, designs
    and relations, stages that pass the cap. Returns {accession: record} with 'd' added."""
    path = Path(path or DATASETS)
    doc = yaml.safe_load(path.read_text())
    if not isinstance(doc, dict) or set(doc) - _TOP:
        raise Refused(f"{path}: unknown top-level keys {sorted(set(doc or {}) - _TOP)}")
    out = {}
    for ds in doc.get("datasets") or []:
        acc = ds.get("accession")
        bad = set(ds) - _KEYS
        if bad:
            raise Refused(f"{path}: {acc}: unknown keys {sorted(bad)} (the screen has no switch)")
        if ds.get("status") not in STATUSES:
            raise Refused(f"{path}: {acc}: status {ds.get('status')!r} not in {STATUSES}")
        if ds.get("tier") not in TIERS:
            raise Refused(f"{path}: {acc}: tier {ds.get('tier')!r} not in {TIERS}")
        if ds.get("design") not in DESIGNS or ds.get("control_relation") not in RELATIONS:
            raise Refused(f"{path}: {acc}: design/control_relation not in {DESIGNS}/{RELATIONS}")
        if ds["design"] == "chimera" and ds["status"] == "cleared" and not ds.get("chimera", {}).get(
                "second_difference"):
            raise Refused(f"{path}: {acc}: a chimera without a wild-type-chimera second difference is not cleared")
        for s in ds.get("stages_used") or []:
            stage_ok(s)
        if acc in out:
            raise Refused(f"{path}: {acc} listed twice")
        from t3x_screen import check_review_cleared

        check_review_cleared(ds.get("review_cleared") or [])
        out[acc] = {**ds, "d": design_factor(ds)}
    return out


def dataset(acc: str, path: Path | None = None, need_cleared: bool = True) -> dict:
    """One record; refused if the accession is not listed (or not cleared when need_cleared)."""
    allds = load_datasets(path)
    if acc not in allds:
        raise Refused(f"{acc} is not in datasets.yaml: refused")
    ds = allds[acc]
    if need_cleared and ds["status"] != "cleared":
        raise Refused(f"{acc} has status {ds['status']} ({ds.get('status_reason', '')}): refused")
    return ds


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
