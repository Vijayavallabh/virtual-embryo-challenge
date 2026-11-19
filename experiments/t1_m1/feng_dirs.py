"""Feng 2022 directions for M1 (docs/PLAN.md 4.1): the pre-registered readout's primary predictor, per group.

The T1 direction readout (experiments/t1_direction, verdict "build" on 2026-10-28, 3c5cc66) set M1's direction source
to Feng's CD1 hearts, E9.5 -> E14.5 (approved by the captain on 2026-10-21). This module computes exactly the
readout's vectors: per harmonised group, the log ratio of mean CP10k between Feng's E9.5 and E14.5 cells of the
group's Feng labels (pseudocount, minimum cells and labels from experiments/t1_direction/config.yaml), on a T1 gene
order (symbols used more than once in Feng are dropped; unmapped genes are NaN). Within-dataset ratios only: no
external level enters a cell. Feng's statistics are checked against the filter's manifest before use. They come from
the clean reduction of the raw GEO files (experiments/external/feng2022_geo.py, run 2026-11-01 at 6093147), whose
group totals equal the readout's 2026-10-18 statistics to float rounding.

The cell label -> group map is the readout's (config regexes), with the ventricular cardiomyocyte labels inside CM, as
in the readout's population prediction (`cm_includes`).

`bracket` (default today's E9.5 -> E14.5, byte-identical) names the two Feng stages a direction spans; the late-proxy
tries (experiments/harness/TRIES.md T13-T18, the X002 amendment of 2026-11-10) use E14.5 -> E18, the one bracket
legal on every T1 pseudo-task. `direction_halves` and `eb_shrink` serve the same tries.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "experiments" / "t1_direction"))
sys.path.insert(0, str(REPO / "experiments" / "external"))
import dirlib as dl  # noqa: E402
import filterio as fio  # noqa: E402

CONFIG = REPO / "experiments" / "t1_direction" / "config.yaml"
SOURCES = REPO / "experiments" / "external" / "sources.yaml"
BRACKET = ("E9.5", "E14.5")  # the readout's primary predictor, feng_E14.5
BRACKET_DAYS = (9.5, 14.5)
REGISTER = REPO / "EXTERNAL_DATA.md"
X002_AMENDMENT = "E14.5 -> E18"  # the key phrase of X002's late-proxy amendment (PLAN decision log 2026-11-10 (a))
# The cell groups Feng's direction is formed on (build_test.py t1 --feng-groups; off by default). "authors": Feng's own
# labels, every file as built. "none": no Feng direction on any group and no Feng file read, a file with no outside
# data in it (docs/PLAN.md decision log 2026-11-19 (a) to (d): the authors' labels were made with every stage pooled,
# the in-window ones included, and an organiser's condition of 2026-11-16 excludes such annotations).
FENG_GROUPS = ("authors", "none")


def load_config(path: Path = CONFIG) -> dict:
    return yaml.safe_load(Path(path).read_text())


def x002_amended(path: Path | None = None) -> bool:
    """Whether EXTERNAL_DATA.md's X002 row already carries the late-proxy amendment: a use is logged before it enters a
    model (AGENTS.md; experiments/t1_direction/readout.py refuses on its register keys the same way). Only X002's own
    row counts: the phrase elsewhere in the file does not."""
    text = Path(REGISTER if path is None else path).read_text(encoding="utf-8")
    return any(line.startswith("| X002 |") and X002_AMENDMENT in line for line in text.splitlines())


def bracket_fraction(ref_day: float, target_day: float, bracket_days: tuple[float, float] = BRACKET_DAYS) -> float:
    """The share of Feng's bracket a step covers (default E9.5 -> E14.5: 0.4 for E10.5 -> E12.5, 0.2 for E9.5 ->
    E10.5; for E14.5 -> E18, a one-day step is 2/7)."""
    return (target_day - ref_day) / (bracket_days[1] - bracket_days[0])


def label_groups(labels, cfg: dict | None = None) -> dict[str, str]:
    """Cell label -> readout group (config regexes; ventricular CM labels count as CM; unmatched labels 'other:...')."""
    cfg = cfg or load_config()
    m = dl.harmonise(labels, cfg["harmonise"])
    inc = set(cfg.get("cm_includes", []))
    return {lab: ("CM" if g in inc else g) for lab, g in m.items()}


def clean_stats_path() -> Path:
    """cd1_stats.npz of the clean GEO reduction (sources.yaml, feng2022_geo)."""
    return REPO / yaml.safe_load(SOURCES.read_text())["feng2022_geo"]["dst"] / "cd1_stats.npz"


def directions_from_stats(st: dict, genes, cfg: dict,
                          bracket: tuple[str, str] = BRACKET) -> tuple[dict[str, np.ndarray], dict, np.ndarray]:
    """From loaded statistics (`dirlib.load_stats`), group -> log ratio of mean CP10k across `bracket` on `genes` (NaN
    where unmapped), for every group with Feng labels and at least min_cells_ext cells at both ends; the groups' cells
    and the gene map. No file is read, so a late-proxy candidate can call it inside the harness's read guard."""
    idx = dl.sym_index(st["genes"], list(genes))
    eps, mc = cfg["metrics"]["pseudocount_cp10k"], cfg["metrics"]["min_cells_ext"]
    out, cells = {}, {}
    for group, labs in cfg["external_labels"].items():
        fl = labs[0] if labs else None
        if not fl or group in cfg.get("cm_includes", []):  # ventricularCM is inside CM (as in the population design)
            continue
        a, b = dl.feng_sel(st, bracket[0], fl), dl.feng_sel(st, bracket[1], fl)
        if a is None or b is None or min(a["n"], b["n"]) < mc:
            continue
        out[group] = dl.on_t1(dl.lfc_ratio(a, b, eps), idx)
        cells[group] = {bracket[0]: int(a["n"]), bracket[1]: int(b["n"]), "feng_labels": list(fl)}
    return out, cells, idx


def feng_directions(genes, cfg: dict | None = None, stats_path: Path | None = None,
                    bracket: tuple[str, str] = BRACKET,
                    groups: str = FENG_GROUPS[0]) -> tuple[dict[str, np.ndarray], dict]:
    """group -> Feng log ratio across `bracket` (default E9.5 -> E14.5, M1's) on `genes` (NaN where unmapped), for every
    group with Feng labels and at least min_cells_ext cells at both ends; and a record (file, sha256, bracket, groups,
    cells, mapped genes). stats_path: default the clean GEO reduction's statistics. groups: "authors" (the default,
    every file as built) forms the directions on Feng's own labels; "none" forms no direction on any group and reads
    nothing, neither the config nor any Feng file (the record then names no file)."""
    if groups not in FENG_GROUPS:
        raise ValueError(f"groups must be one of {FENG_GROUPS}, got {groups!r}")
    if groups == "none":
        return {}, {"file": None, "sha256": None, "bracket": list(bracket), "groups": {}, "feng_groups": "none"}
    cfg = cfg or load_config()
    path = Path(stats_path) if stats_path is not None else clean_stats_path()
    entry = fio.verify_manifest_output(path)
    st = dl.load_stats(path)
    out, cells, idx = directions_from_stats(st, genes, cfg, bracket)
    eps, mc = cfg["metrics"]["pseudocount_cp10k"], cfg["metrics"]["min_cells_ext"]
    record = {"file": str(path.relative_to(REPO)), "sha256": entry["sha256"], "bracket": list(bracket),
              "pseudocount_cp10k": eps, "min_cells_ext": mc, "groups": cells,
              "mapped_genes": int((idx >= 0).sum()), "of_genes": len(idx)}
    return out, record


def bracket_stats(st: dict, bracket: tuple[str, str]) -> dict:
    """The statistics rows whose stage is one of the bracket's two ends, so that a candidate handed them can reach no
    other Feng stage. A bracket end with no row is refused: nothing else may stand in for it."""
    stages = [k.split("|")[0] for k in st["keys"]]
    missing = sorted(set(bracket) - set(stages))
    if missing:
        raise ValueError(f"Feng's statistics hold no stage {missing}: the bracket {list(bracket)} cannot be computed")
    keep = np.isin(stages, list(bracket))
    rows = len(st["keys"])
    return {k: (v[keep] if k != "genes" and np.shape(v)[:1] == (rows,) else v) for k, v in st.items()}


def direction_halves(st: dict, genes, cfg: dict, bracket: tuple[str, str] = BRACKET,
                     groups=None) -> dict[str, dict[str, np.ndarray]]:
    """Per group with a direction across `bracket` (or among `groups`), on `genes`' order (NaN where unmapped): the log
    ratio from half A alone and from half B alone of the group's Feng cells (the statistics' per-cell halves), and the
    mean CP10k averaged over the two ends (the decile basis of the shrinkage and of the permuted control, N12's). A
    group with a direction but an empty half is refused: its split cannot be read."""
    idx = dl.sym_index(st["genes"], list(genes))
    eps, mc = cfg["metrics"]["pseudocount_cp10k"], cfg["metrics"]["min_cells_ext"]
    out = {}
    for group, labs in cfg["external_labels"].items():
        fl = labs[0] if labs else None
        if not fl or group in cfg.get("cm_includes", []) or (groups is not None and group not in groups):
            continue
        a, b = dl.feng_sel(st, bracket[0], fl), dl.feng_sel(st, bracket[1], fl)
        if a is None or b is None or min(a["n"], b["n"]) < mc:
            continue
        parts = [dl.feng_sel(st, stage, fl, halves=(h,)) for stage in bracket for h in ("A", "B")]
        if any(p is None for p in parts):
            raise ValueError(f"{group}: a half of its Feng cells is empty at a bracket end: no split to read")
        a_a, a_b, b_a, b_b = parts
        out[group] = {"d_a": dl.on_t1(dl.lfc_ratio(a_a, b_a, eps), idx),
                      "d_b": dl.on_t1(dl.lfc_ratio(a_b, b_b, eps), idx),
                      "mean_cp": dl.on_t1((a["cp"] + b["cp"]) / 2.0, idx)}
    return out


def eb_shrink(d, d_a, d_b, basis, bins: int = 10) -> tuple[np.ndarray, list[float]]:
    """Empirical-Bayes reliability shrinkage of a per-gene log ratio toward 0, by decile of expression (`basis`): within
    each decile of the genes with a value, the noise variance of d is the mean of ((d_a - d_b) / 2) ** 2 (if each half
    has noise variance s, d_a - d_b has variance 2 s and d, which pools both halves, about s / 2, a quarter of it), the
    signal variance is the mean of d ** 2 less that noise (floored at 0), and every d in the decile is multiplied by
    the same factor signal / (signal + noise): a gene is not judged by its own halves' disagreement, so genes keep their
    order inside a decile. NaN stays NaN. Returns the shrunk d and each decile's factor, lowest expression first."""
    d, d_a, d_b, basis = (np.asarray(x, dtype=np.float64) for x in (d, d_a, d_b, basis))
    out = d.copy()
    idx = np.flatnonzero(np.isfinite(d) & np.isfinite(d_a) & np.isfinite(d_b) & np.isfinite(basis))
    factors: list[float] = []
    for run in np.array_split(idx[np.argsort(basis[idx], kind="stable")], bins):
        noise = float(np.mean(((d_a[run] - d_b[run]) / 2.0) ** 2)) if run.size else 0.0
        signal = max(float(np.mean(d[run] ** 2)) - noise, 0.0) if run.size else 0.0
        f = signal / (signal + noise) if signal + noise > 0 else 1.0
        out[run] = f * d[run]
        factors.append(f)
    return out, factors
