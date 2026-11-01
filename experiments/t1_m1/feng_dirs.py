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


def load_config(path: Path = CONFIG) -> dict:
    return yaml.safe_load(Path(path).read_text())


def bracket_fraction(ref_day: float, target_day: float) -> float:
    """The share of Feng's E9.5 -> E14.5 bracket a step covers: 0.4 for E10.5 -> E12.5, 0.2 for E9.5 -> E10.5."""
    return (target_day - ref_day) / (BRACKET_DAYS[1] - BRACKET_DAYS[0])


def label_groups(labels, cfg: dict | None = None) -> dict[str, str]:
    """Cell label -> readout group (config regexes; ventricular CM labels count as CM; unmatched labels 'other:...')."""
    cfg = cfg or load_config()
    m = dl.harmonise(labels, cfg["harmonise"])
    inc = set(cfg.get("cm_includes", []))
    return {lab: ("CM" if g in inc else g) for lab, g in m.items()}


def clean_stats_path() -> Path:
    """cd1_stats.npz of the clean GEO reduction (sources.yaml, feng2022_geo)."""
    return REPO / yaml.safe_load(SOURCES.read_text())["feng2022_geo"]["dst"] / "cd1_stats.npz"


def feng_directions(genes, cfg: dict | None = None,
                    stats_path: Path | None = None) -> tuple[dict[str, np.ndarray], dict]:
    """group -> Feng E9.5 -> E14.5 log ratio on `genes` (NaN where unmapped), for every group with Feng labels and at
    least min_cells_ext cells at both ends; and a record (file, sha256, groups, cells, mapped genes). stats_path:
    default the clean GEO reduction's statistics."""
    cfg = cfg or load_config()
    path = Path(stats_path) if stats_path is not None else clean_stats_path()
    entry = fio.verify_manifest_output(path)
    st = dl.load_stats(path)
    idx = dl.sym_index(st["genes"], list(genes))
    eps, mc = cfg["metrics"]["pseudocount_cp10k"], cfg["metrics"]["min_cells_ext"]
    out, cells = {}, {}
    for group, labs in cfg["external_labels"].items():
        fl = labs[0] if labs else None
        if not fl or group in cfg.get("cm_includes", []):  # ventricularCM is inside CM (as in the population design)
            continue
        a, b = dl.feng_sel(st, BRACKET[0], fl), dl.feng_sel(st, BRACKET[1], fl)
        if a is None or b is None or min(a["n"], b["n"]) < mc:
            continue
        out[group] = dl.on_t1(dl.lfc_ratio(a, b, eps), idx)
        cells[group] = {BRACKET[0]: int(a["n"]), BRACKET[1]: int(b["n"]), "feng_labels": list(fl)}
    record = {"file": str(path.relative_to(REPO)), "sha256": entry["sha256"], "bracket": list(BRACKET),
              "pseudocount_cp10k": eps, "min_cells_ext": mc, "groups": cells,
              "mapped_genes": int((idx >= 0).sum()), "of_genes": len(idx)}
    return out, record
