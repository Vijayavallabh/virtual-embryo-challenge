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

`groups` (default "authors", byte-identical) names the cell groups. "none" forms no direction at all and reads no Feng
file (main's option, unchanged). "recomputed" forms them from the permitted-only relabelling of
experiments/external/feng_relabel.py (feng_relabel.yaml's group map and fallback; the organiser's condition of
2026-11-16 on annotations learned using protected cells, register row X002). That path reads the recomputed statistics
file alone, refuses it unless its manifest says it was written by feng_relabel through its gate (frozen yaml as it
stands, no unsourced gene, committed code, a clean tree, row N54 and the ruling), passed the yaml's unassigned cap and
is committed at HEAD as the labeller's RECORD byte for byte, and refuses any key whose label is not a recomputed one.
Its record carries the yaml's terminal route.
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
# labels, every file as built. "none": no Feng direction on any group and no Feng file read (docs/PLAN.md decision
# log 2026-11-19 (a) to (e): the authors' labels were made with every stage pooled, the in-window ones included, and
# an organiser's condition of 2026-11-16 excludes such annotations). "recomputed": the groups
# experiments/external/feng_relabel.py recomputes from permitted cells alone (register rows X002 and X010), off unless
# asked for; only the captain's ruling puts it in an upload. "authors" stays first (the default).
FENG_GROUPS = ("authors", "none", "recomputed")
RELABEL_CONFIG = REPO / "experiments" / "external" / "feng_relabel.yaml"
RELABEL_SOURCE = "feng_relabel"  # the `source` of the recomputed file's manifest
RELABEL_PREFIX = "r_"  # every recomputed label starts with it (feng_relabel.LABEL_PREFIX); no author label does


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
    nothing, neither the config nor any Feng file (the record then names no file); "recomputed"
    (`recomputed_directions`) forms them on the permitted-only groups of experiments/external/feng_relabel.py."""
    if groups not in FENG_GROUPS:
        raise ValueError(f"groups must be one of {FENG_GROUPS}, got {groups!r}")
    if groups == "none":
        return {}, {"file": None, "sha256": None, "groups": {}, "feng_groups": "none"}  # no step, so no bracket
    if groups == "recomputed":
        return recomputed_directions(genes, cfg, stats_path, bracket)
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


# ------------------------------------------------------------------ the recomputed groups (off by default)


def load_relabel_config(path: Path | None = None) -> dict:
    """experiments/external/feng_relabel.yaml, structurally checked by the labeller's own check."""
    import feng_relabel

    return feng_relabel.load_config(RELABEL_CONFIG if path is None else path)


def recomputed_stats_path(rcfg: dict | None = None) -> Path:
    """The recomputed statistics file the yaml names (output.dir / output.stats)."""
    rcfg = rcfg or load_relabel_config()
    return REPO / rcfg["output"]["dir"] / rcfg["output"]["stats"]


def groups_config(cfg: dict, rcfg: dict) -> dict:
    """The readout's config with its Feng label lists (each group's first external list) replaced by the recomputed
    labels of feng_relabel.yaml's group map; a group the map does not name (Blood, Hepatocyte, ventricularCM) gets no
    Feng labels and so no direction. The JAX lists, the harmonisation and the metrics are untouched; the file
    experiments/t1_direction/config.yaml itself is never edited."""
    out = dict(cfg)
    out["external_labels"] = {g: [list(rcfg["groups"][g]) if g in rcfg["groups"] else None, *labs[1:]]
                              for g, labs in cfg["external_labels"].items()}
    return out


def record_matches(man_path: Path) -> bool:
    """The labeller's manifest at `man_path` is committed at HEAD as its RECORD, byte for byte (so the run is on the
    record before any reader takes its file)."""
    import feng_relabel

    head = feng_relabel.committed_at_head(feng_relabel.RECORD)
    return head is not None and head == Path(man_path).read_bytes()


def check_recomputed(path: Path, relabel_config: Path | None = None, require_usable: bool = True,
                     require_record: bool = True) -> tuple[dict, dict]:
    """The recomputed statistics file's checks, all before it is read: listed with a matching sha256 in its manifest
    (filterio.verify_manifest_output), written by feng_relabel through its gate (the yaml frozen, no gene unsourced,
    the labeller's code committed, a clean tree, row N54 and a ruling recorded), marked usable (the yaml's per-stage
    unassigned cap passed; `require_usable` False only for N54's listing, which forms no direction from an unusable
    file), written from the yaml as it stands now (its sha256), and its manifest committed at HEAD as the labeller's
    RECORD (`require_record`; final review, lesser point 1). Returns the manifest entry and the manifest."""
    entry = fio.verify_manifest_output(path)
    man_path = next(d / "MANIFEST.json" for d in Path(path).resolve().parents if (d / "MANIFEST.json").exists())
    man = yaml.safe_load(man_path.read_text())  # JSON is YAML
    if man.get("source") != RELABEL_SOURCE:
        raise ValueError(f"{path} was not written by feng_relabel (manifest source {man.get('source')!r}): the "
                         "recomputed option reads only the permitted-only relabelling")
    gate = man.get("gate") or {}
    if man.get("config", {}).get("status") != "frozen":
        raise ValueError(f"{path} was written from a yaml that was not frozen")
    if man.get("unsourced") or man.get("allow_unsourced"):
        raise ValueError(f"{path} was written with unsourced genes")
    if man.get("provenance", {}).get("git_uncommitted"):
        raise ValueError(f"{path} was written by uncommitted labeller code")
    if not (gate.get("clean") is True and gate.get("row") == "N54" and gate.get("ruling")):
        raise ValueError(f"{path} was written without the labeller's gate (clean tree, row N54, the ruling)")
    if require_usable and entry.get("usable") is not True:
        raise ValueError(f"{path} is marked unusable by its manifest (the per-stage unassigned cap failed: "
                         f"{man.get('cap')}); the terminal rule's route is none")
    want = fio.sha256_file(RELABEL_CONFIG if relabel_config is None else relabel_config)
    if man.get("config", {}).get("sha256") != want:
        raise ValueError(f"{path} was written from another feng_relabel.yaml (sha256 "
                         f"{man.get('config', {}).get('sha256')}, now {want}): relabel again under a new row")
    if require_record and not record_matches(man_path):
        raise ValueError(f"{path}: the relabelling's manifest is not committed at HEAD as the labeller's RECORD (byte "
                         "for byte): commit it first, so the run is on the record before anything reads its file")
    return entry, man


def check_recomputed_keys(st: dict) -> None:
    """Every key's label part is a recomputed label: no author label can come back through this file."""
    bad = sorted({k.split("|")[2] for k in st["keys"] if not k.split("|")[2].startswith(RELABEL_PREFIX)})
    if bad:
        raise ValueError(f"the recomputed statistics hold labels that are not recomputed: {bad[:5]}")


def apply_fallback(dirs: dict[str, np.ndarray], rcfg: dict) -> tuple[dict[str, np.ndarray], dict]:
    """feng_relabel.yaml's group-level fallback: a group of its map without its own direction (fewer than min_cells_ext
    cells at an end) takes its fallback group's own direction, or none ('none', or a fallback without its own
    direction). One level: only own directions are handed on. Returns the directions and, per group, what it took."""
    own = set(dirs)
    out = dict(dirs)
    record = {}
    for g in rcfg["groups"]:
        if g in own:
            record[g] = {"own": True, "fired": False, "takes": g}
            continue
        fb = rcfg["fallback"].get(g, "none")
        takes = fb if fb != "none" and fb in own else None
        if takes is not None:
            out[g] = np.array(dirs[takes], copy=True)
        record[g] = {"own": False, "fired": True, "takes": takes}
    return out, record


def recomputed_directions(genes, cfg: dict | None = None, stats_path: Path | None = None,
                          bracket: tuple[str, str] = BRACKET, relabel_config: Path | None = None
                          ) -> tuple[dict[str, np.ndarray], dict]:
    """As feng_directions, on the recomputed groups: the recomputed file (default the yaml's) after check_recomputed
    and check_recomputed_keys, the readout's pseudocount and minimum (the yaml's minimum must equal it), the yaml's
    group map in place of the authors' labels, then the yaml's fallback. The returned directions hold only the groups
    with a direction after the fallback (a group with none is absent, so no cell of it moves). The record names the
    file, its sha256, the yaml's sha256 and the route, and says which groups fell back."""
    cfg = cfg or load_config()
    rpath = RELABEL_CONFIG if relabel_config is None else Path(relabel_config)
    rcfg = load_relabel_config(rpath)
    if int(rcfg["min_cells_ext"]) != int(cfg["metrics"]["min_cells_ext"]):
        raise ValueError("feng_relabel.yaml's min_cells_ext differs from the readout's")
    path = Path(stats_path) if stats_path is not None else recomputed_stats_path(rcfg)
    entry, man = check_recomputed(path, rpath)
    st = dl.load_stats(path)
    check_recomputed_keys(st)
    out, cells, idx = directions_from_stats(st, genes, groups_config(cfg, rcfg), bracket)
    out, fired = apply_fallback(out, rcfg)
    import feng_relabel

    route = feng_relabel.route(entry.get("usable") is True, fired, rcfg)
    eps, mc = cfg["metrics"]["pseudocount_cp10k"], cfg["metrics"]["min_cells_ext"]
    try:
        shown = str(path.resolve().relative_to(REPO))
    except ValueError:
        shown = str(path)
    record = {"file": shown, "sha256": entry["sha256"], "bracket": list(bracket), "pseudocount_cp10k": eps,
              "min_cells_ext": mc, "groups": cells, "mapped_genes": int((idx >= 0).sum()), "of_genes": len(idx),
              "feng_groups": "recomputed", "relabel_config_sha256": fio.sha256_file(rpath), "fallbacks": fired,
              "with_direction": sorted(out), "route": route, "relabel_gate": man.get("gate")}
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
