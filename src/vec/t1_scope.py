"""T1 scope re-weighting and count-space momentum shifts: the B0 floor-plus and the M1 base (docs/PLAN.md section 4.1).

The E12.5 T1 target is a heart-chamber dissection, while the released stages also hold surrounding tissue. Re-weighting
real reference-stage cells to a chamber-only composition is the largest single T1 lever (screened research digest,
gap-fill 1). Predictions are built in four steps, each usable on its own:

1. `classify`: each cell's chamber category from its released label (configs/t1_chamber_classes.yaml). A label that
   is not in the config stops the run; new labels are mapped by the rule in the config's header before use.
2. `reweight`: real cells sampled per category to target shares. Shares are always an input (a dict or a YAML file):
   nothing in this module or its config describes the target stage. Cells are drawn uniformly without replacement
   within a category; a category smaller than its quota reuses cells (as few repeats as possible) and the report says
   so. A category with no cells of its own can be filled by stand-in labels (`fallback`), also reported.
3. `momentum_shift`: the per-type change between two stages, damped, applied to the re-weighted cells in count space.
   Factors are log ratios of mean CP10k with a pseudocount, times `damping`. A decrease is binomial thinning and an
   increase a Poisson top-up of the cell's own counts, so counts stay non-negative integers, exact zeros stay zero and
   held-flat genes keep their counts. By default each cell's other factors are rescaled so that its expected library
   size is unchanged (the change is compositional, like the CP10k means it comes from).
4. `build_submission`: the board's gene order and cell limits, float32 log1p CP10k; `write_checked` writes the file and
   runs `vec.io.check_submission` on it.

Only real cells, count-space edits and log1p CP10k output (PLAN principle 1); never decoded means.

This module is the T1 head of the shared core, vec.core (docs/PLAN.md section 5): it decides the categories, the
shares, the stand-ins and the momentum factors, and the core does the cell work: the allocation (`quotas`), the draws
(`draw_least_used`) and the count-space edits (`counts_of`, `mean_cp10k`, `apply_factors`, `log1p_cp10k`). Those
functions moved there from here and stay importable from here under their old names (`_draw` is `draw_least_used`).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from vec.boards import RELEASED_T1_STAGES, REPO_ROOT, Board
from vec.core import apply_factors, counts_of, draw_least_used, herd_rows, log1p_cp10k, mean_cp10k, quotas
from vec.io import CheckReport, check_submission, make_submission, write_submission
from vec.prep import TARGET_SUM, lineage_map, recover_counts  # noqa: F401 (kept importable here; the core uses them)

CHAMBER_CLASSES = REPO_ROOT / "configs" / "t1_chamber_classes.yaml"
HELD_FLAT_GENES = REPO_ROOT / "experiments" / "t1_direction" / "results" / "held_flat_genes.csv"
CATEGORY_KEY = "category"
LABEL_KEY = "celltype"
SHARE_TOL = 1e-6
RELEASED_STEMS = frozenset(f"{stage}_RNA" for stage in RELEASED_T1_STAGES)  # released T1 file stems, e.g. E8.5_RNA

Seed = int | Sequence[int]


def _seed_record(seed: Seed):
    return int(seed) if np.isscalar(seed) else [int(s) for s in seed]


# ---------------------------------------------------------------------------------------------------------------------
# Chamber categories


@dataclass(frozen=True)
class ChamberSpec:
    """configs/t1_chamber_classes.yaml after validation (see that file for the meaning of each field)."""

    categories: dict[str, dict]
    heart: tuple[str, ...]
    group_base: dict[str, str | None]
    labels: dict[str, str]
    exceptions: dict[str, str]
    file_overrides: dict[str, dict[str, str]]
    file_fallbacks: dict[str, dict[str, tuple[str, ...]]]
    path: str = ""

    def in_scope(self, scope: str) -> list[str]:
        """Categories with this `scope` ("chamber", "adjacent" or "surrounding"), in config order."""
        return [c for c, v in self.categories.items() if v.get("scope") == scope]

    def check_source(self, source: str | None) -> str:
        """The key of a file's per-file rules: "" for None or "" (no per-file rules), else `source` itself when it is
        a released T1 file stem (e.g. "E8.5_RNA", "E10.5_RNA") or a key of file_overrides / file_fallbacks. Anything
        else raises, since it would silently select no rules (e.g. the stage name "E8.5" for the file E8.5_RNA).
        Config keys that are not released stems are accepted for custom configs; `check_spec` reports them."""
        if not source:
            return ""
        if source in RELEASED_STEMS or source in self.file_overrides or source in self.file_fallbacks:
            return source
        keys = sorted(RELEASED_STEMS | set(self.file_overrides) | set(self.file_fallbacks))
        raise ValueError(f"unknown source {source!r}: pass the file stem (one of {keys}), or \"\" to apply no "
                         "per-file rules")

    def mapping(self, source: str | None = None) -> dict[str, str]:
        """label -> category for a file (its stem selects the per-file overrides; see `check_source`)."""
        return {**self.labels, **self.file_overrides.get(self.check_source(source), {})}

    def fallbacks(self, source: str | None = None) -> dict[str, tuple[str, ...]]:
        """category -> stand-in labels for a file that lacks the category's own cells (see `check_source`)."""
        return dict(self.file_fallbacks.get(self.check_source(source), {}))


@lru_cache
def load_spec(path: str | Path = CHAMBER_CLASSES) -> ChamberSpec:
    """Read and validate the chamber-category config: every category it names must be defined."""
    raw = yaml.safe_load(Path(path).read_text())
    cats = {str(c): dict(v or {}) for c, v in raw["categories"].items()}

    def known(cat, where):
        if cat not in cats:
            raise ValueError(f"{path}: {where} names an undefined category {cat!r}")

    labels = {str(k): str(v) for k, v in raw["labels"].items()}
    for lab, cat in labels.items():
        known(cat, f"label {lab!r}")
    overrides = {str(src): {str(k): str(v) for k, v in (m or {}).items()}
                 for src, m in (raw.get("file_overrides") or {}).items()}
    for src, m in overrides.items():
        for lab, cat in m.items():
            known(cat, f"file_overrides[{src}][{lab}]")
    fallbacks = {}
    for src, m in (raw.get("file_fallbacks") or {}).items():
        for cat, labs in (m or {}).items():
            known(cat, f"file_fallbacks[{src}]")
            if not labs:
                raise ValueError(f"{path}: file_fallbacks[{src}][{cat}] lists no labels")
        fallbacks[str(src)] = {str(c): tuple(str(x) for x in labs) for c, labs in (m or {}).items()}
    base = {str(g): (None if c is None else str(c)) for g, c in (raw.get("group_base") or {}).items()}
    for g, cat in base.items():
        if cat is not None:
            known(cat, f"group_base[{g}]")
    heart = tuple(str(c) for c in raw["heart"])
    for cat in heart:
        known(cat, "heart")
    return ChamberSpec(categories=cats, heart=heart, group_base=base, labels=labels,
                       exceptions={str(k): str(v) for k, v in (raw.get("exceptions") or {}).items()},
                       file_overrides=overrides, file_fallbacks=fallbacks, path=str(path))


def check_spec(spec: ChamberSpec | None = None) -> list[str]:
    """Problems with the config (empty when consistent). Against configs/lineage_map.yaml: labels missing from the
    lineage map, labels whose category differs from their group's base without an entry under `exceptions`, labels
    of a group with a null base that are not explained there either, and stale exceptions. Per-file sections
    (file_overrides, file_fallbacks) keyed by anything but a released T1 file stem, such as the stage name "E10.5" for
    the file E10.5_RNA: the file's stem would select none of their rules. (`ChamberSpec.check_source` still accepts
    such keys, so a custom config passed explicitly may key its sections as it likes.)"""
    spec = spec or load_spec()
    lm = lineage_map()
    problems = []
    for lab, cat in spec.labels.items():
        group = lm["labels"].get(lab)
        if group is None:
            problems.append(f"{lab!r} is not in configs/lineage_map.yaml (map it there first)")
            continue
        base = spec.group_base.get(group)
        if group not in spec.group_base:
            problems.append(f"{lab!r}: lineage group {group!r} has no entry under group_base")
        elif cat != base and lab not in spec.exceptions:
            problems.append(f"{lab!r}: category {cat!r} differs from group {group!r} ({base!r}) and has no exception")
    for lab in spec.exceptions:
        group = lm["labels"].get(lab)
        if lab not in spec.labels:
            problems.append(f"exception for {lab!r}, which is not a listed label")
        elif group is not None and spec.group_base.get(group) == spec.labels[lab]:
            problems.append(f"exception for {lab!r} is stale: its category matches its group's base")
    for section, per_file in (("file_overrides", spec.file_overrides), ("file_fallbacks", spec.file_fallbacks)):
        for key in per_file:
            if key not in RELEASED_STEMS:
                problems.append(f"{section} key {key!r} is not a released T1 file stem ({sorted(RELEASED_STEMS)}): "
                                "the file's stem would select none of its rules")
    return problems


class UnmappedLabelError(ValueError):
    """Labels without a chamber category. Map them in configs/t1_chamber_classes.yaml (rule in its header) first."""

    def __init__(self, labels: Mapping[str, int], source: str | None = None):
        self.labels = dict(labels)
        where = f" in {source}" if source else ""
        shown = ", ".join(f"{lab!r} ({n} cells)" for lab, n in sorted(self.labels.items()))
        super().__init__(f"{len(self.labels)} label(s) without a chamber category{where}: {shown}. Map them in "
                         "configs/t1_chamber_classes.yaml by the rule in its header before using this file.")


def classify(adata: ad.AnnData, label_key: str = LABEL_KEY, source: str | None = None,
             spec: ChamberSpec | None = None, key_added: str | None = CATEGORY_KEY) -> np.ndarray:
    """Chamber category of every cell, from its label. Raises `UnmappedLabelError` if any label is not in the config.

    source: the file stem (e.g. "E8.5_RNA"), which selects that file's overrides; None or "" selects none, and any
    other value that is neither a released T1 file stem nor a key of the config's per-file sections raises ValueError
    (`ChamberSpec.check_source`). key_added: obs column written with the categories (None writes nothing).
    """
    spec = spec or load_spec()
    if label_key not in adata.obs:
        raise KeyError(f"obs has no {label_key!r} column (columns: {list(adata.obs.columns)})")
    labels = adata.obs[label_key].astype(str).to_numpy()
    mapping = spec.mapping(source)
    uniq, counts = np.unique(labels, return_counts=True)
    missing = {str(u): int(n) for u, n in zip(uniq, counts, strict=True) if u not in mapping}
    if missing:
        raise UnmappedLabelError(missing, source)
    cats = np.array([mapping[lab] for lab in labels], dtype=object)
    if key_added:
        adata.obs[key_added] = cats
    return cats


def composition(cells: ad.AnnData | Iterable[str], key: str = CATEGORY_KEY) -> dict[str, float]:
    """Share of each category among `cells` (an AnnData's obs[key], or the categories themselves), sorted by name."""
    cats = cells.obs[key].astype(str).to_numpy() if isinstance(cells, ad.AnnData) else np.asarray(list(cells), str)
    if cats.size == 0:
        return {}
    uniq, counts = np.unique(cats, return_counts=True)
    return {str(c): float(n) / cats.size for c, n in zip(uniq, counts, strict=True)}


# ---------------------------------------------------------------------------------------------------------------------
# Target shares


def load_shares(shares: Mapping[str, float] | str | Path, spec: ChamberSpec | None = None,
                tol: float = SHARE_TOL) -> dict[str, float]:
    """Target share of every chamber category, from a dict or a YAML file (a mapping, or one under the key `shares`).

    Every category named must be defined, every share must lie in [0, 1] and the shares must sum to 1 within `tol`:
    nothing is renormalised silently. Categories left out get 0. There is deliberately no default, because target
    shares describe the target stage and are decided after the release.
    """
    spec = spec or load_spec()
    if isinstance(shares, str | Path):
        raw = yaml.safe_load(Path(shares).read_text())
        shares = raw["shares"] if isinstance(raw, dict) and "shares" in raw else raw
    if not isinstance(shares, Mapping) or not shares:
        raise ValueError("target shares must be a non-empty mapping of category -> share")
    unknown = sorted(set(map(str, shares)) - set(spec.categories))
    if unknown:
        raise ValueError(f"unknown categories in the target shares: {unknown}; defined: {list(spec.categories)}")
    out = {c: float(shares.get(c, 0.0)) for c in spec.categories}
    bad = {c: v for c, v in out.items() if not (np.isfinite(v) and 0.0 <= v <= 1.0)}
    if bad:
        raise ValueError(f"shares must lie in [0, 1]: {bad}")
    total = sum(out.values())
    if abs(total - 1.0) > tol:
        raise ValueError(f"target shares sum to {total:.6f}, not 1 (they are never renormalised silently)")
    return out


def share_of(shares: Mapping[str, float], category: str, within: Iterable[str]) -> float:
    """`category`'s share of the joint share of the `within` categories (e.g. the CM share of heart cells)."""
    group = list(within)
    total = sum(float(shares.get(c, 0.0)) for c in group)
    return float(shares.get(category, 0.0)) / total if total > 0 else float("nan")


def with_share(shares: Mapping[str, float] | str | Path, category: str, share: float,
               within: Iterable[str] | None = None, spec: ChamberSpec | None = None) -> dict[str, float]:
    """Shares with `category` at `share` of the `within` categories (default: the config's heart categories).

    The other `within` categories are rescaled in proportion, so their joint total with `category` is unchanged, and
    every category outside `within` keeps its share. The cardiomyocyte hedge is with_share(s, "cm", x), with x
    set after the release.
    """
    spec = spec or load_spec()
    base = load_shares(shares, spec)
    group = list(spec.heart if within is None else within)
    if category not in group:
        raise ValueError(f"{category!r} is not one of the categories it is a share of ({group})")
    if not 0.0 <= share <= 1.0:
        raise ValueError(f"share must lie in [0, 1], got {share}")
    total = sum(base[c] for c in group)
    if total <= 0:
        raise ValueError(f"the categories {group} have no share to divide")
    others = [c for c in group if c != category]
    rest = sum(base[c] for c in others)
    if rest <= 0 and share < 1.0:
        raise ValueError(f"no other category in {group} has a share to take up the remainder")
    out = dict(base)
    out[category] = share * total
    for c in others:
        out[c] = base[c] / rest * (1.0 - share) * total if rest > 0 else 0.0
    return load_shares(out, spec)


# ---------------------------------------------------------------------------------------------------------------------
# Re-weighting


@dataclass
class ReweightReport:
    """What `reweight` drew per category. A category with `reused` > 0 was smaller than its quota (or than what was
    left of a shared stand-in pool) and repeats real cells."""

    n_cells: int
    seed: object
    rows: list[dict] = field(default_factory=list)

    @property
    def reused(self) -> list[str]:
        return [r["category"] for r in self.rows if r["reused"] > 0]

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)

    def as_dict(self) -> dict:
        return {"n_cells": self.n_cells, "seed": self.seed, "reused": self.reused, "rows": self.rows}


def _take(adata: ad.AnnData, idx: np.ndarray) -> ad.AnnData:
    """Rows `idx` of an AnnData (repeats allowed) as a new AnnData with unique obs names cell_0, cell_1, ..."""
    obs = adata.obs.iloc[idx].reset_index(drop=True)
    obs.index = [f"cell_{i}" for i in range(len(idx))]
    layers = {k: v[idx] for k, v in adata.layers.items() if k is not None}  # anndata >= 0.13 lists X as layers[None]
    return ad.AnnData(X=adata.X[idx], obs=obs, var=adata.var.copy(), layers=layers)


_draw = draw_least_used  # its name here before it moved to vec.core (feat-020)


def reweight(adata: ad.AnnData, target_shares: Mapping[str, float] | str | Path, n_cells: int, seed: Seed = 0, *,
             key: str = CATEGORY_KEY, label_key: str = LABEL_KEY, within: Mapping[str, Sequence[str]] | None = None,
             fallback: Mapping[str, Sequence[str]] | None = None,
             spec: ChamberSpec | None = None, herd: bool = False) -> tuple[ad.AnnData, ReweightReport]:
    """Sample `n_cells` real cells so that each category holds exactly its quota of the target shares.

    adata: classified cells (obs[key], see `classify`). target_shares: a dict or YAML file (`load_shares`). Quotas are
    `quotas(shares, n_cells)`, so realised shares equal the targets to within 1/n_cells. within: category -> labels;
    that category is drawn only from these labels (e.g. foregut-type cells as a small surrounding contamination).
    fallback: category -> stand-in labels, used only when the category has no cells of its own (after `within`).
    Categories with own cells draw first, stand-ins after them, so a stand-in pool shared with another category
    repeats cells only once its unused cells run out. Drawn cells keep their obs columns; obs[key] becomes the
    category they fill and obs["source_category"] / obs["source_obs"] record where they came from. Rows are shuffled.
    herd (T1 try T5, experiments/harness/TRIES.md): each category's draw is then herded towards its pool's mean
    expression (`vec.core.herd_rows`, from its own random stream, so the draws it starts from are unchanged).
    """
    spec = spec or load_spec()
    shares = load_shares(target_shares, spec)
    if key not in adata.obs:
        raise KeyError(f"obs has no {key!r} column: run classify() first")
    cats = adata.obs[key].astype(str).to_numpy()
    needs_labels = bool(within) or bool(fallback)
    if needs_labels and label_key not in adata.obs:
        raise KeyError(f"obs has no {label_key!r} column, which `within` and `fallback` need")
    labels = adata.obs[label_key].astype(str).to_numpy() if label_key in adata.obs else None
    for name, m in (("within", within), ("fallback", fallback)):
        unknown = sorted(set(m or {}) - set(spec.categories))
        if unknown:
            raise ValueError(f"{name} names undefined categories {unknown}")
    q = quotas(shares, n_cells)
    rng = np.random.default_rng(seed)
    plan = []
    for c in spec.categories:
        if q[c] == 0:
            continue
        own = np.flatnonzero(cats == c)
        pool, how, labs = own, "own", sorted(set(labels[own])) if labels is not None else []
        if within and c in within:
            pool = own[np.isin(labels[own], list(within[c]))]
            how, labs = "own, restricted", sorted(set(labels[pool]))
        if len(pool) == 0 and fallback and c in fallback:
            pool = np.flatnonzero(np.isin(labels, list(fallback[c])))
            how, labs = "stand-in", sorted(set(labels[pool]))
        if len(pool) == 0:
            raise ValueError(f"category {c!r} has a quota of {q[c]} cells but no cells to draw from "
                             f"(own cells: {len(own)}; give stand-in labels with `fallback`, or set its share to 0)")
        plan.append((how == "stand-in", c, pool, how, labs, len(own)))
    plan.sort(key=lambda p: p[0])  # stable: own pools first in config order, then stand-ins
    used = np.zeros(adata.n_obs, dtype=np.int64)
    picks, roles, rows = [], [], []
    hrng = np.random.default_rng([*np.atleast_1d(seed).tolist(), 9]) if herd else None
    herded = {}
    for _, c, pool, how, labs, n_own in plan:
        idx, n_reused = draw_least_used(pool, q[c], used, rng)
        if herd:
            new, herded[c] = herd_rows(adata.X, pool, idx, hrng)
            np.subtract.at(used, idx, 1)
            np.add.at(used, new, 1)
            idx = new
        picks.append(idx)
        roles.append(np.full(len(idx), c, dtype=object))
        rows.append({"category": c, "target_share": shares[c], "quota": q[c], "own_cells": n_own,
                     "pool": int(len(pool)), "drawn_from": how, "labels": ";".join(labs),
                     "unique": int(len(np.unique(idx))), "reused": n_reused,
                     "realised_share": q[c] / n_cells if n_cells else float("nan"),
                     **({"herd": herded[c]} if herd else {})})
    idx = np.concatenate(picks) if picks else np.array([], dtype=np.int64)
    role = np.concatenate(roles) if roles else np.array([], dtype=object)
    perm = rng.permutation(len(idx))
    idx, role = idx[perm], role[perm]
    out = _take(adata, idx)
    out.obs["source_category"] = cats[idx]
    out.obs["source_obs"] = np.asarray(adata.obs_names)[idx]
    out.obs[key] = role
    order = {c: i for i, c in enumerate(spec.categories)}
    rows.sort(key=lambda r: order[r["category"]])
    return out, ReweightReport(n_cells=int(n_cells), seed=_seed_record(seed), rows=rows)


# ---------------------------------------------------------------------------------------------------------------------
# Count-space momentum shift


def held_flat_genes(path: str | Path = HELD_FLAT_GENES) -> list[str]:
    """The pre-registered held-flat list: sex-linked, Y-linked, strain and pool-effect genes (experiments/t1_direction,
    computed from the two released stages only)."""
    return pd.read_csv(path)["gene"].astype(str).tolist()


def _type_rows(adata: ad.AnnData, t: str, key: str, label_key: str, within, fallback,
               min_cells: int) -> tuple[np.ndarray, str]:
    """Rows standing for type t at one stage: its own cells (restricted to `within[t]` labels when given), or its
    stand-in labels when it has fewer than min_cells own cells and the stand-ins reach min_cells."""
    rows = np.flatnonzero(adata.obs[key].astype(str).to_numpy() == t)
    how = "own"
    labels = adata.obs[label_key].astype(str).to_numpy() if label_key in adata.obs else None
    if within and t in within:
        rows = rows[np.isin(labels[rows], list(within[t]))]
        how = "own, restricted"
    if len(rows) < min_cells and fallback and t in fallback and labels is not None:
        fb = np.flatnonzero(np.isin(labels, list(fallback[t])))
        if len(fb) >= min_cells:
            rows, how = fb, "stand-in " + ";".join(fallback[t])
    return rows, how


@dataclass
class ShiftReport:
    """Per-type factors `momentum_shift` applied. `restored_cells` lost every count to thinning and were restored."""

    damping: float
    eps: float
    min_cells: int
    keep_depth: bool
    max_abs_log_factor: float | None
    n_held_flat: int
    standin_momentum: bool = False
    rows: list[dict] = field(default_factory=list)
    restored_cells: int = 0

    @property
    def skipped(self) -> list[str]:
        return [r["type"] for r in self.rows if r["status"] != "shifted"]

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)

    def as_dict(self) -> dict:
        return {"damping": self.damping, "eps": self.eps, "min_cells": self.min_cells, "keep_depth": self.keep_depth,
                "max_abs_log_factor": self.max_abs_log_factor, "n_held_flat": self.n_held_flat,
                "standin_momentum": self.standin_momentum, "restored_cells": self.restored_cells,
                "skipped": self.skipped, "rows": self.rows}


def momentum_shift(start: ad.AnnData, ref: ad.AnnData, target_like: ad.AnnData, damping: float,
                   held_flat: Iterable[str] = (), *, key: str = CATEGORY_KEY, label_key: str = LABEL_KEY,
                   eps: float = 0.5, min_cells: int = 20, max_abs_log_factor: float | None = None,
                   keep_depth: bool = True, within: Mapping[str, Sequence[str]] | None = None,
                   start_fallback: Mapping[str, Sequence[str]] | None = None,
                   ref_fallback: Mapping[str, Sequence[str]] | None = None,
                   standin_momentum: bool = False, seed: Seed = 0) -> tuple[ad.AnnData, ShiftReport]:
    """Move `target_like` cells by the damped per-type change from `start` to `ref`, in count space.

    start, ref: the two stages whose change is the momentum (start is the earlier one; ref is normally the stage the
    cells came from). target_like: the real cells to edit, normally `reweight`'s output; each cell gets the factors
    of its type (obs[key], all three classified the same way). With m the mean CP10k over a type's cells at a stage,
        log f[type, gene] = damping * log((m_ref + eps) / (m_start + eps)),
    clipped to +-max_abs_log_factor when given, and f = 1 for held-flat genes and for types with fewer than
    `min_cells` cells at either stage. A type's cells at a stage are its own cells, restricted to `within[type]` labels
    when given (as in `reweight`), or the stand-in labels of `start_fallback` / `ref_fallback` when it has fewer than
    `min_cells` own cells there. A type drawn from stand-in labels at either stage gets no momentum (f = 1) unless
    standin_momentum is True: the change between a stand-in and the type's own cells measures a difference of cell
    identity, not development (e.g. E8.5 JCF standing in for epicardium). The factors are applied by `apply_factors`
    (binomial thinning / Poisson top-up of each cell's own counts, zeros kept, held-flat counts unchanged, expected
    depth kept with keep_depth).

    Returns a copy of target_like with layers["counts"] edited and X = log1p CP10k of the new counts (float32 CSR),
    and a per-type report.
    """
    if not np.isfinite(damping):
        raise ValueError(f"damping must be finite, got {damping}")
    if eps <= 0:
        raise ValueError(f"eps must be > 0, got {eps}")
    for name, a in (("start", start), ("ref", ref), ("target_like", target_like)):
        if key not in a.obs:
            raise KeyError(f"{name} has no obs[{key!r}]: run classify() first")
    genes = list(target_like.var_names)
    C = counts_of(target_like)
    S, R = counts_of(start, genes), counts_of(ref, genes)
    held_set = set(map(str, held_flat))
    held = np.array([g in held_set for g in genes], dtype=bool)
    t_types = target_like.obs[key].astype(str).to_numpy()
    types = sorted(set(t_types))
    tindex = {t: i for i, t in enumerate(types)}
    F = np.ones((len(types), len(genes)))
    rows = []
    for t in types:
        s_rows, s_how = _type_rows(start, t, key, label_key, within, start_fallback, min_cells)
        r_rows, r_how = _type_rows(ref, t, key, label_key, within, ref_fallback, min_cells)
        row = {"type": t, "n_target": int((t_types == t).sum()), "n_start": int(len(s_rows)),
               "n_ref": int(len(r_rows)), "start_cells": s_how, "ref_cells": r_how}
        if len(s_rows) < min_cells or len(r_rows) < min_cells:
            row.update(status=f"skipped: fewer than {min_cells} cells at a stage", n_up=0, n_down=0,
                       mean_abs_log_f=0.0, max_abs_log_f=0.0)
            rows.append(row)
            continue
        if not standin_momentum and (s_how.startswith("stand-in") or r_how.startswith("stand-in")):
            row.update(status="skipped: stand-in cells at a stage (identity, not development)", n_up=0, n_down=0,
                       mean_abs_log_f=0.0, max_abs_log_f=0.0)
            rows.append(row)
            continue
        logf = damping * (np.log(mean_cp10k(R[r_rows]) + eps) - np.log(mean_cp10k(S[s_rows]) + eps))
        if max_abs_log_factor is not None:
            logf = np.clip(logf, -max_abs_log_factor, max_abs_log_factor)
        logf[held] = 0.0
        F[tindex[t]] = np.exp(logf)
        row.update(status="shifted", n_up=int((logf > np.log(1.1)).sum()), n_down=int((logf < -np.log(1.1)).sum()),
                   mean_abs_log_f=float(np.abs(logf).mean()), max_abs_log_f=float(np.abs(logf).max()))
        rows.append(row)
    tidx = np.array([tindex[t] for t in t_types], dtype=np.int64)
    new, scale, dead = apply_factors(C, F, tidx, held, keep_depth=keep_depth, seed=seed)
    for row in rows:
        m = t_types == row["type"]
        row["median_depth_scale"] = float(np.median(scale[m])) if m.any() else float("nan")
    out = target_like.copy()
    out.layers["counts"] = new.astype(np.int32)
    out.X = log1p_cp10k(new)
    report = ShiftReport(damping=float(damping), eps=float(eps), min_cells=int(min_cells), keep_depth=bool(keep_depth),
                         max_abs_log_factor=max_abs_log_factor, n_held_flat=int(held.sum()),
                         standin_momentum=bool(standin_momentum), rows=rows, restored_cells=int(dead.sum()))
    return out, report


def external_shift(target_like: ad.AnnData, directions: Mapping[str, np.ndarray], fraction: float,
                   group_of: Mapping[str, str], held_flat: Iterable[str] = (), *, label_key: str = LABEL_KEY,
                   max_abs_log_factor: float | None = None, keep_depth: bool = True,
                   seed: Seed = 0) -> tuple[ad.AnnData, ShiftReport]:
    """Move cells along per-group external directions, in count space (M1's direction step, docs/PLAN.md 4.1).

    directions: group -> log ratio of mean CP10k between the two ends of an external bracket (e.g. Feng 2022 E9.5 ->
    E14.5), on target_like's gene order; NaN where the gene has no external value. fraction: the share of the bracket
    the step covers, a time fraction (E10.5 -> E12.5 in the E9.5 -> E14.5 bracket: 2/5). group_of: cell label
    (obs[label_key]) -> group. A cell whose label maps to a group with a direction gets
        log f[gene] = fraction * direction[group][gene]   (0 where NaN and for held-flat genes),
    clipped to +-max_abs_log_factor when given, applied by `apply_factors` like `momentum_shift`. Every other cell is
    returned with its counts unchanged. Within-dataset ratios only: no external level enters a cell.

    Returns a copy of target_like with layers["counts"] edited and X = log1p CP10k of the new counts (float32 CSR),
    and a per-group report (ShiftReport, damping = fraction).
    """
    if not np.isfinite(fraction):
        raise ValueError(f"fraction must be finite, got {fraction}")
    if label_key not in target_like.obs:
        raise KeyError(f"target_like has no obs[{label_key!r}]")
    genes = list(target_like.var_names)
    for g, v in directions.items():
        if np.shape(v) != (len(genes),):
            raise ValueError(f"direction for {g!r} has shape {np.shape(v)}, expected ({len(genes)},)")
    C = counts_of(target_like)
    held_set = set(map(str, held_flat))
    held = np.array([g in held_set for g in genes], dtype=bool)
    labels = target_like.obs[label_key].astype(str).to_numpy()
    cell_group = np.array([group_of.get(lab, "") for lab in labels], dtype=object)
    groups = sorted(g for g in directions if (cell_group == g).any())
    F = np.ones((len(groups) + 1, len(genes)))  # last row: cells without a direction (f = 1)
    rows = []
    for i, g in enumerate(groups):
        logf = fraction * np.nan_to_num(np.asarray(directions[g], dtype=np.float64), nan=0.0)
        if max_abs_log_factor is not None:
            logf = np.clip(logf, -max_abs_log_factor, max_abs_log_factor)
        logf[held] = 0.0
        F[i] = np.exp(logf)
        rows.append({"type": g, "n_target": int((cell_group == g).sum()), "status": "shifted",
                     "labels": ";".join(sorted(set(labels[cell_group == g]))),
                     "n_no_value": int(np.isnan(np.asarray(directions[g], dtype=np.float64)).sum()),
                     "n_up": int((logf > np.log(1.1)).sum()), "n_down": int((logf < -np.log(1.1)).sum()),
                     "mean_abs_log_f": float(np.abs(logf).mean()), "max_abs_log_f": float(np.abs(logf).max())})
    gindex = {g: i for i, g in enumerate(groups)}
    tidx = np.array([gindex.get(g, len(groups)) for g in cell_group], dtype=np.int64)
    untouched = tidx == len(groups)
    rows.append({"type": "(no direction)", "n_target": int(untouched.sum()), "status": "unchanged",
                 "labels": ";".join(sorted(set(labels[untouched])))})
    new, scale, dead = apply_factors(C, F, tidx, held, keep_depth=keep_depth, seed=seed)
    for row in rows:
        m = (tidx == gindex[row["type"]]) if row["type"] in gindex else untouched
        row["median_depth_scale"] = float(np.median(scale[m])) if m.any() else float("nan")
    out = target_like.copy()
    out.layers["counts"] = new.astype(np.int32)
    out.X = log1p_cp10k(new)
    report = ShiftReport(damping=float(fraction), eps=float("nan"), min_cells=0, keep_depth=bool(keep_depth),
                         max_abs_log_factor=max_abs_log_factor, n_held_flat=int(held.sum()), rows=rows,
                         restored_cells=int(dead.sum()))
    return out, report


def _stack_rows(parts: Sequence[tuple[np.ndarray, ad.AnnData]], n: int) -> ad.AnnData:
    """Reassemble cell subsets (original row positions, AnnData) into one AnnData in the original row order."""
    order = np.concatenate([idx for idx, _ in parts])
    if len(order) != n or len(np.unique(order)) != n:
        raise ValueError("the parts do not cover every row exactly once")
    whole = ad.concat([a for _, a in parts], merge="first")
    return whole[np.argsort(order)].copy()


# ---------------------------------------------------------------------------------------------------------------------
# Submission


def build_submission(ref: ad.AnnData, board: Board, target_shares: Mapping[str, float] | str | Path, n_cells: int,
                     seed: int = 0, *, start: ad.AnnData | None = None, damping: float = 0.0,
                     held_flat: Iterable[str] | None = None, ref_source: str | None = None,
                     start_source: str | None = None, label_key: str = LABEL_KEY,
                     within: Mapping[str, Sequence[str]] | None = None,
                     fallback: Mapping[str, Sequence[str]] | None = None, spec: ChamberSpec | None = None,
                     directions: Mapping[str, np.ndarray] | None = None, direction_fraction: float = 0.0,
                     direction_groups: Mapping[str, str] | None = None, herd: bool = False,
                     **shift_kw) -> tuple[ad.AnnData, dict]:
    """A T1 prediction: real `ref` cells re-weighted to `target_shares`, optionally moved by the damped `start` ->
    `ref` momentum, in the board's gene order as float32 log1p CP10k. Returns (AnnData, report).

    directions (M1, docs/PLAN.md 4.1): group -> external log ratio on ref's gene order (see `external_shift`), with
    direction_groups (label -> group) and direction_fraction (the share of the external bracket the step covers).
    Cells whose label maps to a group with a direction move along it; the other cells get the momentum shift (when
    `start` and a damping are given) or stay as drawn. held_flat applies to both.

    ref: the reference stage (loaded with vec.prep.load_stage; counts are recovered from X if absent). n_cells must lie
    within the board's limits (it is an input, like the shares). ref_source / start_source: file stems selecting the
    config's per-file overrides and stand-ins ("E8.5_RNA"; "" applies none; any other value that is neither a released
    T1 file stem nor a key of the config's per-file sections, such as the stage name "E8.5", raises). fallback:
    stand-in labels for categories without own cells in ref (default: the config's for ref_source). Stand-ins fill
    re-weighting quotas; a type drawn from stand-ins at either stage gets no momentum unless standin_momentum=True.
    held_flat: genes the shift leaves alone (default: the pre-registered list, `held_flat_genes`). shift_kw go to
    `momentum_shift` (eps, min_cells, max_abs_log_factor, keep_depth, standin_momentum). Without a shift the rows are
    the released values of the drawn cells. `classify` writes obs["category"] on ref and start. herd: herd each
    category's draw towards its pool's mean expression before any shift (`reweight`).
    """
    if board.task != "T1":
        raise ValueError(f"{board.key} is a {board.task} board; this builder makes T1 predictions")
    if n_cells < board.min_cells or (board.max_cells is not None and n_cells > board.max_cells):
        raise ValueError(f"n_cells {n_cells} is outside the board's limits [{board.min_cells}, {board.max_cells}]")
    if damping != 0 and start is None:
        raise ValueError("a momentum shift (damping != 0) needs the earlier stage `start`")
    spec = spec or load_spec()
    per_file = bool(spec.file_overrides or spec.file_fallbacks)
    if per_file and ref_source is None or (per_file and start is not None and damping != 0 and start_source is None):
        raise ValueError("the chamber config has per-file overrides and stand-ins: pass ref_source (and start_source "
                         "with a shift) as the file stem, or \"\" to apply none")
    for source in (ref_source, start_source):  # before any work, and also a start_source that goes unused
        spec.check_source(source)
    shares = load_shares(target_shares, spec)
    classify(ref, label_key, ref_source, spec)
    fb = spec.fallbacks(ref_source) if fallback is None else dict(fallback)
    pop, rw = reweight(ref, shares, n_cells, seed, label_key=label_key, within=within, fallback=fb, spec=spec,
                       herd=herd)
    shift = ext = None
    momentum = start is not None and damping != 0
    held = (held_flat_genes() if held_flat is None else list(held_flat)) if momentum or directions else []

    def with_momentum(cells: ad.AnnData) -> tuple[ad.AnnData, ShiftReport]:
        classify(start, label_key, start_source, spec)
        return momentum_shift(start, ref, cells, damping, held, label_key=label_key, within=within,
                              start_fallback=spec.fallbacks(start_source), ref_fallback=fb, seed=[seed, 1],
                              **shift_kw)

    if directions:
        if direction_groups is None:
            raise ValueError("directions need direction_groups (cell label -> group)")
        labels = pop.obs[label_key].astype(str).to_numpy()
        covered = np.array([direction_groups.get(lab) in directions for lab in labels], dtype=bool)
        parts = []
        if covered.any():
            idx = np.flatnonzero(covered)
            ext_kw = {k: shift_kw[k] for k in ("max_abs_log_factor", "keep_depth") if k in shift_kw}
            moved, ext = external_shift(pop[idx].copy(), directions, direction_fraction, direction_groups, held,
                                        label_key=label_key, seed=[seed, 2], **ext_kw)
            parts.append((idx, moved))
        if not covered.all():
            idx = np.flatnonzero(~covered)
            rest = pop[idx].copy()
            if momentum:
                rest, shift = with_momentum(rest)
            parts.append((idx, rest))
        pop = _stack_rows(parts, pop.n_obs)
    elif momentum:
        pop, shift = with_momentum(pop)
    cols = pd.Index(pop.var_names).get_indexer(board.genes)
    if (cols < 0).any():
        missing = [g for g, i in zip(board.genes, cols, strict=True) if i < 0]
        raise KeyError(f"{len(missing)} board genes are missing from the reference, e.g. {missing[:5]}")
    X = pop.X[:, cols]
    sub = make_submission(X.toarray() if sparse.issparse(X) else np.asarray(X), board.genes)
    report = {"board": board.key, "n_cells": int(n_cells), "seed": int(seed), "ref_source": ref_source,
              "start_source": start_source, "target_shares": shares, "realised": composition(pop),
              "within": {k: list(v) for k, v in (within or {}).items()},
              "fallback": {k: list(v) for k, v in fb.items()}, "reweight": rw.as_dict(),
              "shift": None if shift is None else shift.as_dict(), "config": spec.path}
    if directions:  # only then, so reports without directions are unchanged (tests/test_core.py golden table)
        report["external"] = None if ext is None else ext.as_dict() | {"fraction": float(direction_fraction)}
    return sub, report


def write_checked(sub: ad.AnnData, path: str | Path, board: Board) -> CheckReport:
    """Write a submission and run the free format check on the written file; raises ValueError on any error."""
    p = write_submission(sub, path)
    rep = check_submission(p, board)
    if not rep.ok:
        raise ValueError(f"{p}: format check failed: {rep.errors}")
    return rep
