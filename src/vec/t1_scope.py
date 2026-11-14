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
5. Strong-null down-tail prior (T1 tries T12a and T12b, experiments/harness/TRIES.md; off unless a build asks for it):
   `apply_down_prior` promotes the most highly expressed genes the model calls flat into the predicted bottom tail of
   the change, by binomial thinning through `apply_factors` with each cell's depth kept, in chosen cells only.

Only real cells, count-space edits and log1p CP10k output (PLAN principle 1); never decoded means.

This module is the T1 head of the shared core, vec.core (docs/PLAN.md section 5): it decides the categories, the
shares, the stand-ins and the momentum factors, and the core does the cell work: the allocation (`quotas`), the draws
(`draw_least_used`) and the count-space edits (`counts_of`, `mean_cp10k`, `apply_factors`, `log1p_cp10k`). Those
functions moved there from here and stay importable from here under their old names (`_draw` is `draw_least_used`).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
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
from vec.io import CheckReport, check_submission, dense, make_submission, write_submission
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
                   standin_momentum: bool = False, seed: Seed = 0,
                   logf_filter: Callable[[np.ndarray, sparse.csr_matrix, np.ndarray], tuple[np.ndarray, dict]]
                   | None = None,
                   logf_by_type: Mapping[str, np.ndarray] | None = None,
                   own_factor_cells: np.ndarray | None = None,
                   kernel: str = "independent",
                   kernel_record: dict | None = None) -> tuple[ad.AnnData, ShiftReport]:
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

    logf_filter (default None: the factors as above; heart try H6 of experiments/harness/TRIES.md,
    vec.t2_build.covariance_filter): called once for each type whose factors are computed, after the clip and the
    held-flat zeros and before `apply_factors`, as logf_filter(logf, counts, held) with counts the type's
    reference-stage cells (in `genes` order, as `ref`'s rows of this type) and held the held-flat genes (boolean, by
    gene); it returns (new logf, a record stored as the type's report row's "filter"). Held-flat genes stay at f = 1
    whatever it returns.

    logf_by_type (default None: every type's factor as above, so every existing output is unchanged; heart candidate
    b0_trend, vec.harness.B0_TREND_RULE): type -> a log factor (already damped, in `genes` order, float64) that replaces
    `damping * log((m_ref + eps) / (m_start + eps))` for that type. Everything else applies to it as to the own factor:
    the skip rules (a type skipped for too few cells or stand-ins stays skipped even if the mapping names it), the clip,
    the held-flat zeros, logf_filter and the report row. Types the mapping does not name keep their own factor.

    own_factor_cells (default None: every cell as above; heart arm b0_trend_mix, vec.harness.B0_TREND_MIX_RULE): a
    boolean mask over target_like's cells; a masked cell of a type the mapping names gets the type's own factor (as if
    the mapping did not name it: the same clip and held-flat zeros) instead of the mapped one. All cells go through one
    `apply_factors` call on one stream, so a mapped factor equal to the own one gives the default output bit for bit.
    The shifted types' report rows then add "n_own_factor_cells". Not defined together with logf_filter (raises).

    kernel (default "independent"): `apply_factors`' kernel (heart arms b0_trend_xt: "exact_total", the repair form;
    b0_trend_mn: "multinomial", the whole-cell redraw; T1 arm b0_xt: "exact_total").

    kernel_record (default None: not passed on, every existing call unchanged; T1 arm b0_xt): handed to `apply_factors`
    as its repair_record (with kernel "exact_total": each cell's d and the cells with a moved entry).

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
    if own_factor_cells is not None:
        if logf_filter is not None:
            raise ValueError("own_factor_cells is not defined together with logf_filter")
        own_factor_cells = np.asarray(own_factor_cells)
        if own_factor_cells.dtype != bool or own_factor_cells.shape != (len(t_types),):
            raise ValueError(f"own_factor_cells must be a boolean mask of shape ({len(t_types)},)")
    own_rows: dict[str, np.ndarray] = {}  # type -> its own factor, for masked cells of a mapped type
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
        if logf_by_type is not None and t in logf_by_type:
            logf = np.array(logf_by_type[t], dtype=np.float64)
            if logf.shape != (len(genes),):
                raise ValueError(f"logf_by_type[{t!r}] has shape {logf.shape}, expected ({len(genes)},)")
            if own_factor_cells is not None:
                n_own = int((own_factor_cells & (t_types == t)).sum())
                row["n_own_factor_cells"] = n_own
                if n_own:
                    own = damping * (np.log(mean_cp10k(R[r_rows]) + eps) - np.log(mean_cp10k(S[s_rows]) + eps))
                    if max_abs_log_factor is not None:
                        own = np.clip(own, -max_abs_log_factor, max_abs_log_factor)
                    own[held] = 0.0
                    own_rows[t] = np.exp(own)
        else:
            logf = damping * (np.log(mean_cp10k(R[r_rows]) + eps) - np.log(mean_cp10k(S[s_rows]) + eps))
        if max_abs_log_factor is not None:
            logf = np.clip(logf, -max_abs_log_factor, max_abs_log_factor)
        logf[held] = 0.0
        if logf_filter is not None:
            logf, row["filter"] = logf_filter(logf, R[r_rows], held)
            logf = np.array(logf, dtype=np.float64)
            logf[held] = 0.0
        F[tindex[t]] = np.exp(logf)
        row.update(status="shifted", n_up=int((logf > np.log(1.1)).sum()), n_down=int((logf < -np.log(1.1)).sum()),
                   mean_abs_log_f=float(np.abs(logf).mean()), max_abs_log_f=float(np.abs(logf).max()))
        rows.append(row)
    tidx = np.array([tindex[t] for t in t_types], dtype=np.int64)
    if own_rows:
        F = np.vstack([F, *own_rows.values()])
        for j, t in enumerate(own_rows):
            tidx[own_factor_cells & (t_types == t)] = len(types) + j
    record_kw = {} if kernel_record is None else {"repair_record": kernel_record}
    new, scale, dead = apply_factors(C, F, tidx, held, keep_depth=keep_depth, seed=seed, kernel=kernel, **record_kw)
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
# The strong-null down-tail prior (T1 tries T12a and T12b, experiments/harness/TRIES.md; the captain's rulings of
# 2026-11-10 (a) and (d))
#
# de_score fills the predicted bottom n_dn slots with the most negative dp = pb(prediction) - pb(reference). When a
# transition is "strong-null" (E9.5 -> E10.5: mostly highly expressed genes going down) many highly expressed genes the
# model calls flat are in the truth's down list. The prior promotes the most highly expressed of the genes the model
# calls near the tail or flat into the bottom tail in count space. Genes are ordered by the reference's level alone
# (`DownPriorTable.level_only`, in the test build as in the backtests): nothing is fitted on any transition's truth,
# and the backtests measure the selection the test build uses. Every constant below is the row's; a change is a new row.

DP_LEVEL_EDGES = (100, 200, 400, 800, 1600, 3200)  # bins of the reference's level: rank edges, highest level first
DP_CHANGE_EDGES = (400, 1600)  # bins of the predicted change: rank edges of dp, most negative first (the tail, near it)
DP_CHANGE_TOP = 1600  # the top ranks of dp, predicted up
DP_PROMOTED_LEVELS = 4  # a promoted gene's level bin: one of the first four (level rank below 800), where the
# thinning's first-order solve meets the band at T1's depth (sparser genes' realised change falls short of it; the row)
DP_PROMOTED_CHANGE = (1, 2)  # a promoted gene's change bin: near the tail or flat, never the tail or predicted up
DP_BUDGET = 200  # the genes promoted: a constant of the rule (half the strong-null transition's down list)
DP_BAND = -0.15  # the change of the promoted genes' pseudobulk against the reference: between ranks 100 and 200 of s1's
# committed dp tail and between ranks 200 and 400 of s2's (experiments/research_2026_11/results/t1-tail.json and
# tails.json), so the promoted block sits inside s1's bottom 400 and only partly inside s2's (the row)
DP_F_MIN = 0.5  # no promoted gene is thinned below this factor (its target is then out of reach: it is skipped)
DP_MIN_PB = 0.25  # the scorer's MIN_LFC: a truth down-gene needs pb(reference) >= this, since log1p values are >= 0
DP_ITERATIONS = 3  # fixed-point steps that account for the depth-keeping constant apply_factors applies on top
DP_TOL = 0.02  # a promoted gene's realised dp may sit this far above its band; every other gene may fall by half this
DP_DEPTH_TOL = 0.01  # the covered cells' realised library ratio stays this close to 1 on average
DP_BOTTOM = 400  # the size of the bottom tail the band is meant to enter, reported with the checks (no count of test)
DP_MMD_BOUND = 0.0035  # raw mmd_u: the file is written only if the prior moves it by no more than this (the scorer's
# mmd_unbiased between the file before and after the prior; the same number as the stop rule's S2)
DP_MMD_CELLS = 5000  # at most this many cells of the file enter that check, the same cells before and after


@dataclass(frozen=True)
class DownPriorTable:
    """The promotion score per (level bin, change bin) of T1 tries T12a and T12b. `p` is levels x changes, the levels
    and changes those of `down_prior_bins`, the last level bin holding every gene below the last edge (never promoted).
    Both the test build and the backtests use `level_only`, which depends on the level alone."""

    p: np.ndarray
    level_edges: tuple[int, ...] = DP_LEVEL_EDGES
    change_edges: tuple[int, ...] = DP_CHANGE_EDGES
    change_top: int = DP_CHANGE_TOP

    @classmethod
    def level_only(cls, level_edges: tuple[int, ...] = DP_LEVEL_EDGES, change_edges: tuple[int, ...] = DP_CHANGE_EDGES,
                   change_top: int = DP_CHANGE_TOP) -> DownPriorTable:
        """The rule's order, with no calibration at all (it reads no truth, so it is legal on every pseudo-task and at
        test): the score falls with the level bin and does not depend on the predicted change; the last level bin is
        0."""
        n_l, n_c = len(level_edges) + 1, len(change_edges) + 2
        p = np.tile(1.0 - np.arange(n_l, dtype=np.float64) / max(n_l - 1, 1), (n_c, 1)).T
        return cls(p=p, level_edges=tuple(level_edges), change_edges=tuple(change_edges), change_top=int(change_top))


def down_prior_bins(pb_ref, dp, table: DownPriorTable) -> tuple[np.ndarray, np.ndarray]:
    """Each gene's level bin (its rank by the reference's mean log1p CP10k, highest first, ties by position, against
    the table's level_edges; the last bin holds every rank past the last edge) and change bin (its rank by dp, most
    negative first: 0 the first change edge's ranks, 1 the next, then the middle; the last bin the top change_top ranks
    of dp, predicted up)."""
    pb_ref, dp = np.asarray(pb_ref, dtype=np.float64), np.asarray(dp, dtype=np.float64)
    n = pb_ref.size
    level_rank = np.empty(n, dtype=np.int64)
    level_rank[np.argsort(-pb_ref, kind="stable")] = np.arange(n)
    level = np.searchsorted(table.level_edges, level_rank, side="right")
    change_rank = np.empty(n, dtype=np.int64)
    change_rank[np.argsort(dp, kind="stable")] = np.arange(n)
    change = np.searchsorted(table.change_edges, change_rank, side="right")
    change[change_rank >= n - table.change_top] = len(table.change_edges) + 1
    return level, change


def down_prior_order(table: DownPriorTable, level, change, pb_ref, dp, genes: Sequence[str], held) -> np.ndarray:
    """Positions of the genes the prior may promote, best first (its table's score, then the most negative dp, then the
    gene's name): a level bin among the first DP_PROMOTED_LEVELS and never the table's last, change bin
    DP_PROMOTED_CHANGE, not held, and pb_ref >= DP_MIN_PB."""
    level, change = np.asarray(level), np.asarray(change)
    ok = ((level < min(len(table.level_edges), DP_PROMOTED_LEVELS)) & np.isin(change, DP_PROMOTED_CHANGE)
          & ~np.asarray(held, dtype=bool) & (np.asarray(pb_ref, dtype=np.float64) >= DP_MIN_PB))
    idx = np.flatnonzero(ok)
    score = table.p[level[idx], change[idx]]
    names = np.array([str(genes[i]) for i in idx], dtype=object).astype(str)
    return idx[np.lexsort((names, np.asarray(dp, dtype=np.float64)[idx], -score))]


@dataclass(frozen=True)
class DownPrior:
    """A build's strong-null prior: the promotion order and which cells it touches (`covered_labels`: None for the cells
    a direction step covers, M1 at test; otherwise the labels of those cells). The other fields are the row's
    constants."""

    table: DownPriorTable
    covered_labels: frozenset[str] | None = None
    budget: int = DP_BUDGET
    band: float = DP_BAND
    f_min: float = DP_F_MIN
    tol: float = DP_TOL
    mmd_bound: float = DP_MMD_BOUND


def _prepost_distance(pre, post, covered: np.ndarray, seed: Seed) -> tuple[float, float]:
    """How far the prior moves the file, with no truth read (T12's part (d)): the scorer's own mmd_unbiased(post, pre)
    on the panel, and the norm of the covered cells' mean shift in 30 PCs fitted on `pre`. Both on at most DP_MMD_CELLS
    cells, the same ones before and after; pre and post are the panel's log1p CP10k, cells x genes. Returns (mmd, shift)
    and records no cell count."""
    from common.core_metrics import mmd_unbiased  # the pinned scorer, as the harness scores with it
    from sklearn.decomposition import PCA

    n = pre.shape[0]
    pick = np.arange(n)
    if n > DP_MMD_CELLS:
        pick = np.sort(np.random.default_rng([int(s) for s in np.atleast_1d(seed)] + [5]).choice(n, DP_MMD_CELLS,
                                                                                              replace=False))
    a, b = dense(pre[pick]).astype(np.float32), dense(post[pick]).astype(np.float32)
    n_pc = min(30, a.shape[0] - 1, a.shape[1])
    mmd = float(mmd_unbiased(b, a, n_pc=n_pc))
    pca = PCA(n_components=n_pc, random_state=0).fit(a)
    cov = np.asarray(covered, dtype=bool)[pick]
    shift = float(np.linalg.norm(pca.components_ @ (b[cov].mean(0) - a[cov].mean(0)))) if cov.any() else 0.0
    return mmd, shift


def apply_down_prior(pop: ad.AnnData, pb_ref: np.ndarray, genes: Sequence[str], covered: np.ndarray,
                     held: Iterable[str], prior: DownPrior, seed: Seed) -> tuple[ad.AnnData, dict]:
    """Promote the prior's genes into the predicted bottom tail of the change, in the covered cells.

    genes: the panel the scorer ranks (the board's genes); pb_ref: the reference's mean log1p CP10k per panel gene (its
    full file, what the scorer's reference copy samples); covered: which rows of `pop` the prior may edit; held: genes
    left exactly as they are, outside the depth-keeping constant. dp = pb(pop) - pb_ref over the panel, all cells. The
    genes promoted are the first `budget` of `down_prior_order` (the table's order: the reference's level alone, then
    the most negative dp) whose target is in reach. A promoted gene moves to the band: in the covered cells its counts
    are thinned by one factor, solved so that the population's mean log1p CP10k falls by band - dp, its dp landing at
    the prior's `band` (the slope is the exact first-order change of the mean log1p under binomial thinning at
    f = 1; the fixed-point steps account for the depth-keeping constant `apply_factors` adds). Factors below `f_min` are
    out of reach and the gene is skipped for the next. Nothing else is chosen: other genes in covered cells move only by
    that constant. Each covered cell's expected library is kept (apply_factors keep_depth), uncovered cells are returned
    unchanged, held genes are untouched.
    Returns (pop edited, report). The report holds gene names, the band and the checks, the prior's distance between the
    file before and after among them (`_prepost_distance`); no cell count."""
    genes = [str(g) for g in genes]
    cols = pd.Index(pop.var_names).get_indexer(genes)
    if (cols < 0).any():
        raise KeyError(f"{int((cols < 0).sum())} panel genes are missing from the prediction")
    covered = np.asarray(covered, dtype=bool)
    if not covered.any():
        raise ValueError("the down-tail prior's covered cells are empty")
    cov_idx, rest_idx = np.flatnonzero(covered), np.flatnonzero(~covered)
    C = counts_of(pop)
    n = C.shape[0]
    held_set = {str(g) for g in held}
    held_all = np.array([g in held_set for g in pop.var_names], dtype=bool)
    held_panel = held_all[cols]
    lib = np.asarray(C.sum(1)).ravel().astype(np.float64)
    a = np.divide(TARGET_SUM, lib, out=np.zeros(n), where=lib > 0)

    def panel_mean(x) -> np.ndarray:
        return np.asarray(sparse.csr_matrix(x[:, cols], dtype=np.float64).sum(axis=0)).ravel() / n

    pb_before = panel_mean(pop.X)
    dp = pb_before - np.asarray(pb_ref, dtype=np.float64)
    level, change = down_prior_bins(pb_ref, dp, prior.table)
    order = down_prior_order(prior.table, level, change, pb_ref, dp, genes, held_panel)
    target = float(prior.band)
    # slope of the population mean log1p of each ordered gene with respect to log f in the covered cells, at f = 1
    M = C[cov_idx][:, cols[order]].tocsc()
    entry_col = np.repeat(np.arange(M.shape[1]), np.diff(M.indptr))
    c = M.data.astype(np.float64)
    term = c * (np.log1p(a[cov_idx][M.indices] * c) - np.log1p(a[cov_idx][M.indices] * (c - 1.0)))
    slope = np.bincount(entry_col, weights=term, minlength=M.shape[1]) / n
    chosen, lam = [], []
    for pos, g in enumerate(order):
        if len(chosen) == prior.budget:
            break
        if slope[pos] <= 0 or target >= dp[g]:
            continue
        need = (target - dp[g]) / slope[pos]
        if need >= np.log(prior.f_min):
            chosen.append(pos)
            lam.append(need)
    promoted = order[chosen]
    lam = np.asarray(lam, dtype=np.float64)
    # the input factors: apply_factors multiplies them by the cell's depth constant s, so solve for f x s (fixed point)
    C_cov = C[cov_idx]
    Cs = M[:, chosen].toarray().astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        B = Cs * (np.log1p(a[cov_idx][:, None] * Cs) - np.log1p(a[cov_idx][:, None] * (Cs - 1.0)))
    B = np.nan_to_num(np.where(Cs > 0, B, 0.0))
    held_cov = np.asarray(C_cov[:, held_all].sum(1)).ravel().astype(np.float64)
    free = lib[cov_idx] - held_cov
    fin = np.exp(lam)
    for _ in range(DP_ITERATIONS):
        freed = (Cs * (1.0 - fin)).sum(1)
        s = np.divide(free, free - freed, out=np.ones_like(free), where=free - freed > 0)
        adj = (B.T @ np.log(s)) / np.maximum(B.sum(0), 1e-12)
        fin = np.clip(np.exp(lam - adj), prior.f_min, 1.0)
    F = np.ones((1, C.shape[1]))
    F[0, cols[promoted]] = fin
    new_cov, scale, dead = apply_factors(C_cov, F, np.zeros(len(cov_idx), dtype=np.int64), held_all, keep_depth=True,
                                         seed=seed)
    cov_part = pop[cov_idx].copy()
    cov_part.layers["counts"] = new_cov.astype(np.int32)
    cov_part.X = log1p_cp10k(new_cov)
    parts = [(cov_idx, cov_part)]
    if len(rest_idx):
        rest_part = pop[rest_idx].copy()
        rest_part.layers["counts"] = C[rest_idx].astype(np.int32)
        parts.append((rest_idx, rest_part))
    out = _stack_rows(parts, n)
    # the target-free checks (T12's part (d)): the tail moves as designed, the depth stays, only the promoted genes move
    pb_after = panel_mean(out.X)
    dp_after = pb_after - np.asarray(pb_ref, dtype=np.float64)
    new = sparse.csr_matrix(out.layers["counts"])
    band_share = float((dp_after[promoted] <= target + prior.tol).mean()) if len(promoted) else 1.0
    bottom = np.argsort(np.argsort(dp_after, kind="stable"), kind="stable") < DP_BOTTOM
    bottom_share = float(bottom[promoted].mean()) if len(promoted) else 1.0
    new_lib = np.asarray(new[cov_idx].sum(1)).ravel().astype(np.float64)
    freed = (Cs * (1.0 - fin)).sum(1)
    expected_total = scale * (free - freed) + held_cov
    with np.errstate(divide="ignore", invalid="ignore"):
        expected_err = float(np.nanmax(np.abs(expected_total / lib[cov_idx] - 1.0))) if len(cov_idx) else 0.0
        depth_ratio = float(np.nanmean(new_lib / lib[cov_idx])) if len(cov_idx) else 1.0
    other = ~held_panel
    other[promoted] = False
    change_other = pb_after[other] - pb_before[other]
    max_log_s = float(np.log(scale).max()) if len(scale) else 0.0
    other_ok = bool((change_other >= -prior.tol / 2).all() and (change_other <= 1.5 * max_log_s + prior.tol / 2).all())
    uncovered_same = (new[rest_idx] != C[rest_idx]).nnz == 0
    held_same = (new[:, held_all] != C[:, held_all]).nnz == 0
    mmd_pre_post, pc_shift = _prepost_distance(pop.X[:, cols], out.X[:, cols], covered, seed)
    checks = {"n_promoted_ok": len(promoted) == prior.budget,  # the rule's constant, not a floor: a shortfall refuses
              "mmd_pre_post": round(mmd_pre_post, 6), "mmd_pre_post_ok": mmd_pre_post <= prior.mmd_bound,
              "pc_mean_shift_covered": round(pc_shift, 4),  # reported, deciding nothing
              "band_share": round(band_share, 4), "band_ok": band_share >= 0.95,
              "promoted_in_bottom_400_share": round(bottom_share, 4),  # reported, deciding nothing (it moves with s2)
              "depth_ratio_mean": round(depth_ratio, 6), "depth_ok": abs(depth_ratio - 1.0) <= DP_DEPTH_TOL,
              "expected_depth_max_error": expected_err, "expected_depth_ok": expected_err < 1e-6,
              "uncovered_cells_unchanged": bool(uncovered_same), "held_genes_unchanged": bool(held_same),
              "other_genes_change_min": round(float(change_other.min()), 6) if change_other.size else 0.0,
              "other_genes_change_max": round(float(change_other.max()), 6) if change_other.size else 0.0,
              "other_genes_ok": other_ok, "no_cell_was_restored": not bool(dead.any())}
    checks["passed"] = bool(checks["n_promoted_ok"] and checks["mmd_pre_post_ok"] and checks["band_ok"]
                            and checks["depth_ok"] and checks["expected_depth_ok"] and uncovered_same and held_same
                            and other_ok)
    report = {"promoted": [genes[i] for i in promoted], "factors": {genes[i]: round(float(f), 4) for i, f in zip(
        promoted, fin, strict=True)}, "n_promoted": int(len(promoted)), "budget": int(prior.budget),
              "band_target_dp": round(target, 4), "f_min": float(prior.f_min), "keep_depth": True, "checks": checks}
    return out, report


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
                     down_prior: DownPrior | None = None,
                     post_cells: Callable[[ad.AnnData], ad.AnnData] | None = None,
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
    `momentum_shift` (eps, min_cells, max_abs_log_factor, keep_depth, standin_momentum; kernel and kernel_record for
    the T1 arm b0_xt). Without a shift the rows are
    the released values of the drawn cells. `classify` writes obs["category"] on ref and start. herd: herd each
    category's draw towards its pool's mean expression before any shift (`reweight`). down_prior (T1 tries T12a and
    T12b): after every shift, the strong-null prior of `apply_down_prior` (`DownPrior`) on the cells the directions
    cover or on its own labels; the report then holds its record under "down_prior". Off (None): nothing here changes.
    category's draw towards its pool's mean expression before any shift (`reweight`). post_cells (default None:
    nothing): a last edit of the cells, called with the drawn and shifted population (counts in layers["counts"], obs
    from `reweight`) after the directions and the momentum and before the board's panel is taken; it returns the
    population (experiments/builds/build_test.py t1 --g-mom, within-type factors at test, try T11).
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
            elif down_prior is not None:  # R1's held cells stay as drawn; the prior's step reads their counts
                rest.layers["counts"] = counts_of(rest).astype(np.int32)
            parts.append((idx, rest))
        pop = _stack_rows(parts, pop.n_obs)
    elif momentum:
        pop, shift = with_momentum(pop)
    prior_record = None
    if down_prior is not None:
        pop_labels = pop.obs[label_key].astype(str).to_numpy()
        if down_prior.covered_labels is not None:
            cov = np.isin(pop_labels, sorted(down_prior.covered_labels))
        elif directions and direction_groups is not None:
            cov = np.array([direction_groups.get(lab) in directions for lab in pop_labels], dtype=bool)
        else:
            raise ValueError("a down-tail prior without covered labels needs directions and direction_groups")
        ref_cols = pd.Index(ref.var_names).get_indexer(board.genes)
        pb_ref = np.asarray(sparse.csr_matrix(ref.X[:, ref_cols], dtype=np.float64).mean(axis=0)).ravel()
        pop, prior_record = apply_down_prior(pop, pb_ref, list(board.genes), cov, held, down_prior, seed=[seed, 4])
    if post_cells is not None:
        pop = post_cells(pop)
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
    if down_prior is not None:  # likewise
        report["down_prior"] = prior_record
    return sub, report


def write_checked(sub: ad.AnnData, path: str | Path, board: Board) -> CheckReport:
    """Write a submission and run the free format check on the written file; raises ValueError on any error."""
    p = write_submission(sub, path)
    rep = check_submission(p, board)
    if not rep.ok:
        raise ValueError(f"{p}: format check failed: {rep.errors}")
    return rep
