#!/usr/bin/env python
"""Tier tables, the all-knockouts table and leave-one-study-out tables from processed outside-knockout units.

  combine --processed D/processed --refstats D/reference/refstats.npz --out D/tables

Inputs: every <ACC>/QC.json under --processed whose accession is in datasets.yaml with status cleared. A QC that
failed, a unit file whose sha256 differs from its QC record, a unit computed against other reference statistics, a
QC made under other pins than the committed datasets.yaml holds, or any output of a NEVER_FETCH accession refuses the
run. Only the committed datasets.yaml is read (no flag takes another). A dataset or unit with status 'left_out' (no
group passed its guards) is left out and recorded.

Weights, fixed 2026-11-20 (none from a score), all in log space before each table's finish (t3x_tables.finish):
  * inside a study: a knockout gene's value per (group, gene) is the mean over every (dataset, stage) entry of that
    gene in the study (its datasets' passing stages, already averaged per dataset); the study's value is the mean
    over its genes holding the group (a double knockout is its own gene);
  * a study's pooled value per gene is the reference-group-share-weighted mean over its groups;
  * a table averages STUDIES with equal weight, per (group, gene) for the per-group form and per gene for the pooled
    form, over the studies holding it.
Sign rule (the challenger's change 4), in a table of 2 or more studies: the mean is kept only if at least 2 studies
share its sign and no study of the opposite sign exceeds the mean in magnitude; otherwise the value is 0 (no change).
A cell held by one study in such a table is therefore 0. A table of one study is that study's values.
Tables: wnt, cardiac, unrelated (by datasets.yaml tier), all (every cleared study), one per study, and for every
table of 2 or more studies one leave-one-study-out table per study ('<table>__loso_<study>'). Each gets a manifest
with its studies, genes, per-gene sign agreement and the inputs' sha256. Nothing gene-level is printed.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import t3x_common as cm  # noqa: E402
import t3x_tables as tb  # noqa: E402

Refused = cm.Refused


def sign_rule_mean(V: np.ndarray, n_table_studies: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per column, the mean over the rows holding a value with the sign rule; NaN where no row holds one.
    Returns (values, same-sign counts, holding counts)."""
    V = np.atleast_2d(np.asarray(V, dtype=np.float64))
    held = np.isfinite(V)
    n = held.sum(0)
    m = np.where(n > 0, np.nansum(np.where(held, V, 0.0), 0) / np.maximum(n, 1), np.nan)
    sgn = np.sign(m)
    same = (held & (np.sign(np.where(held, V, 0.0)) == sgn) & (np.where(held, V, 0.0) != 0)).sum(0)
    if n_table_studies < 2:
        return m, same, n
    opp = np.where(held & (np.sign(np.where(held, V, 0.0)) == -sgn), np.abs(np.where(held, V, 0.0)), 0.0).max(0)
    keep = (same >= 2) & (opp <= np.abs(m)) & (sgn != 0)
    out = np.where(n > 0, np.where(keep, m, 0.0), np.nan)
    return out, same, n


def load_units(processed: Path, datasets: Path | None, ref: dict) -> tuple[list[dict], dict]:
    allds = cm.load_datasets(datasets)
    units, record = [], {"left_out": [], "not_cleared": [], "inputs": {}}
    for qcp in sorted(Path(processed).glob("*/QC.json")):
        qc = json.loads(qcp.read_text())
        acc = qc["accession"]
        if acc in cm.NEVER_FETCH:
            raise Refused(f"{acc} is in t3x_common.NEVER_FETCH: a processed output of it must not exist: refused")
        if acc not in allds or allds[acc]["status"] != "cleared":
            record["not_cleared"].append(acc)
            continue
        if qc.get("status") == "failed":
            raise Refused(f"{acc}: its QC failed ({qc.get('reason')}): refused")
        if not allds[acc]["pinned"] or qc.get("pins_sha256") != cm.pins_sha256(allds[acc]):
            raise Refused(f"{acc}: processed under other pins than datasets.yaml holds now: refused")
        record["inputs"][acc] = {"qc": str(qcp), "qc_sha256": cm.sha256_file(qcp)}
        for gene, u in (qc.get("units") or {}).items():
            if u.get("status") == "failed":
                raise Refused(f"{acc} {gene}: unit QC failed: refused")
            if u.get("status") != "ok":
                record["left_out"].append(f"{acc}:{gene}")
                continue
            p = Path(u["unit_file"])
            if not p.is_file() or cm.sha256_file(p) != u["unit_sha256"]:
                raise Refused(f"{acc} {gene}: unit file missing or changed since its QC: refused")
            un = json.loads(p.read_text())
            if un["refstats_sha256"] != ref["file_sha256"] or un["genes"] != ref["genes"]:
                raise Refused(f"{acc} {gene}: computed against other reference statistics: refused")
            L = {g: np.array([np.nan if v is None else v for v in vals], dtype=np.float64)
                 for g, vals in un["L"].items()}
            units.append({"accession": acc, "gene": gene, "study": allds[acc]["study"],
                          "tier": allds[acc]["tier"], "L": L, "sha256": u["unit_sha256"]})
    return units, record


def _nanmean_rows(A: np.ndarray) -> np.ndarray:
    """Column means over the rows holding a value; NaN where none does (no warning)."""
    held = np.isfinite(A)
    n = held.sum(0)
    return np.where(n > 0, np.where(held, A, 0.0).sum(0) / np.maximum(n, 1), np.nan)


def study_tables(units: list[dict], groups: list[str], P: int) -> dict[str, dict]:
    """{study: {'tier', 'genes', 'L': {group: (P,)}}} with equal weight per gene inside the study."""
    out: dict[str, dict] = {}
    for st in sorted({u["study"] for u in units}):
        us = [u for u in units if u["study"] == st]
        tiers = {u["tier"] for u in us}
        if len(tiers) != 1:
            raise Refused(f"study {st} spans tiers {sorted(tiers)}")
        genes = sorted({u["gene"] for u in us})
        Lst = {}
        for g in groups:
            per_gene = []
            for gene in genes:
                vals = [u["L"][g] for u in us if u["gene"] == gene and g in u["L"]]
                if vals:
                    per_gene.append(_nanmean_rows(np.vstack(vals)))
            if per_gene:
                Lst[g] = _nanmean_rows(np.vstack(per_gene))
        out[st] = {"tier": tiers.pop(), "genes": genes, "L": Lst,
                   "units": [f"{u['accession']}:{u['gene']}" for u in us]}
    return out


def table_from_studies(studies: dict[str, dict], names: list[str], ref: dict) -> tuple[dict, np.ndarray, dict]:
    """(per-group log table, pooled log vector, agreement record) over the named studies with the sign rule."""
    P = len(ref["genes"])
    n = len(names)
    per_group = {}
    for g in ref["groups"]:
        rows = [studies[s]["L"].get(g, np.full(P, np.nan)) for s in names]
        if any(np.isfinite(r).any() for r in rows):
            v, _, _ = sign_rule_mean(np.vstack(rows), n)
            if np.isfinite(v).any():
                per_group[g] = v
    pooled_rows = [tb.pool_groups(studies[s]["L"], ref["shares"])[0] if studies[s]["L"] else np.full(P, np.nan)
                   for s in names]
    pooled, same, held = sign_rule_mean(np.vstack(pooled_rows), n)
    agree = {"same_sign_studies": same.tolist(), "holding_studies": held.tolist()}
    return per_group, pooled, agree


def finish_pooled_log(pooled: np.ndarray, per_group: dict, studies: dict, names: list[str], ref: dict) -> np.ndarray:
    used = sorted({g for s in names for g, L in studies[s]["L"].items() if np.isfinite(L).any()})
    if not used:
        return np.ones(len(ref["genes"]))
    den = sum(ref["shares"][g] * ref["mu"][g] for g in used)
    return tb.finish(pooled, den)


def combine(processed: Path, refstats: Path, out: Path, datasets: Path | None = None) -> dict:
    out = cm.outside_repo(out)
    ref = tb.load_ref_stats(refstats)
    units, record = load_units(processed, datasets, ref)
    if not units:
        raise Refused("no unit passed: nothing to combine")
    studies = study_tables(units, ref["groups"], len(ref["genes"]))
    tables = {f"study_{s}": [s] for s in studies}
    for t in cm.TIERS:
        names = [s for s in studies if studies[s]["tier"] == t]
        if names:
            tables[t] = names
    tables["all"] = list(studies)
    for t in [*cm.TIERS, "all"]:
        if t in tables and len(tables[t]) >= 2:
            for s in tables[t]:
                tables[f"{t}__loso_{s}"] = [x for x in tables[t] if x != s]
    written = {}
    for name, names in tables.items():
        per_group_L, pooled_L, agree = table_from_studies(studies, names, ref)
        phi = finish_pooled_log(pooled_L, per_group_L, studies, names, ref)
        pg = tb.finish_per_group(per_group_L, ref)
        man = {"table": name, "studies": names, "tiers": sorted({studies[s]["tier"] for s in names}),
               "genes_by_study": {s: studies[s]["genes"] for s in names},
               "units_by_study": {s: studies[s]["units"] for s in names}, "sign_agreement": agree,
               "sign_rule": "applied" if len(names) >= 2 else "single study", "weights": "equal per study",
               "unit_sha256": {f"{u['accession']}:{u['gene']}": u["sha256"] for u in units if u["study"] in names},
               "refstats_sha256": ref["file_sha256"], "inputs": record,
               "created": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
               "checked": "never" if name.startswith("wnt") else "by try K21 only, after its row",
               "note": ("Wnt tier: biology only, never checked; no released truth can test it" if name.startswith("wnt")
                        else "")}
        written[name] = tb.write_table(out / name, name, ref["genes"], phi, pg, man)
    cm.dump_json({"tables": written, "record": record}, out / "COMBINE.json")
    return {"tables": sorted(written), "record": record}


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--processed", type=Path, required=True)
    ap.add_argument("--refstats", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    try:
        r = combine(a.processed, a.refstats, a.out)
    except Refused as e:
        print(f"refused: {e}")
        return 2
    print(f"combined: {len(r['tables'])} tables: {r['tables']}")
    print(f"left out: {r['record']['left_out']}; not cleared: {r['record']['not_cleared']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
