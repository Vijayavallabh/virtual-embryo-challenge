#!/usr/bin/env python
"""Tier tables, the all-knockouts table and leave-one-study-out tables from processed outside-knockout units.

  combine --processed D/processed --refstats D/reference/refstats.npz --out D/tables --dest D

Inputs: every <ACC>/QC.json under --processed whose accession is in datasets.yaml with status cleared. A QC that
failed, a unit file whose sha256 differs from its QC record, a unit computed against other reference statistics, a
QC made under other pins than the committed datasets.yaml holds, or any output of a NEVER_FETCH accession refuses the
run. Stage 5: so does a QC whose plan is gone or changed: D/meta/<ACC>/plan.json must exist, pass (status ok, the
committed datasets.yaml, its pins verified) and hash to the QC's plan_sha256. Only the committed datasets.yaml is
read (no flag takes another). A dataset or unit with status 'left_out' (no group passed its guards) is left out and
recorded.

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

Stage 9 (README "Bulk whole-embryo units"): bulk units (key '__whole__', modality bulk) build their own family only,
bulk_all, bulk_study_<s> and bulk_all__loso_<s> (pooled form, finished against the reference's whole-embryo
pseudo-bulk), with the 'two_thirds' sign rule for a table of 3 or more studies; the single-cell family above reads
single-cell units only and is unchanged. A study spanning modalities, a '__whole__' key outside a bulk unit, and one
count file sha256 under two bulk datasets are refused. Each family's manifests carry a record limited to its own
accessions: COMBINE.json is the single-cell family's (byte-identical with and without bulk units, as are its
manifests apart from their 'created' time), COMBINE_bulk.json the bulk family's.
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
import fetch as fe  # noqa: E402  (experiments/t3_ext/fetch.py: the plan loader and the pin check)
import t3x_common as cm  # noqa: E402
import t3x_tables as tb  # noqa: E402

Refused = cm.Refused


SIGN_RULES = ("two_study", "two_thirds")


def sign_rule_mean(V: np.ndarray, n_table_studies: int, branch: str = "two_study"
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per column, the mean over the rows holding a value with the sign rule; NaN where no row holds one.
    Returns (values, same-sign counts, holding counts).

    Stage 9: branch 'two_thirds' is the bulk tables' rule for a table of 3 or more studies: the mean is kept iff the
    studies sharing its non-zero sign number at least max(2, ceil(2h/3)), h the studies holding the gene; else 0. (The
    design's median clause is struck: it is vacuous once two thirds share the mean's sign.) A bulk table of 2 or fewer
    studies, and every single-cell table, keeps the 'two_study' rule."""
    if branch not in SIGN_RULES:
        raise Refused(f"sign rule {branch!r} not in {SIGN_RULES}")
    V = np.atleast_2d(np.asarray(V, dtype=np.float64))
    held = np.isfinite(V)
    n = held.sum(0)
    m = np.where(n > 0, np.nansum(np.where(held, V, 0.0), 0) / np.maximum(n, 1), np.nan)
    sgn = np.sign(m)
    same = (held & (np.sign(np.where(held, V, 0.0)) == sgn) & (np.where(held, V, 0.0) != 0)).sum(0)
    if n_table_studies < 2:
        return m, same, n
    if branch == "two_thirds" and n_table_studies >= 3:
        need = np.maximum(2, np.ceil(2.0 * n / 3.0))
        keep = (same >= need) & (sgn != 0)
        return np.where(n > 0, np.where(keep, m, 0.0), np.nan), same, n
    opp = np.where(held & (np.sign(np.where(held, V, 0.0)) == -sgn), np.abs(np.where(held, V, 0.0)), 0.0).max(0)
    keep = (same >= 2) & (opp <= np.abs(m)) & (sgn != 0)
    out = np.where(n > 0, np.where(keep, m, 0.0), np.nan)
    return out, same, n


def check_plan(acc: str, qc: dict, dest: Path, datasets: Path | None) -> None:
    """Stage 5: a QC counts only while the plan it ran under is still the stored plan: dest/meta/<ACC>/plan.json
    exists, says ok, was made from the committed datasets.yaml, its pins still verify, and it hashes to the QC's
    plan_sha256. Refused otherwise."""
    pp = Path(dest) / "meta" / acc / "plan.json"
    if pp.is_symlink() or not pp.is_file():
        raise Refused(f"{acc}: no stored plan at {pp}: its processed output is refused")
    want = qc.get("plan_sha256")
    if not isinstance(want, str) or cm.sha256_file(pp) != want:
        raise Refused(f"{acc}: its plan no longer passes or changed since processing (plan_sha256): refused")
    try:
        pl = fe._load_plan(acc, Path(dest), datasets)
        fe.verify_pins(cm.dataset(acc, datasets), pl)
    except Refused as e:
        raise Refused(f"{acc}: its plan no longer passes or changed since processing ({e}): refused") from None


def load_units(processed: Path, datasets: Path | None, ref: dict, dest: Path) -> tuple[list[dict], dict]:
    allds = cm.load_datasets(datasets)
    units, record = [], {"left_out": [], "not_cleared": [], "inputs": {}}
    bulk_files: dict[str, str] = {}   # stage 9, G8: a bulk count file's sha256 -> the dataset it was read under
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
        check_plan(acc, qc, dest, datasets)
        mod = cm.modality(allds[acc])
        if (qc.get("modality") or "single_cell") != mod:
            raise Refused(f"{acc}: its QC's modality is not the dataset's ({mod}): refused")
        if mod == "bulk":
            for key, sha in sorted((qc.get("count_files_sha256") or {}).items()):
                other = bulk_files.setdefault(sha, acc)
                if other != acc:
                    raise Refused(f"{acc}: the count file {key} has the sha256 of a count file of {other}: one file "
                                  "under two datasets, refused (G8)")
        record["inputs"][acc] = {"qc": str(qcp), "qc_sha256": cm.sha256_file(qcp), "plan_sha256": qc["plan_sha256"]}
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
            # stage 9: '__whole__' is a bulk unit's only key, and only a bulk unit carries it
            umod = un.get("modality") or "single_cell"
            if umod != mod:
                raise Refused(f"{acc} {gene}: a {umod} unit under a {mod} dataset: refused")
            if (mod == "bulk") != ("__whole__" in un["L"]) or (mod == "bulk" and set(un["L"]) != {"__whole__"}):
                raise Refused(f"{acc} {gene}: the group key '__whole__' is a bulk unit's only key, and a bulk unit's "
                              "only key: refused")
            L = {g: np.array([np.nan if v is None else v for v in vals], dtype=np.float64)
                 for g, vals in un["L"].items()}
            entry = {"accession": acc, "gene": gene, "study": allds[acc]["study"],
                     "tier": allds[acc]["tier"], "L": L, "sha256": u["unit_sha256"]}
            if mod == "bulk":
                entry.update(modality="bulk", k=un.get("k"), m=un.get("m"),
                             libraries={"ko": un.get("ko_samples"), "control": un.get("control_samples")},
                             guards={k: v.get("skipped", "passed") if isinstance(v, dict) else v
                                     for k, v in (u.get("stages") or {}).items()})
            units.append(entry)
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
        mods = {u.get("modality", "single_cell") for u in us}
        if len(mods) != 1:   # stage 9: a study is single-cell or bulk, never both
            raise Refused(f"study {st} spans modalities {sorted(mods)}")
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


def combine(processed: Path, refstats: Path, out: Path, datasets: Path | None = None, *, dest: Path) -> dict:
    out = cm.outside_repo(out)
    ref = tb.load_ref_stats(refstats)
    all_units, record = load_units(processed, datasets, ref, cm.outside_repo(dest))
    if not all_units:
        raise Refused("no unit passed: nothing to combine")
    # stage 9: bulk units go to their own family (bulk_tables); the single-cell family below reads only single-cell
    # units, so its names, rule and bytes are as before
    units = [u for u in all_units if u.get("modality", "single_cell") == "single_cell"]
    bulk_units = [u for u in all_units if u.get("modality") == "bulk"]
    mixed = sorted({u["study"] for u in units} & {u["study"] for u in bulk_units})
    if mixed:
        raise Refused(f"studies {mixed} span modalities ['bulk', 'single_cell']: a study is one or the other")
    # round 2 (both reviewers' lesser points): each family gets a record limited to its own accessions, so the
    # single-cell manifests and COMBINE.json are byte-identical with and without bulk units; the bulk family's record
    # goes to COMBINE_bulk.json
    mod_of = {acc: cm.modality(d) for acc, d in cm.load_datasets(datasets).items()}
    mod_of.update({u["accession"]: u.get("modality", "single_cell") for u in all_units})
    sc_record, bulk_record = (family_record(record, mod_of, m) for m in ("single_cell", "bulk"))
    written, bulk_written = {}, {}
    if units:
        written.update(single_cell_tables(units, ref, sc_record, out))
    if bulk_units:
        bulk_written.update(bulk_tables(bulk_units, ref, bulk_record, out))
    cm.dump_json({"tables": written, "record": sc_record}, out / "COMBINE.json")
    (out / "COMBINE_bulk.json").unlink(missing_ok=True)
    if bulk_units:
        cm.dump_json({"tables": bulk_written, "record": bulk_record}, out / "COMBINE_bulk.json")
    return {"tables": sorted({**written, **bulk_written}), "record": record}


def family_record(record: dict, mod_of: dict[str, str], mod: str) -> dict:
    """The combine record limited to one family's accessions (an accession of unknown modality is single-cell, as
    before stage 9)."""
    def mine(acc: str) -> bool:
        return mod_of.get(acc, "single_cell") == mod

    return {"left_out": [x for x in record["left_out"] if mine(x.split(":", 1)[0])],
            "not_cleared": [a for a in record["not_cleared"] if mine(a)],
            "inputs": {a: v for a, v in record["inputs"].items() if mine(a)}}


def single_cell_tables(units: list[dict], ref: dict, record: dict, out: Path) -> dict:
    """The single-cell family (unchanged since stage 1): tier tables, all, one per study, leave-one-study-out."""
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
    return written


def bulk_study_tables(units: list[dict], P: int) -> dict[str, dict]:
    """{study: {'tier', 'genes', 'L': {'__whole__': (P,)}, 'units', 'line_agreement', 'lines_held'}} with equal weight
    per gene inside the study; line_agreement per panel gene is the share of the study's genes whose value has the sign
    of the study mean, and lines_held the number of the study's genes (lines) holding a value for it, the share's
    denominator (records, for reading only)."""
    out: dict[str, dict] = {}
    for st in sorted({u["study"] for u in units}):
        us = [u for u in units if u["study"] == st]
        tiers = {u["tier"] for u in us}
        if len(tiers) != 1:
            raise Refused(f"study {st} spans tiers {sorted(tiers)}")
        genes = sorted({u["gene"] for u in us})
        per_gene = np.vstack([_nanmean_rows(np.vstack([u["L"]["__whole__"] for u in us if u["gene"] == g]))
                              for g in genes])
        L = _nanmean_rows(per_gene)
        held = np.isfinite(per_gene)
        agree = (held & (np.sign(np.where(held, per_gene, 0.0)) == np.sign(L)) & (np.where(held, per_gene, 0.0) != 0)
                 ).sum(0)
        share = np.where(held.sum(0) > 0, agree / np.maximum(held.sum(0), 1), np.nan)
        out[st] = {"tier": tiers.pop(), "genes": genes, "L": {"__whole__": L},
                   "units": [f"{u['accession']}:{u['gene']}" for u in us],
                   "line_agreement": [None if not np.isfinite(v) else float(v) for v in share],
                   # row K24's restated (b4): the number of the study's lines holding each panel gene, the same
                   # count the share above divides by (integers, panel order)
                   "lines_held": [int(v) for v in held.sum(0)]}
    return out


def bulk_tables(units: list[dict], ref: dict, record: dict, out: Path) -> dict:
    """The bulk family (stage 9): bulk_all, bulk_study_<s>, and bulk_all__loso_<s> for every study when bulk_all holds
    2 or more studies; pooled form only (no per-group file), finished against the reference's whole-embryo pseudo-bulk
    (t3x_tables.whole_den). Equal weight per gene inside a study and per study across studies; the sign rule is today's
    for 2 or fewer studies and 'two_thirds' for 3 or more. No bulk tier table and no table mixing the two families."""
    studies = bulk_study_tables(units, len(ref["genes"]))
    tables = {f"bulk_study_{s}": [s] for s in studies}
    tables["bulk_all"] = list(studies)
    if len(studies) >= 2:
        for s in studies:
            tables[f"bulk_all__loso_{s}"] = [x for x in studies if x != s]
    den = tb.whole_den(ref)
    written = {}
    for name, names in tables.items():
        branch = "two_thirds" if len(names) >= 3 else "two_study"
        pooled, same, held = sign_rule_mean(np.vstack([studies[s]["L"]["__whole__"] for s in names]), len(names),
                                            branch)
        phi = tb.finish(pooled, den)
        us = [u for u in units if u["study"] in names]
        man = {"table": name, "studies": names, "modality": "bulk",
               "modality_by_study": {s: "bulk" for s in names},
               "tiers": sorted({studies[s]["tier"] for s in names}),
               "genes_by_study": {s: studies[s]["genes"] for s in names},
               "units_by_study": {s: studies[s]["units"] for s in names},
               "sign_agreement": {"same_sign_studies": same.tolist(), "holding_studies": held.tolist()},
               "sign_rule": ("single study" if len(names) < 2 else branch),
               "weights": "equal per gene inside a study, equal per study",
               "k_by_unit": {f"{u['accession']}:{u['gene']}": u.get("k") for u in us},
               "libraries_by_unit": {f"{u['accession']}:{u['gene']}": u.get("libraries") for u in us},
               "guards_by_unit": {f"{u['accession']}:{u['gene']}": u.get("guards") for u in us},
               "line_agreement_by_study": {s: studies[s]["line_agreement"] for s in names
                                           if len(studies[s]["genes"]) > 1},
               "lines_held_by_study": {s: studies[s]["lines_held"] for s in names if len(studies[s]["genes"]) > 1},
               "unit_sha256": {f"{u['accession']}:{u['gene']}": u["sha256"] for u in us},
               "refstats_sha256": ref["file_sha256"], "inputs": record,
               "created": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
               "checked": "by the practice read named in its TRIES row only",
               "note": "bulk whole-embryo family: composition and developmental delay are guarded and recorded, not "
                       "removed"}
        written[name] = tb.write_table(out / name, name, ref["genes"], phi, None, man)
    return written


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--processed", type=Path, required=True)
    ap.add_argument("--refstats", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--dest", type=Path, required=True, help="the fetch destination holding meta/<ACC>/plan.json")
    a = ap.parse_args(argv)
    try:
        r = combine(a.processed, a.refstats, a.out, dest=a.dest)
    except Refused as e:
        print(f"refused: {e}")
        return 2
    print(f"combined: {len(r['tables'])} tables: {r['tables']}")
    print(f"left out: {r['record']['left_out']}; not cleared: {r['record']['not_cleared']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
