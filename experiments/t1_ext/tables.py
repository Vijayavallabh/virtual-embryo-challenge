#!/usr/bin/env python
"""Row T72's tables step (candidate MX2s; Task 1 only): the MEMBER TABLES of the new atlases, and nothing else. The
consensus, the rescale and the two controls are the practice-read script's (experiments/t1_m1/mx2s_read.py `pregate`,
which needs the member drops); no consensus, variant, composition, level or share is computed here.

At the row's commit, after every guard of cell 8 (row_t72.guards) and with data/derived/t1_ext_t72/ absent (created
exclusively), for each member whose load passed (data/external/t1_ext_t72/reduced/LOAD.json; a failed load is a
failed (p1), recorded and carried to the manifest):
  LABELS: groups.label_configs and groups.label_atlas (row N54's frozen sets and thresholds unchanged, its
    label-to-group map, the skm class on both new atlases, the derived classes on X020's E8.5 samples only); per cell
    the label and the skm score;
  TABLES (S step 5), per atlas a and group g of the atlas's groups_formed with at least MIN_CELLS cells at each end
    (and at least one in every sample of the end, so that each sample has a mean):
      m_end = the equal-weight mean over the end's samples of each sample's mean over the group's cells of the cell's
              CP10k (counts x 10^4 over the cell's total on Gene Expression rows, load.py's total_ge); no downsampling;
      D_a,g = log((m_late + 0.5) / (m_early + 0.5)) x 5 / (late - early) = 5/8 for E8.5 -> E16.5;
      mapped to T1's genes (the released E9.5 file's gene names, its var_names only) by unique symbol: a symbol T1
      holds twice gets no value (NaN), a gene the atlas lacks is NaN; ^Hb[ab]- and ^mt- genes are NaN in D;
  the skm purity figure of (m2) for each new CM member: of the late end's r_cm cells, the share whose skm score is at
    least 0.3 (the clause itself is the pre-gate's).
Outputs (engineer_notes section B): data/derived/t1_ext_t72/tables/<ACC>/<group>.npz with genes (T1 order),
m_early, m_late, D, n_cells_early, n_cells_late, samples_early, samples_late, per_sample_early, per_sample_late
(samples x genes, T1 order); data/derived/t1_ext_t72/MANIFEST.json with, per atlas and group, the table's sha256,
the cell counts per end and sample (or why no table), the skm purity figure (late, CM), the cell-call counts and load
QC per sample, each labeller config's sha256, (p1) per atlas, the yaml's sha256, the commit, clean, tasks ["T1"]. It
prints each table's path and sha256 and the manifest's sha256.

  uv run --frozen --no-sync python -I -B experiments/t1_ext/tables.py --commit <row commit> --released-dir data/raw/T1
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import t1x_common as cm  # noqa: E402

Refused = cm.Refused
PSEUDOCOUNT = 0.5           # X002's pseudocount (experiments/t1_direction/config.yaml)
FENG_BRACKET_DAYS = 5.0     # D x 5 / (late - early): X002's E9.5 -> E14.5 equivalent
MIN_CELLS = 10              # S step 5: a group with fewer cells at an end gets no table
HELD = (r"^Hb[ab]-", r"^mt-")
CP10K = 1e4
T1_GENES_FILE = "E9.5_RNA.h5ad"  # a reference stage of every practice task, never a truth; var_names only
TABLE_KEYS = ("genes", "m_early", "m_late", "D", "n_cells_early", "n_cells_late", "samples_early", "samples_late",
              "per_sample_early", "per_sample_late")


# ------------------------------------------------------------------------------------------------ pure rules


def cp10k_mean(C: sparse.spmatrix, total: np.ndarray) -> np.ndarray:
    """Mean over the rows of C of each row's CP10k (counts x 10^4 over the row's total on Gene Expression rows)."""
    total = np.asarray(total, dtype=np.float64)
    if total.size == 0 or not np.all(np.isfinite(total)) or (total <= 0).any():
        raise Refused("a cell has no Gene Expression total")
    X = sparse.diags(CP10K / total) @ sparse.csr_matrix(C, dtype=np.float64)
    return np.asarray(X.mean(0)).ravel()


def end_mean(C, total, gsm, in_group, samples: Sequence[str]) -> tuple[np.ndarray, np.ndarray, dict]:
    """(m_end, per-sample means, cells per sample): the equal-weight mean over `samples` of each sample's mean over its
    group cells' CP10k."""
    per, n = [], {}
    for s in samples:
        m = np.flatnonzero((gsm == s) & in_group)
        n[s] = int(m.size)
        if m.size == 0:
            return None, None, n
        per.append(cp10k_mean(sparse.csr_matrix(C)[m], total[m]))
    P = np.vstack(per)
    return P.mean(0), P, n


def direction(m_early: np.ndarray, m_late: np.ndarray, early_day: float, late_day: float) -> np.ndarray:
    cm.bracket_ok(early_day, late_day)
    return np.log((np.asarray(m_late) + PSEUDOCOUNT) / (np.asarray(m_early) + PSEUDOCOUNT)) * \
        (FENG_BRACKET_DAYS / (late_day - early_day))


def to_t1(values: np.ndarray, genes: Sequence[str], t1_genes: Sequence[str]) -> np.ndarray:
    """Atlas values (on its own unique symbols) on T1's gene order: a symbol T1 holds twice, or the atlas holds twice
    or lacks, is NaN."""
    genes = [str(g) for g in genes]
    vals, counts = np.unique(genes, return_counts=True)
    atl_dup = set(vals[counts > 1])
    pos = {g: i for i, g in enumerate(genes) if g not in atl_dup}
    tv, tc = np.unique(np.asarray(list(t1_genes), dtype=str), return_counts=True)
    t1_dup = set(tv[tc > 1])
    v = np.asarray(values, dtype=np.float64)
    return np.array([v[pos[g]] if g in pos and g not in t1_dup else np.nan for g in t1_genes], dtype=np.float64)


def held_mask(t1_genes: Sequence[str]) -> np.ndarray:
    rx = [re.compile(p) for p in HELD]
    return np.array([any(r.search(str(g)) for r in rx) for g in t1_genes], dtype=bool)


def member_tables(red: Mapping, lab: Mapping, groups_map: Mapping[str, Sequence[str]], groups: Sequence[str],
                  early: str, late: str, t1_genes: Sequence[str]) -> dict:
    """{group: {"table": {...} or None, "record": {...}}} for one atlas."""
    stage, gsm = np.asarray(red["stage"], dtype=str), np.asarray(red["gsm"], dtype=str)
    labels = np.asarray(lab["labels"], dtype=str)
    total = np.asarray(red["total_ge"], dtype=np.float64)
    ends = {e: sorted(set(gsm[stage == s])) for e, s in (("early", early), ("late", late))}
    held = held_mask(t1_genes)
    out = {}
    for g in groups:
        in_g = np.isin(labels, list(groups_map[g]))
        res = {e: end_mean(red["C"], total, gsm, in_g & (stage == s), ends[e])
               for e, s in (("early", early), ("late", late))}
        rec = {"cells": {e: res[e][2] for e in res}, "cells_per_end": {e: int(sum(res[e][2].values())) for e in res}}
        why = [f"{e}: {rec['cells_per_end'][e]} cells (fewer than {MIN_CELLS})" for e in res
               if rec["cells_per_end"][e] < MIN_CELLS]
        why += [f"{e}: a sample without a cell of the group" for e in res if res[e][0] is None]
        if not ends["early"] or not ends["late"]:
            why.append("an end has no sample")
        if why:
            out[g] = {"table": None, "record": {**rec, "no_table": why}}
            continue
        me, ml = res["early"][0], res["late"][0]
        D = to_t1(direction(me, ml, cm.label_day(early), cm.label_day(late)), red["genes"], t1_genes)
        D[held] = np.nan
        tab = {"genes": np.array(list(t1_genes), dtype=str),
               "m_early": to_t1(me, red["genes"], t1_genes), "m_late": to_t1(ml, red["genes"], t1_genes), "D": D,
               "n_cells_early": np.int64(rec["cells_per_end"]["early"]),
               "n_cells_late": np.int64(rec["cells_per_end"]["late"]),
               "samples_early": np.array(ends["early"], dtype=str), "samples_late": np.array(ends["late"], dtype=str),
               "per_sample_early": np.vstack([to_t1(x, red["genes"], t1_genes) for x in res["early"][1]]),
               "per_sample_late": np.vstack([to_t1(x, red["genes"], t1_genes) for x in res["late"][1]])}
        fin = np.isfinite(D)
        out[g] = {"table": tab, "record": {**rec, "genes_finite": int(fin.sum()),
                                           "genes_nonzero": int(np.sum(fin & (D != 0.0)))}}
    return out


def skm_purity(lab: Mapping, stage: np.ndarray, late: str, threshold: float) -> dict:
    """(m2)'s figure: of the late end's r_cm cells, the number and share whose skm score is at least `threshold`."""
    m = (np.asarray(stage, dtype=str) == late) & (np.asarray(lab["labels"], dtype=str) == "r_cm")
    n = int(m.sum())
    k = int(np.sum(np.asarray(lab["skm"])[m] >= threshold))
    return {"late_r_cm_cells": n, "late_r_cm_skm_at_least": k, "threshold": threshold,
            "share": (k / n) if n else None}


def write_table(path: Path, tab: Mapping) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(cm.writable(path), "wb") as f:
        np.savez(f, **{k: tab[k] for k in TABLE_KEYS})
    return cm.sha256_file(path)


# ------------------------------------------------------------------------------------------------ the step


def t1_genes_from(path: Path) -> list[str]:
    import anndata as ad

    a = ad.read_h5ad(path, backed="r")
    try:
        return [str(g) for g in a.var_names]
    finally:
        a.file.close()


def row_tables(commit: str, released_dir: Path, repo: Path | None = None, guard=None, cfg: Mapping | None = None,
               t1_genes: Sequence[str] | None = None, extra: Mapping | None = None) -> dict:
    """ROW T72'S TABLES STEP (see the module docstring). Returns the manifest."""
    rt, ld, gr = cm.sibling("row_t72"), cm.sibling("load"), cm.sibling("groups")
    repo = Path(repo or cm.REPO)
    g0 = (guard or rt.guards)("tables", commit, repo)
    cfg = cfg or cm.load_atlases(repo / rt.ATLASES_REL)
    ext, der = rt.row_external(repo), rt.row_derived(repo)
    lp = ext / "reduced" / "LOAD.json"
    if not lp.is_file():
        raise Refused(f"{lp} is missing: the load step did not finish")
    lj = json.loads(lp.read_text())
    if lj.get("commit") != commit:
        raise Refused(f"{lp} was written at another commit")
    gf = Path(released_dir) / T1_GENES_FILE
    if t1_genes is None:
        t1_genes = t1_genes_from(gf)
    extra = extra if extra is not None else gr.load_extra(cfg)
    groups_map = gr.base_config(cfg)["groups"]
    rt.make_exclusive(der)
    man = {**g0, "source": "t1_ext_tables_t72", "load_sha256": cm.sha256_file(lp), "p1": lj["p1"],
           "dropped_by_screens": lj.get("dropped_by_screens") or {}, "t1_genes": {"file": T1_GENES_FILE,
                                                                                    "count": len(t1_genes)},
           "constants": {"pseudocount": PSEUDOCOUNT, "feng_bracket_days": FENG_BRACKET_DAYS, "min_cells": MIN_CELLS,
                         "held_in_D": list(HELD), "cp10k": CP10K, "skm_purity_score": gr.SKM_PURITY_SCORE,
                         "skm_genes": list(gr.SKM_GENES)},
           "relabel_base_sha256": cfg["relabel"]["base_sha256"],
           "extra_classes_sha256": cfg["extra_competing"]["derived_sha256"], "atlases": {}}
    if t1_genes is not None and gf.is_file():
        man["t1_genes"]["sha256"] = cm.sha256_file(gf)
    for acc, p1 in sorted(lj["p1"].items()):
        if not p1["passed"]:
            man["atlases"][acc] = {"p1": p1, "groups": {}}
            continue
        a = cm.atlas(acc, cfg)
        red = ld.read_reduced(acc, ext, cfg)  # checked against the pins in force (the known blocker, closed)
        if red["manifest_sha256"] != lj["load"][acc]["manifest_sha256"]:
            raise Refused(f"{acc}: the reduction is not the one LOAD.json records")
        ch = red["manifest"]["chosen"]
        runs = gr.label_configs(a, cfg, extra, ch["early"], ch["late"])
        lab = gr.label_atlas(red, runs)
        tabs = member_tables(red, lab, groups_map, a["groups_formed"], ch["early"], ch["late"], t1_genes)
        rec = {"p1": p1, "register": a.get("register"), "early": ch["early"], "late": ch["late"],
               "factor": FENG_BRACKET_DAYS / (cm.label_day(ch["late"]) - cm.label_day(ch["early"])),
               "labeller_runs": lab["runs"], "reduction_manifest_sha256": red["manifest_sha256"],
               "load_qc": lj["load"][acc]["qc"]["per_sample"], "groups": {}}
        if "CM" in a["groups_formed"]:
            rec["skm_purity_late"] = skm_purity(lab, red["stage"], ch["late"], gr.SKM_PURITY_SCORE)
        for grp, t in tabs.items():
            e = dict(t["record"])
            if t["table"] is not None:
                p = der / "tables" / acc / f"{grp}.npz"
                e.update({"table": str(p.relative_to(der)), "sha256": write_table(p, t["table"])})
            rec["groups"][grp] = e
        man["atlases"][acc] = rec
    cm.dump_json(man, der / "MANIFEST.json")
    return man


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--commit", required=True, help="the row's commit, full sha (HEAD must be it)")
    ap.add_argument("--released-dir", required=True, type=Path, help="the checkout's data/raw/T1 (gene names only)")
    a = ap.parse_args(argv)
    try:
        man = row_tables(a.commit, a.released_dir)
    except Refused as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    rt = cm.sibling("row_t72")
    der = rt.row_derived()
    for acc, x in sorted(man["atlases"].items()):
        if not x["p1"]["passed"]:
            print(f"table\t{acc}\tnone: (p1) failed at the load")
        for grp, e in sorted(x["groups"].items()):
            print(f"table\t{acc}\t{grp}\t" + (f"{e['table']}\tsha256 {e['sha256']}" if "table" in e else
                                             f"none: {'; '.join(e['no_table'])}"))
    print(f"manifest\t{der / 'MANIFEST.json'}\tsha256 {cm.sha256_file(der / 'MANIFEST.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
