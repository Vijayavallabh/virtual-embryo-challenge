#!/usr/bin/env python
"""The per-atlas direction tables and the consensus direction (Task 1 only). The strict form (ruling row R27, candidate
MX2s) deleted the composition block: no share (level or change form) is computed or written by this step.

PER ATLAS (an outside atlas loaded by load.py, labelled by groups.py; or the atlas in use through its ruled file):
  - statistics in the clean reduction's layout (feng2022_geo.group_stats; keys 'stage|<GSM>|<recomputed label>|half'),
    directions by feng_dirs.directions_from_stats on the atlas's derived group map and bracket, the readout's
    pseudocount and min_cells_ext, for the groups atlases.yaml lets it form and the checks pass;
  - held genes set to NaN before the shrinkage: R2's modules (consensus.held_modules) in every member, and, in an
    atlas whose early end is a whole embryo, every competing class's genes and every gene the embryo gene gate refuses;
  - feng_dirs.eb_shrink on the group's own halves (decile basis: mean CP10k over both ends; consensus.eb_bins);
  - rescaled to the atlas in use's 5-day bracket: d * 5 / (late - early). The builder applies a table at
    bracket_fraction(10.5, 12.5) = 0.4, so each atlas's own per-day rate enters at 2 / (late - early) at the test and
    1 / (late - early) on the E9.5 -> E10.5 practice step; both are recorded, and both assume linear calendar time.
THE ATLAS IN USE: its ruled direction (feng_dirs.recomputed_directions, unshrunk, after its fallback) is kept byte for
byte for every group no other atlas forms (k = 1 at group level), and for every group of a variant with fewer than two
members. Its member in a consensus is its own eb_shrunk direction (factor 1, its bracket is 5 days).
CONSENSUS per group g and gene i, fixed now: A = the members with a finite, NON-ZERO value (an exact zero, as eb_shrink
gives a whole decile it floors, counts as no value, like NaN); k = |A|; the majority sign is held by more than half of
A; c = (1/k) * the sum of the agreeing values if at least two agree and they are more than k/2, else 0; NaN where no
member has a finite value. Equal weights per atlas.
VARIANTS, all written, fixed now: 'deciding' = c x (norm of the ruled vector over genes where c is finite and
non-zero) / (norm of c over the same genes), the second review's rule (the challenger's change: the shrinkage and the
1/k mean must not decide by size alone); the record holds the factor, the share of genes at zero and the largest
|deciding| / max over members |member| per group;
'raw' = the consensus as computed; 'noX002' = deciding without the atlas in use; 'noE85start' = deciding without the
members whose early end is before E9.5 (their per-day rate holds E8.5 -> E9.5 heart-tube maturation). Only 'deciding'
is meant to decide a practice row; the others are counted sensitivities. (The specification's own consensus and
rescale, its steps 4 to 6, are other engineers' work after this one.)
READING A REDUCTION: load.read_reduced checks each reduction against the atlas entry and pins in force, the clearance
and the screens record before it is used (the parked branch's known blocker, closed).

Every output goes under --dest (outside the repository) with a MANIFEST.json (sha256 of each output; atlases, GSMs,
stages, derived configs, check outcomes; tasks ["T1"]). Nothing prints a gene-level value.

  uv run --frozen --no-sync python -I -B experiments/t1_ext/tables.py --dest ~/v/vec-ext/t1 --atlases GSE244783 \\
      GSE230531 --released-dir data/raw/T1 --name consensus_v1
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import t1x_common as cm  # noqa: E402

Refused = cm.Refused
VARIANTS = ("deciding", "raw", "noX002", "noE85start")
IN_USE = "GSE193346"


def _imports():
    for p in (cm.REPO / "experiments" / "t1_m1", cm.REPO / "experiments" / "t1_direction"):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    import dirlib as dl
    import feng2022_geo
    import feng_dirs as fd

    return fd, dl, feng2022_geo


# ------------------------------------------------------------------------------------------------ pure rules


def consensus(members: Mapping[str, np.ndarray]) -> dict:
    """The sign rule on one group's members (name -> vector on one gene order)."""
    names = sorted(members)
    X = np.vstack([np.asarray(members[n], dtype=np.float64) for n in names])
    valid = np.isfinite(X) & (X != 0.0)
    k = valid.sum(0)
    pos, neg = (valid & (X > 0)).sum(0), (valid & (X < 0)).sum(0)
    maj = np.where(pos > neg, 1.0, np.where(neg > pos, -1.0, 0.0))
    agree = np.where(maj > 0, pos, np.where(maj < 0, neg, 0))
    ok = (agree >= 2) & (2 * agree > k)
    S = np.where(valid & (np.sign(np.nan_to_num(X)) == maj), X, 0.0).sum(0)
    c = np.where(ok, S / np.maximum(k, 1), 0.0)
    c[~np.isfinite(X).any(0)] = np.nan
    return {"c": c, "k": k, "agree": agree, "members": names}


def norm_match(c: np.ndarray, ref: np.ndarray, members: Mapping[str, np.ndarray] | None = None
               ) -> tuple[np.ndarray, dict]:
    """The second review's rule, fixed before any run: deciding = c x (norm of the ruled vector over the genes where c
    is finite and non-zero, and the ruled vector finite) / (norm of c over the same genes). The record holds the
    factor, the share of finite consensus genes at zero, and the largest |deciding| / max over members |member| (over
    genes where a member is non-zero), so a step piled onto a few genes shows."""
    c, ref = np.asarray(c, dtype=np.float64), np.asarray(ref, dtype=np.float64)
    m = np.isfinite(c) & (c != 0.0) & np.isfinite(ref)
    nc, nr = float(np.linalg.norm(c[m])), float(np.linalg.norm(ref[m]))
    fin = np.isfinite(c)
    rec = {"norm_c": nc, "norm_ref": nr, "genes_in_norm": int(m.sum()),
           "zero_share": float(np.mean(c[fin] == 0.0)) if fin.any() else None, "factor": None,
           "max_ratio_to_largest_member": None}
    if nc == 0.0:
        return c.copy(), rec
    d = c * (nr / nc)
    rec["factor"] = nr / nc
    if members:
        X = np.abs(np.vstack([np.asarray(v, dtype=np.float64) for v in members.values()]))
        with np.errstate(invalid="ignore"):
            big = np.nanmax(np.where(np.isfinite(X), X, np.nan), axis=0) if np.isfinite(X).any() else \
                np.full(X.shape[1], np.nan)
        ok = np.isfinite(d) & np.isfinite(big) & (big > 0)
        rec["max_ratio_to_largest_member"] = float(np.max(np.abs(d[ok]) / big[ok])) if ok.any() else None
    return d, rec


def group_tables(members: Mapping[str, Mapping[str, np.ndarray]], ruled: Mapping[str, np.ndarray],
                 early_day: Mapping[str, float], in_use: str = IN_USE) -> dict:
    """Every variant for every group of the ruled direction. members: atlas -> group -> vector (held, shrunk,
    rescaled); ruled: the atlas in use's ruled (unshrunk) direction. A group with fewer than two members in a variant
    is the ruled vector unchanged."""
    out: dict[str, dict] = {v: {} for v in VARIANTS}
    rec: dict[str, dict] = {v: {} for v in VARIANTS}
    for g, ref in ruled.items():
        for v in VARIANTS:
            mem = {a: d[g] for a, d in members.items() if g in d}
            if v == "noX002":
                mem.pop(in_use, None)
            if v == "noE85start":
                mem = {a: x for a, x in mem.items() if early_day[a] >= 9.5}
            if len(mem) < 2:
                out[v][g] = np.array(ref, dtype=np.float64, copy=True)
                rec[v][g] = {"passthrough": True, "members": sorted(mem)}
                continue
            r = consensus(mem)
            if v == "raw":
                out[v][g], nm = r["c"], None
            else:
                out[v][g], nm = norm_match(r["c"], np.asarray(ref, dtype=np.float64), mem)
            fin = np.isfinite(r["c"])
            rec[v][g] = {"passthrough": False, "members": r["members"], "k_group": len(mem),
                         "genes_with_value": int(fin.sum()), "genes_nonzero": int(np.sum(fin & (r["c"] != 0))),
                         # the second review's point 2: genes set to zero because only one member has a value
                         "genes_zero_one_member": int(np.sum(fin & (r["k"] == 1))),
                         "agree_counts": {str(i): int(np.sum(r["agree"] == i)) for i in range(len(mem) + 1)},
                         "norm_match": nm}
    return {"tables": out, "record": rec}


# ------------------------------------------------------------------------------------------------ members


def held_mask(genes: Sequence[str], patterns: Sequence[str], extra: set[str] = frozenset()) -> np.ndarray:
    rx = [re.compile(p) for p in patterns]
    return np.array([any(r.search(g) for r in rx) or g in extra for g in genes], dtype=bool)


def shrink_member(d: np.ndarray, h: Mapping[str, np.ndarray], held: np.ndarray, factor: float, bins: int) -> tuple:
    fd, _, _ = _imports()
    d, da, db = (np.where(held, np.nan, np.asarray(x, dtype=np.float64)) for x in (d, h["d_a"], h["d_b"]))
    s, factors = fd.eb_shrink(d, da, db, h["mean_cp"], bins=bins)
    fin = np.isfinite(s)
    return s * factor, {"eb_factors": factors, "zero_share": float(np.mean(s[fin] == 0.0)) if fin.any() else None,
                        "held_genes": int(held.sum()), "rescale": factor}


def atlas_member(acc: str, cfg: Mapping, red: Mapping, lab: Mapping, cfg_a: Mapping, t1_genes: Sequence[str],
                 checks_in: Mapping) -> dict:
    """One outside atlas's member: its directions on the groups it may form that pass the checks (no share)."""
    fd, dl, fg = _imports()
    a = cm.atlas(acc, cfg)
    early, late = cfg_a["stages"]
    e_day, l_day = cm.label_day(early), cm.label_day(late)
    cm.bracket_ok(e_day, l_day)
    dir_cfg = fd.load_config()
    st = fg.group_stats(red["C"].tocsr(), red["stage"], red["gsm"], lab["labels"], red["cell"], red["genes"])
    gcfg = fd.groups_config(dir_cfg, cfg_a)
    dirs, cells, _ = fd.directions_from_stats(st, t1_genes, gcfg, (early, late))
    halves = fd.direction_halves(st, t1_genes, gcfg, (early, late))
    groups = [g for g in a["groups_formed"] if g in dirs]
    checks: dict[str, dict] = {}
    G = cfg_a["groups"]

    def prof(stage, labels):
        x = dl.feng_sel(st, stage, labels)
        return None if x is None else {"genes": list(st["genes"]), "cp": x["cp"], "log": x["log"], "n": x["n"]}

    cats = checks_in["categories"]
    own_cat = cfg["checks"]["own_categories"]
    keep = []
    for g in groups:
        c = {}
        early_p, late_p = prof(early, G[g]), prof(late, G[g])
        refs_e = checks_in["released_profiles"].get(early, {})
        own_e = [lab_ for lab_ in refs_e if cats.get(lab_) in own_cat[g]]
        c["a_profile_early"] = (cm.sibling("groups").profile_check(early_p, refs_e, own_e,
                                                                    cfg["checks"]["profile_min_cp10k"])
                                if refs_e else {"passed": False, "reason": "no released profiles at this stage"})
        c["b_profile_late"] = cm.sibling("groups").profile_check(late_p, checks_in["feng_late_profiles"], [g],
                                                                  cfg["checks"]["profile_min_cp10k"])
        mask = (red["stage"] == early) & np.isin(lab["labels"], G[g])
        pg = checks_in.get("purity", {}).get(early)
        if pg is None:
            c["purity"] = {"passed": False, "reason": "released purity means not computed"}
        else:
            c["purity"] = cm.sibling("groups").purity_gate(
                cm.sibling("groups").class_means(lab["scores"], mask, cfg_a), pg["cm"], pg["own"])
        c["min_cells"] = {"cells": cells[g], "passed": min(cells[g][early], cells[g][late]) >=
                          int(cfg["checks"]["min_cells_ext"])}
        c["passed"] = all(c[k]["passed"] for k in ("a_profile_early", "b_profile_late", "purity", "min_cells"))
        checks[g] = c
        if c["passed"]:
            keep.append(g)
    extra_held: set[str] = set()
    gate = None
    if a.get("early_tissue") == "whole_embryo":
        extra_held = {gene for s in cfg_a["sets"].values() if s["role"] == "competing" for gene in s["genes"]}
        heart = set(cm.sibling("groups").TARGET_LABELS)
        nh_labels = sorted(set(lab["labels"]) - heart)
        cm_e, nh = dl.feng_sel(st, early, G["CM"]), dl.feng_sel(st, early, nh_labels)
        idx = dl.sym_index(st["genes"], list(t1_genes))
        gate = (cm.sibling("groups").embryo_gene_gate(dl.on_t1(cm_e["cp"], idx), dl.on_t1(nh["cp"], idx))
                if cm_e is not None and nh is not None else np.zeros(len(t1_genes), dtype=bool))
    held = held_mask(t1_genes, cfg["consensus"]["held_modules"], extra_held)
    factor = float(cfg["consensus"]["feng_bracket_days"]) / (l_day - e_day)
    out, rec = {}, {}
    for g in keep:
        hg = held | (~gate if (gate is not None and g == "CM") else False)
        out[g], rec[g] = shrink_member(dirs[g], halves[g], hg, factor, int(cfg["consensus"]["eb_bins"]))
    return {"atlas": acc, "early": early, "late": late, "early_day": e_day, "late_day": l_day, "directions": out,
            "record": {"groups_formed": a["groups_formed"], "groups_kept": keep, "checks": checks, "shrink": rec,
                       "fractions": {"rescale_to_5_days": factor, "test": 2.0 / (l_day - e_day),
                                     "practice": 1.0 / (l_day - e_day), "assumes": "linear calendar time"},
                       "embryo_gene_gate_genes_kept": None if gate is None else int(gate.sum()),
                       # the second review's point 6: competing classes the purity gate left unchecked, per group
                       "purity_unchecked_classes": {
                           g: sum(1 for x in (c["purity"].get("classes") or {}).values() if not x.get("checked"))
                           for g, c in checks.items()},
                       # per sample: cells, the unassigned split and score quantiles; never the label counts
                       "per_sample": {g: {k: v for k, v in x.items() if k != "labels"}
                                      for g, x in (lab.get("per_sample") or {}).items()}}}


def feng_member(cfg: Mapping, t1_genes: Sequence[str]) -> dict:
    """The atlas in use: its ruled direction (every check of feng_dirs.recomputed_directions), its shrunk member and
    its E14.5 group profiles (check (b)). No label count or share."""
    fd, dl, _ = _imports()
    dir_cfg = fd.load_config()
    ruled, rec = fd.recomputed_directions(t1_genes, dir_cfg)
    rcfg = fd.load_relabel_config()
    st = dl.load_stats(fd.recomputed_stats_path(rcfg))
    halves = fd.direction_halves(st, t1_genes, fd.groups_config(dir_cfg, rcfg))
    held = held_mask(t1_genes, cfg["consensus"]["held_modules"])
    out, srec = {}, {}
    for g in ruled:
        own = rec["fallbacks"].get(g, {}).get("takes") or g
        out[g], srec[g] = shrink_member(ruled[g], halves[own], held, 1.0, int(cfg["consensus"]["eb_bins"]))
    late_profiles = {}
    for g, labs in rcfg["groups"].items():
        x = dl.feng_sel(st, "E14.5", labs)
        if x is not None:
            late_profiles[g] = {"genes": list(st["genes"]), "cp": x["cp"], "log": x["log"], "n": x["n"]}
    return {"atlas": IN_USE, "early": "E9.5", "late": "E14.5", "early_day": 9.5, "late_day": 14.5, "ruled": ruled,
            "directions": out, "late_profiles": late_profiles,
            "record": {"ruled_record": {k: rec[k] for k in ("file", "sha256", "relabel_config_sha256", "fallbacks",
                                                            "with_direction", "route")},
                       "shrink": srec,
                       "fractions": {"rescale_to_5_days": 1.0, "test": 0.4, "practice": 0.2,
                                     "assumes": "linear calendar time"}}}


# ------------------------------------------------------------------------------------------------ the run


def t1_genes_from(path: Path) -> list[str]:
    import anndata as ad

    a = ad.read_h5ad(path, backed="r")
    try:
        return [str(g) for g in a.var_names]
    finally:
        a.file.close()


def released_inputs(cfg: Mapping, released_dir: Path, cfgs: Mapping[str, Mapping]) -> dict:
    """Released T1 training stages E8.5 and E9.5 (never a window stage): label profiles (check (a)) and, per atlas
    config, the released labels' mean class scores (the purity gate). No share."""
    g = cm.sibling("groups")
    cats = g.chamber_categories()
    out = {"categories": cats, "released_profiles": {}, "purity": {}}
    for stage in ("E8.5", "E9.5"):
        C, genes, labels = g.released_stage(released_dir / f"{stage}_RNA.h5ad")
        out["released_profiles"][stage] = g.label_profiles(C, genes, labels)
        for acc, cfg_a in cfgs.items():
            if cfg_a["stages"][0] != stage:
                continue
            import copy as _c

            import feng_relabel as fr

            cr = _c.deepcopy(cfg_a)
            cr["stages"] = [stage]
            res = fr.relabel(C, genes, np.full(len(labels), stage), cr)
            cm_mask = np.array([cats.get(lab) == "cm" for lab in labels])
            own = {}
            for n, s in cr["sets"].items():
                if s["role"] != "competing":
                    continue
                m = np.isin(labels, list(s.get("released_labels") or []))
                own[n] = float(np.mean(res["scores"][n][m])) if m.any() else None
            out["purity"].setdefault(stage, {})[acc] = {
                "cm": {n: float(np.mean(res["scores"][n][cm_mask])) for n in own}, "own": own}
    return out


def run(cfg: Mapping, dest: Path, accs: Sequence[str], released_dir: Path, name: str) -> dict:
    cm.require_frozen(cfg)
    g, ld = cm.sibling("groups"), cm.sibling("load")
    extra = g.load_extra(cfg)
    t1_genes = t1_genes_from(released_dir / "E10.5_RNA.h5ad")
    reds, cfgs, shas = {}, {}, {}
    for acc in accs:
        a = cm.atlas(acc, cfg)
        if a["status"] not in cm.FETCHABLE:
            raise Refused(f"{acc}: status {a['status']!r} is not tabled")
        reds[acc] = ld.read_reduced(acc, dest, cfg)  # checked against the pins in force (the known blocker)
        ch = reds[acc]["manifest"]["chosen"]
        cfgs[acc], shas[acc] = g.derived_config(a, cfg, extra, ch["early"], ch["late"])
    rin = released_inputs(cfg, released_dir, cfgs)
    feng = feng_member(cfg, t1_genes)
    checks_in = {**rin, "feng_late_profiles": feng["late_profiles"]}
    members = {IN_USE: feng}
    for acc in accs:
        lab = g.label_atlas(reds[acc], cfgs[acc])
        ci = {**checks_in, "purity": {s: v.get(acc) for s, v in rin["purity"].items() if v.get(acc)}}
        members[acc] = atlas_member(acc, cfg, reds[acc], lab, cfgs[acc], t1_genes, ci)
    gt = group_tables({a: m["directions"] for a, m in members.items()}, feng["ruled"],
                      {a: m["early_day"] for a, m in members.items()})
    fallbacks = feng["record"]["ruled_record"].get("fallbacks") or {}
    for d in gt["record"].values():  # the second review's point 5: is the atlas in use's member a fallback?
        for grp, rr in d.items():
            rr["in_use_member_fallback"] = fallbacks.get(grp)
    out = cm.outside_repo(Path(dest) / "tables" / name)
    out.mkdir(parents=True, exist_ok=True)
    arrays = {f"{v}|{grp}": vec for v, d in gt["tables"].items() for grp, vec in d.items()}
    for a, m in members.items():
        arrays.update({f"member|{a}|{grp}": vec for grp, vec in m["directions"].items()})
    np.savez_compressed(out / "directions.npz", genes=np.array(t1_genes, dtype=str), **arrays)
    record = {"name": name, "tasks": list(cm.TASKS), "variants": list(VARIANTS), "deciding_variant": "deciding",
              "groups_ruled": sorted(feng["ruled"]), "consensus": gt["record"],
              "members": {a: m["record"] for a, m in members.items()},
              "atlases": {a: {"early": m["early"], "late": m["late"]} for a, m in members.items()},
              "derived_config_sha256": shas, "base_config_sha256": cfg["relabel"]["base_sha256"],
              "extra_classes_sha256": cfg["extra_competing"]["derived_sha256"],
              "reductions": {a: {"manifest_sha256": reds[a]["manifest_sha256"],
                                 "fetch_manifest_sha256": reds[a]["manifest"]["inputs"]["fetch_manifest_sha256"],
                                 "gsms": sorted(set(reds[a]["gsm"]))} for a in accs}}
    cm.dump_json(record, out / "directions.json")
    import filterio as fio

    outputs = {n: {"sha256": cm.sha256_file(out / n)} for n in ("directions.npz", "directions.json")}
    fio.write_manifest(out, {"source": "t1_ext_tables", "name": name, "tasks": list(cm.TASKS), "outputs": outputs,
                             "atlases": [IN_USE, *accs], "atlas_entries": {a: cm.entry_sha256(a, cfg) for a in accs},
                             "provenance": fio.provenance([Path(__file__), HERE / "groups.py", HERE / "load.py",
                                                           HERE / "t1x_common.py"])})
    return {"out": str(out), "outputs": outputs, "kept": {a: m["record"].get("groups_kept", sorted(m["directions"]))
                                                          for a, m in members.items()}}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", required=True, type=Path)
    ap.add_argument("--atlases", nargs="+", required=True)
    ap.add_argument("--released-dir", required=True, type=Path, help="the checkout's data/raw/T1 (training stages)")
    ap.add_argument("--name", default="consensus_v1")
    a = ap.parse_args(argv)
    try:
        res = run(cm.load_atlases(), a.dest, a.atlases, a.released_dir, a.name)
    except Refused as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    print(json.dumps(res, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
