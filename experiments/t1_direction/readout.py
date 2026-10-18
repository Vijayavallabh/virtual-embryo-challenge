"""Pre-registered T1 external-direction readout (docs/PLAN.md section 4.1; the rule and every parameter: config.yaml).

Decides on the T1 E9.5 -> E10.5 transition whether M2 (OT bracketing from our E10.5 to external E14.5) is built:
"build" only if the external direction (Feng 2022 CD1 hearts, E9.5 -> E14.5) beats cross-fitted internal momentum
(T1 E8.5 -> E9.5) by >= 0.05 on both population DE metrics (de_score and de_direction, mean of two folds, matched
coverage) and its cardiomyocyte sign agreement is >= 0.70; otherwise "do not build".

One command after the E10.5 release (DGX, repo root, a few minutes; run it in tmux with a log under runs/logs/):
  uv run python experiments/t1_direction/readout.py data/raw/T1/E10.5_RNA.h5ad
If the file is absent it says so and exits with status 0. A real run also needs: the procedure files committed
(refused otherwise; --allow-uncommitted overrides and is recorded), the Feng source registered in EXTERNAL_DATA.md, and
data/external written by experiments/external (each file is checked against its MANIFEST.json). Check the label
harmonisation first with --labels-only (prints the label -> group map and computes nothing else); a --label-map
override must be decided from that list alone. The decision record and aggregate tables go to results/E10.5/, which
is never overwritten without --overwrite.

Before the release:
  --smoke      target T1 E8.5 -> E9.5 with the legal external stages; no momentum (T1 has no earlier stage), so the
               rule is not evaluated. Writes aggregates to results/smoke_E8.5_E9.5/ and the held-flat list to
               results/held_flat_genes.csv. Rows marked held_flat=no use the screened research report's definitions.
  --rehearsal  the full E10.5 code path with the released E9.5 file standing in for E10.5: a null target that shares
               cells with the momentum half, so its numbers mean nothing. Outputs go to runs/t1_direction/rehearsal.

Method (T1 gene order, all genes):
  groups      cell-type labels harmonised by the config's regex rules (E8.5 Endothelium stands in for Endocardium);
              a group needs >= 20 cells at the start and target stages, and at the previous stage for momentum
  target      change in mean log1p-CP10k per group, start -> target stage
  momentum    previous -> start stage per group. Cross-fitted: previous stage (all cells) -> start-stage half M,
              scored on the target built from the other half T (start half T -> all target-stage cells); folds
              (M, T) = (A, B) and (B, A), metrics averaged over the two folds
  external    Feng E9.5 -> E14.5 (and E15.5) per group: log ratio of mean CP10k (pseudocount 0.5), mapped by symbol.
              JAX brackets are secondary diagnostics, used only when registered in EXTERNAL_DATA.md
  population  start population + sum over groups of (start-stage share x per-group change), composition held,
              scored with vector forms of veckit's de_score and de_direction against the population change.
              matched_coverage (decisive): both predictions use only the groups where both predictors exist
  held flat   sex-linked, Y-linked and strain/pool genes (config) are 0 in every prediction
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "experiments" / "external"))
import dirlib as dl  # noqa: E402
import filterio as fio  # noqa: E402
import t1_stats  # noqa: E402

PROCEDURE_FILES = ["experiments/t1_direction/config.yaml", "experiments/t1_direction/dirlib.py",
                   "experiments/t1_direction/readout.py", "experiments/t1_direction/t1_stats.py",
                   "experiments/external/windows.py", "experiments/external/filterio.py",
                   "experiments/external/feng2022.py", "experiments/external/jax_qiu2024.py",
                   "experiments/external/sources.yaml"]
ROW_COLUMNS = ["type", "predictor", "target", "held_flat", "n_start", "n_target", "n_S", "n_S_all", "sign_agree_S",
               "spearman_S", "n_P", "precision_P", "de_direction_mapped", "de_direction_full", "de_score_full",
               "de_score_chance", "n_mapped"]


def repo_path(p) -> Path:
    p = Path(p).expanduser()
    return p if p.is_absolute() else REPO / p


def stage_stats(path: Path, stage: str, cache_dir: Path, label_key: str, seed: int) -> tuple[dict, str, str]:
    """T1 statistics from committed code, cached under runs/ by the input's sha256."""
    sha = t1_stats.sha256_file(path)
    cache = cache_dir / f"t1_{stage}_{sha[:16]}.npz"
    if cache.exists():
        st = dl.load_stats(cache)
        if str(st["source_sha256"]) == sha and str(st["label_key"]) == label_key and int(st["seed"]) == seed:
            return st, sha, "cache"
    t1_stats.save(t1_stats.compute(path, stage, label_key, seed), cache)
    return dl.load_stats(cache), sha, "computed"


def align_genes(st: dict, genes) -> tuple[dict, str]:
    """Put a stage's statistics on the start stage's gene order. Extra genes are dropped; a missing gene stops the run
    (the population metrics need every gene the scorer sees)."""
    if np.array_equal(st["genes"], genes):
        return st, "same order"
    src = pd.Index(st["genes"])
    if not src.is_unique:
        raise SystemExit("duplicate gene names in a stage file")
    idx = src.get_indexer(genes)
    if (idx < 0).any():
        missing = [g for g, i in zip(genes, idx, strict=True) if i < 0]
        raise SystemExit(f"{len(missing)} start-stage genes missing from another stage, e.g. {missing[:5]}")
    out = dict(st)
    for k in ("S_log", "S_cp", "S_nz"):
        out[k] = st[k][:, idx]
    out["genes"] = genes
    return out, f"reordered to the start-stage gene order ({len(src) - len(genes)} extra genes dropped)"


def registered(cfg: dict) -> dict[str, bool]:
    text = (REPO / "EXTERNAL_DATA.md").read_text()
    return {k: v in text for k, v in cfg["external"]["register_keys"].items()}


def jax_embryos(cfg: dict, smoke: bool) -> list[int]:
    es = set(sum(cfg["jax"]["starts"].values(), [])) | set(sum(cfg["jax"]["lates"].values(), []))
    if smoke:
        es |= set(cfg["jax"]["same_interval"]["early"]) | set(cfg["jax"]["same_interval"]["late"])
    return sorted(es)


def load_external(cfg: dict, genes, use_jax: bool, smoke: bool) -> tuple[dict, dict, dict]:
    """Feng (and optionally JAX) statistics from data/external, each checked against the filter's manifest."""
    ext, idx, info = {}, {}, {}
    fpath = repo_path(cfg["external"]["feng_stats"])
    entry = fio.verify_manifest_output(fpath)
    ext["feng"] = dl.load_stats(fpath)
    idx["feng"] = dl.sym_index(ext["feng"]["genes"], genes)
    info["feng"] = {"path": str(fpath.relative_to(REPO)), "sha256": entry["sha256"],
                    "mapped_t1_genes": int((idx["feng"] >= 0).sum())}
    if use_jax:
        jdir, table = repo_path(cfg["external"]["jax_stats_dir"]), repo_path(cfg["external"]["jax_gene_table"])
        embryos = jax_embryos(cfg, smoke)
        for e in embryos:
            fio.verify_manifest_output(jdir / f"embryo_{e}.npz")
        fio.verify_manifest_output(table)
        J = dl.load_jax(jdir, embryos)
        ext["jax"], idx["jax"] = J, dl.jax_index(J["genes"], genes, table)
        info["jax"] = {"embryos": embryos, "mapped_t1_genes": int((idx["jax"] >= 0).sum())}
    return ext, idx, info


def external_vectors(cfg: dict, group: str, ext: dict, idx: dict, smoke: bool) -> dict[str, np.ndarray]:
    labs = cfg["external_labels"].get(group)
    out: dict[str, np.ndarray] = {}
    if not labs:
        return out
    fl, jle, jll = labs
    eps, mc = cfg["metrics"]["pseudocount_cp10k"], cfg["metrics"]["min_cells_ext"]
    if fl and "feng" in ext:
        f0 = dl.feng_sel(ext["feng"], cfg["feng"]["start"], fl)
        for end in cfg["feng"]["ends"]:
            f1 = dl.feng_sel(ext["feng"], end, fl)
            if f0 is not None and f1 is not None and min(f0["n"], f1["n"]) >= mc:
                out[f"feng_{end}"] = dl.on_t1(dl.lfc_ratio(f0, f1, eps), idx["feng"])
    if "jax" in ext:
        pairs = [(f"jax_{s}_{late}", se, jle, le, jll) for s, se in cfg["jax"]["starts"].items()
                 for late, le in cfg["jax"]["lates"].items()]
        if smoke:
            si = cfg["jax"]["same_interval"]
            pairs.append(("jax_same_interval", si["early"], jle, si["late"], jle))
        for name, e0, l0, e1, l1 in pairs:
            j0, j1 = dl.jax_sel(ext["jax"], e0, l0), dl.jax_sel(ext["jax"], e1, l1)
            if j0 is not None and j1 is not None and min(j0["n"], j1["n"]) >= mc:
                out[name] = dl.on_t1(dl.lfc_ratio(j0, j1, eps), idx["jax"])
    return out


def build_groups(cfg: dict, st: dict, H: dict, ext: dict, idx: dict, smoke: bool) -> dict:
    m, inc, fb = cfg["metrics"], cfg["cm_includes"], cfg["fallback"]
    target_groups = set(H["target"].values())
    if target_groups & set(inc):
        target_groups.add("CM")
    groups = {}
    for g in sorted(target_groups):
        Lt = dl.group_labels(H["target"], g, inc)
        Ls = dl.antecedent_labels(H["start"], g, inc, fb)
        if not Lt or not Ls:
            continue
        S_all, T_all = dl.t1_sel(st["start"], Ls), dl.t1_sel(st["target"], Lt)
        if S_all is None or T_all is None or min(S_all["n"], T_all["n"]) < m["min_cells_t1"]:
            continue
        rec = {"Lt": Lt, "Ls": Ls, "S_all": S_all, "T_all": T_all,
               "S_half": {h: dl.t1_sel(st["start"], Ls, halves=(h,)) for h in "AB"},
               "T_half": {h: dl.t1_sel(st["target"], Lt, halves=(h,)) for h in "AB"}}
        if st.get("prev") is not None:
            Lp = dl.antecedent_labels(H["prev"], g, inc, fb)
            P_all = dl.t1_sel(st["prev"], Lp) if Lp else None
            if P_all is not None and P_all["n"] >= m["min_cells_t1"]:
                rec["Lp"], rec["P_all"] = Lp, P_all
        rec["ext"] = external_vectors(cfg, g, ext, idx, smoke)
        groups[g] = rec
    return groups


def held(v: np.ndarray, flat: np.ndarray) -> np.ndarray:
    v = np.array(v, dtype=float)
    v[flat] = np.nan
    return v


def per_type_rows(cfg: dict, groups: dict, flat: np.ndarray, tag: str) -> list[dict]:
    m = cfg["metrics"]
    thr, mg = m["de_threshold"], m["min_genes_sign"]
    rows = []
    for g, r in groups.items():
        base = {"type": g, "held_flat": tag, "n_start": int(r["S_all"]["n"]), "n_target": int(r["T_all"]["n"])}
        d_true, ref = dl.lfc_log(r["S_all"], r["T_all"]), r["S_all"]["log"]
        preds = dict(r["ext"])
        if "P_all" in r:
            preds["mom"] = dl.lfc_log(r["P_all"], r["S_all"])  # shares the start endpoint with the full target
            if "feng_E14.5" in preds:
                preds["mom+feng"] = dl.zscale(preds["mom"]) + dl.zscale(preds["feng_E14.5"])
        for k, p in preds.items():
            rows.append({**base, "predictor": k, "target": "full", **dl.concord(held(p, flat), d_true, ref, thr, mg)})
        folds: dict[str, list[dict]] = {}
        for M, T in m["folds"]:
            dT, refT = dl.lfc_log(r["S_half"][T], r["T_all"]), r["S_half"][T]["log"]
            xp = dict(r["ext"])
            if "P_all" in r:
                xp["mom_xfit"] = dl.lfc_log(r["P_all"], r["S_half"][M])
                if "feng_E14.5" in xp:
                    xp["mom_xfit+feng"] = dl.zscale(xp["mom_xfit"]) + dl.zscale(xp["feng_E14.5"])
            for k, p in xp.items():
                res = dl.concord(held(p, flat), dT, refT, thr, mg)
                folds.setdefault(k, []).append(res)
                rows.append({**base, "predictor": k, "target": f"xfit_{T}", **res})
        for k, rs in folds.items():
            rows.append({**base, "predictor": k, "target": "xfit_mean", **dl.mean_rows(rs)})
        dA = dl.lfc_log(r["S_half"]["A"], r["T_half"]["A"])
        dB = dl.lfc_log(r["S_half"]["B"], r["T_half"]["B"])
        rep = dl.mean_rows([dl.concord(held(dB, flat), dA, r["S_half"]["A"]["log"], thr, mg),
                            dl.concord(held(dA, flat), dB, r["S_half"]["B"]["log"], thr, mg)])
        rows.append({**base, "predictor": "replicate_half", "target": "half_vs_half", **rep})
    return rows


def _score(pred, dpop, ref, thr) -> tuple[float, float, float]:
    s, ch = dl.de_score_vec(pred, dpop, ref, thr)
    return s, dl.de_direction_vec(pred, dpop, ref), ch


def population_rows(cfg: dict, st: dict, groups: dict, flat: np.ndarray, tag: str) -> list[dict]:
    m = cfg["metrics"]
    thr, skip = m["de_threshold"], set(cfg["population_skip"])
    pg = {g: r for g, r in groups.items() if g not in skip}
    pop_t = dl.t1_sel(st["target"], None)
    rows = []

    # full target: start (all cells) -> target (all cells); 'mom' shares the start endpoint with it
    pop_s = dl.t1_sel(st["start"], None)
    dpop, ref = pop_t["log"] - pop_s["log"], pop_s["log"]
    w = {g: r["S_all"]["n"] / pop_s["n"] for g, r in pg.items()}
    V: dict[str, dict[str, np.ndarray]] = {"oracle": {g: dl.lfc_log(r["S_all"], r["T_all"]) for g, r in pg.items()}}
    for g, r in pg.items():
        for k, p in r["ext"].items():
            V.setdefault(k, {})[g] = held(p, flat)
        if "P_all" in r:
            V.setdefault("mom", {})[g] = held(dl.lfc_log(r["P_all"], r["S_all"]), flat)
    for k, vg in V.items():
        s, d, ch = _score(dl.population_prediction(vg, w, vg, flat), dpop, ref, thr)
        rows.append({"target": "full", "design": "own_coverage", "predictor": k, "held_flat": tag,
                     "groups": ",".join(sorted(vg)), "de_score": s, "de_direction": d, "chance": ch,
                     "n_pop_DE": int((np.abs(dpop) >= thr).sum())})
    if not any("P_all" in r for r in pg.values()):
        return rows

    # cross-fitted folds: momentum from start half M, target from start half T -> all target cells
    per_fold: dict[tuple[str, str], list[dict]] = {}
    for M, T in m["folds"]:
        pop_sT = dl.t1_sel(st["start"], None, halves=(T,))
        dpop, ref = pop_t["log"] - pop_sT["log"], pop_sT["log"]
        w = {g: r["S_half"][T]["n"] / pop_sT["n"] for g, r in pg.items()}
        V = {"oracle": {g: dl.lfc_log(r["S_half"][T], r["T_all"]) for g, r in pg.items()}}
        for g, r in pg.items():
            for k, p in r["ext"].items():
                V.setdefault(k, {})[g] = held(p, flat)
            if "P_all" in r:
                mx = dl.lfc_log(r["P_all"], r["S_half"][M])
                V.setdefault("mom_xfit", {})[g] = held(mx, flat)
                if "feng_E14.5" in r["ext"]:
                    V.setdefault("mom_xfit+feng", {})[g] = held(dl.zscale(mx) + dl.zscale(r["ext"]["feng_E14.5"]), flat)
        base = V.get("mom_xfit", {})
        n_de = int((np.abs(dpop) >= thr).sum())
        own_base = _score(dl.population_prediction(base, w, base, flat), dpop, ref, thr)
        for k, vg in V.items():
            s, d, ch = _score(dl.population_prediction(vg, w, vg, flat), dpop, ref, thr)
            row = {"design": "own_coverage", "predictor": k, "groups": ",".join(sorted(vg)), "de_score": s,
                   "de_direction": d, "baseline_de_score": own_base[0], "baseline_de_direction": own_base[1],
                   "chance": ch, "n_pop_DE": n_de}
            per_fold.setdefault(("own_coverage", k), []).append(row)
            if k in ("mom_xfit", "oracle"):
                continue
            gm = sorted(set(vg) & set(base))
            if not gm:
                continue
            sp, dp, _ = _score(dl.population_prediction(vg, w, gm, flat), dpop, ref, thr)
            sb, db, _ = _score(dl.population_prediction(base, w, gm, flat), dpop, ref, thr)
            so, do, _ = _score(dl.population_prediction(V["oracle"], w, gm, flat), dpop, ref, thr)
            per_fold.setdefault(("matched_coverage", k), []).append(
                {"design": "matched_coverage", "predictor": k, "groups": ",".join(gm), "de_score": sp,
                 "de_direction": dp, "baseline_de_score": sb, "baseline_de_direction": db, "oracle_de_score": so,
                 "oracle_de_direction": do, "chance": ch, "n_pop_DE": n_de})
    for (design, k), fold_rows in per_fold.items():
        for (_, T), r in zip(m["folds"], fold_rows, strict=True):
            rows.append({"target": f"xfit_{T}", "held_flat": tag, **r})
        mean = dl.mean_rows(fold_rows)
        groups_seen = sorted({r["groups"] for r in fold_rows})
        rows.append({"target": "xfit_mean", "design": design, "predictor": k, "held_flat": tag,
                     "groups": groups_seen[0] if len(groups_seen) == 1 else " | ".join(groups_seen), **mean})
    for r in rows:
        if "baseline_de_score" in r:
            r["delta_de_score"] = r["de_score"] - r["baseline_de_score"]
            r["delta_de_direction"] = r["de_direction"] - r["baseline_de_direction"]
    return rows


def report_population(cfg: dict, st8: dict, st9: dict, ext: dict, idx: dict, flat: np.ndarray, tag: str) -> list[dict]:
    """The screened report's population construction (its now_pop.py): E8.5 population + per-group changes weighted by
    the mean of the two stage shares, JAX same-interval changes (embryos 1-12 -> 22)."""
    thr, eps = cfg["metrics"]["de_threshold"], cfg["metrics"]["pseudocount_cp10k"]
    si = cfg["jax"]["same_interval"]
    pop8, pop9 = dl.t1_sel(st8, None), dl.t1_sel(st9, None)
    dpop, ref = pop9["log"] - pop8["log"], pop8["log"]
    cardiac = {"CM", "Endo", "SHF", "Pericardium"}
    preds = {k: np.zeros_like(dpop) for k in ("oracle", "jax_same_interval", "oracle_cardiac",
                                              "jax_same_interval_cardiac")}
    used = []
    for g, (a8, a9, jl) in cfg["smoke"]["report_groups"].items():
        A, B = dl.t1_sel(st8, a8), dl.t1_sel(st9, a9)
        ja, jb = dl.jax_sel(ext["jax"], si["early"], jl), dl.jax_sel(ext["jax"], si["late"], jl)
        if A is None or B is None or ja is None or jb is None:
            continue
        wg = 0.5 * (A["n"] / pop8["n"] + B["n"] / pop9["n"])
        d_true = dl.lfc_log(A, B)
        d_jax = np.nan_to_num(dl.on_t1(dl.lfc_ratio(ja, jb, eps), idx["jax"]), nan=0.0)
        preds["oracle"] += wg * d_true
        preds["jax_same_interval"] += wg * d_jax
        if g in cardiac:
            preds["oracle_cardiac"] += wg * d_true
            preds["jax_same_interval_cardiac"] += wg * d_jax
        used.append(g)
    rows = []
    for k, p in preds.items():
        p = p.copy()
        p[flat] = 0.0
        s, d, ch = _score(p, dpop, ref, thr)
        rows.append({"target": "full", "design": "report_now_pop", "predictor": k, "held_flat": tag,
                     "groups": ",".join(used), "de_score": s, "de_direction": d, "chance": ch,
                     "n_pop_DE": int((np.abs(dpop) >= thr).sum())})
    return rows


def readings_for_decision(cfg: dict, type_rows: list[dict], pop_rows: list[dict]) -> dict:
    rule = cfg["decision"]
    prim, design = rule["primary_predictor"], rule["population_design"]
    sel = [r for r in pop_rows if r["target"] == "xfit_mean" and r["design"] == design and r["predictor"] == prim]
    out: dict = {"design": design}
    if sel:
        r = sel[0]
        out["primary"] = {"de_score": r["de_score"], "de_direction": r["de_direction"], "groups": r["groups"]}
        out["baseline"] = {"de_score": r["baseline_de_score"], "de_direction": r["baseline_de_direction"]}
        out["per_fold"] = [{k: x[k] for k in ("target", "de_score", "de_direction", "baseline_de_score",
                                              "baseline_de_direction")}
                           for x in pop_rows if x["design"] == design and x["predictor"] == prim
                           and x["target"].startswith("xfit_") and x["target"] != "xfit_mean"]
    cm = [r for r in type_rows if r["type"] == rule["cm_group"] and r["predictor"] == prim
          and r["target"] == rule["cm_target"]]
    if cm:
        out["cm_sign_agreement"] = cm[0].get("sign_agree_S", float("nan"))
        out["cm_n_S"] = cm[0]["n_S"]
    return out


def write_table(rows: list[dict], path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    if columns:
        df = df.reindex(columns=[c for c in columns if c in df.columns] + [c for c in df.columns if c not in columns])
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, float_format="%.6g")
    return df


def json_default(o):
    if isinstance(o, np.integer | np.bool_):
        return o.item()
    if isinstance(o, np.floating):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def clean(o):
    """NaN -> None recursively, for strict JSON."""
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, list | tuple):
        return [clean(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def print_summary(type_rows: list[dict], pop_rows: list[dict]) -> None:
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    t = pd.DataFrame(type_rows)
    show = t[t.target.isin(["full", "xfit_mean", "half_vs_half"])]
    cols = ["type", "predictor", "target", "held_flat", "n_start", "n_target", "n_S", "sign_agree_S", "spearman_S",
            "de_direction_full", "de_score_full"]
    print(show.reindex(columns=cols).round(3).to_string(index=False))
    p = pd.DataFrame(pop_rows)
    show = p[~p.target.isin(["xfit_A", "xfit_B"])]
    cols = ["target", "design", "predictor", "held_flat", "de_score", "de_direction", "baseline_de_score",
            "baseline_de_direction", "delta_de_score", "delta_de_direction", "groups"]
    print(show.reindex(columns=cols).round(3).to_string(index=False))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("e105", nargs="?", default="data/raw/T1/E10.5_RNA.h5ad", help="released T1 E10.5 file")
    ap.add_argument("--config", default=str(HERE / "config.yaml"))
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--smoke", action="store_true", help="E8.5 -> E9.5 geometry (no E10.5; rule not evaluated)")
    mode.add_argument("--rehearsal", action="store_true", help="code-path check: E9.5 file stands in for E10.5")
    mode.add_argument("--labels-only", action="store_true", help="print the label harmonisation and stop")
    ap.add_argument("--label-map", help="JSON {label: group} overriding single labels (recorded with its sha256)")
    ap.add_argument("--label-key", help="obs column of the cell-type labels (default: config label_key)")
    ap.add_argument("--no-jax", action="store_true", help="skip the JAX diagnostics")
    ap.add_argument("--allow-uncommitted", action="store_true", help="run although procedure files are uncommitted")
    ap.add_argument("--overwrite", action="store_true", help="replace an existing E10.5 decision record")
    ap.add_argument("--out", help="output directory (default from the config, by mode)")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(Path(a.config).read_text())
    run_mode = "smoke" if a.smoke else "rehearsal" if a.rehearsal else "e10.5"
    e105 = repo_path(a.e105)
    if run_mode == "e10.5" and not e105.exists():
        print(f"E10.5 file not found: {e105}\nNothing to do: the readout runs once the validation stage is released "
              f"(2026-10-20). Fetch it with `vec data fetch`, then pass its path.")
        return 0
    label_key = a.label_key or cfg["label_key"]
    override = json.loads(Path(a.label_map).read_text()) if a.label_map else {}

    # pre-registration: the procedure must be committed before a real run. A non-default config or a label map is part
    # of the procedure too, so it must be a committed, unmodified file in the repo.
    extra = [Path(f).resolve() for f in (a.config, a.label_map) if f]
    extra = [f for f in extra if f != (HERE / "config.yaml").resolve()]
    prov = fio.provenance([REPO / p for p in PROCEDURE_FILES if (REPO / p).exists()]
                          + [f for f in extra if REPO in f.parents])
    held_csv = repo_path(cfg["outputs"]["held_flat"])
    for f in [held_csv, *extra]:
        if not f.exists():
            continue
        if REPO not in f.parents:
            prov["git_uncommitted"].append(f"?? {f} (outside the repo)")
            continue
        rel = str(f.relative_to(REPO))
        if not fio._git("ls-files", "--", rel).strip():
            prov["git_uncommitted"].append(f"?? {rel}")
        prov["git_uncommitted"] += [ln for ln in fio._git("status", "--porcelain", "--", rel).splitlines() if ln]
    prov["config_path"] = str(Path(a.config).resolve())
    prov["config_sha256"] = hashlib.sha256(Path(a.config).read_bytes()).hexdigest()
    if run_mode == "e10.5" and prov["git_uncommitted"] and not a.allow_uncommitted and not a.labels_only:
        print("REFUSED: uncommitted changes to the files that define the readout (pre-registration):\n  "
              + "\n  ".join(prov["git_uncommitted"]) + "\nCommit them first, or pass --allow-uncommitted (recorded).",
              file=sys.stderr)
        return 3
    out_dir = Path(a.out) if a.out else repo_path(cfg["outputs"][{"smoke": "smoke", "rehearsal": "rehearsal",
                                                                  "e10.5": "results"}[run_mode]])
    record_path = out_dir / ("decision.json" if run_mode != "smoke" else "summary.json")
    default_record = repo_path(cfg["outputs"]["results"]) / "decision.json"
    if run_mode == "e10.5" and not a.labels_only:
        existing = [r for r in {record_path, default_record} if r.exists()]  # --out cannot sidestep the first record
        if existing and not a.overwrite:
            print(f"REFUSED: {existing[0]} exists (the first readout is the record); pass --overwrite to replace it.",
                  file=sys.stderr)
            return 6
        if record_path.exists():  # --overwrite keeps the replaced record next to the new one
            stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
            record_path.rename(record_path.with_name(f"decision.replaced-{stamp}.json"))

    # T1 statistics from the released files (committed code), cached under runs/
    cache_dir = repo_path(cfg["outputs"]["stats_cache"])
    seed = int(cfg["seed"])
    paths = {s: repo_path(p) for s, p in cfg["t1_files"].items()}
    st_by_stage, inputs = {}, {}
    for s, p in paths.items():
        st_by_stage[s], sha, how = stage_stats(p, s, cache_dir, cfg["label_key"], seed)
        inputs[s] = {"path": str(p.relative_to(REPO)) if REPO in p.parents else str(p), "sha256": sha, "stats": how}
    stg = cfg["stages"]
    if run_mode == "smoke":
        st = {"prev": None, "start": st_by_stage[stg["prev"]], "target": st_by_stage[stg["start"]]}
        stage_names = {"start": stg["prev"], "target": stg["start"]}
    else:
        tpath = paths[stg["start"]] if run_mode == "rehearsal" else e105
        tname = f"{stg['target']}-rehearsal" if run_mode == "rehearsal" else stg["target"]
        st_t, sha, how = stage_stats(tpath, tname, cache_dir, label_key, seed)
        inputs[tname] = {"path": str(tpath.relative_to(REPO)) if REPO in tpath.parents else str(tpath),
                         "sha256": sha, "stats": how}
        st = {"prev": st_by_stage[stg["prev"]], "start": st_by_stage[stg["start"]], "target": st_t}
        stage_names = {"prev": stg["prev"], "start": stg["start"], "target": tname}
    genes = st["start"]["genes"]
    for role, s in st.items():
        if s is not None:
            st[role], how = align_genes(s, genes)
            if how != "same order":
                inputs[stage_names[role]]["genes"] = how

    H = {role: dl.harmonise({k.split("|")[0] for k in s["keys"]}, cfg["harmonise"], override)
         for role, s in st.items() if s is not None}
    harmonisation = {"labels": {stage_names[r]: m for r, m in H.items()}, "override": override,
                     "override_sha256": hashlib.sha256(Path(a.label_map).read_bytes()).hexdigest()
                     if a.label_map else None, "label_key": label_key}
    for r, mapping in H.items():
        print(stage_names[r], json.dumps(mapping))
    if a.labels_only:
        return 0

    # held-flat genes: from the two released stages only, identical in every mode
    hf = dl.held_flat_genes(cfg["held_flat"], genes, st_by_stage)
    flat = dl.flat_mask(genes, hf.gene)
    held_note = "computed"
    if held_csv.exists():
        committed = pd.read_csv(held_csv)
        if sorted(committed.gene) != sorted(hf.gene):
            msg = f"held-flat list differs from {held_csv}: {sorted(set(committed.gene) ^ set(hf.gene))}"
            if run_mode == "e10.5":
                print("REFUSED: " + msg, file=sys.stderr)
                return 5
            print("WARNING: " + msg)
            held_note = f"differs from {held_csv.relative_to(REPO)}"
        else:
            held_note = f"matches {held_csv.relative_to(REPO)}"
    if run_mode == "smoke" or (run_mode == "e10.5" and not held_csv.exists()):
        write_table(hf.to_dict("records"), held_csv)
        held_note += f"; written to {held_csv.relative_to(REPO)}"
    print(f"held flat: {len(hf)} genes ({held_note}):", ", ".join(hf.gene))

    # external data (filter outputs only; registration required for a real run)
    reg = registered(cfg)
    if run_mode == "e10.5" and not reg["feng"]:
        print(f"REFUSED: the Feng source ({cfg['external']['register_keys']['feng']}) has no row in EXTERNAL_DATA.md "
              "(docs/DATA_POLICY.md: register before use).", file=sys.stderr)
        return 4
    jax_present = repo_path(cfg["external"]["jax_stats_dir"]).exists()
    use_jax = not a.no_jax and (reg["jax"] or run_mode != "e10.5") and jax_present
    ext, idx, ext_info = load_external(cfg, genes, use_jax, run_mode == "smoke")
    ext_info["registered_in_EXTERNAL_DATA"] = reg

    groups = build_groups(cfg, st, H, ext, idx, run_mode == "smoke")
    harmonisation["groups"] = {g: {"target": r["Lt"], "start": r["Ls"], "prev": r.get("Lp"),
                                   "external": sorted(r["ext"])} for g, r in groups.items()}
    none = np.zeros(len(genes), dtype=bool)
    type_rows = per_type_rows(cfg, groups, flat, "yes")
    pop_rows = population_rows(cfg, st, groups, flat, "yes")
    if run_mode == "smoke":  # the screened report's definitions: nothing held flat
        type_rows += per_type_rows(cfg, groups, none, "no")
        pop_rows += population_rows(cfg, st, groups, none, "no")
        if "jax" in ext:
            pop_rows += report_population(cfg, st["start"], st["target"], ext, idx, flat, "yes")
            pop_rows += report_population(cfg, st["start"], st["target"], ext, idx, none, "no")

    write_table(type_rows, out_dir / "per_type.csv", ROW_COLUMNS)
    write_table(pop_rows, out_dir / "population.csv")
    (out_dir / "harmonisation.json").write_text(json.dumps(clean(harmonisation), indent=1, default=json_default) + "\n")

    record = {"mode": run_mode, "stages": stage_names, "rule": cfg["decision"], "config_version": cfg["version"],
              "inputs": inputs, "external": ext_info, "held_flat": {"n": int(len(hf)), "genes": list(hf.gene),
                                                                     "note": held_note},
              "provenance": prov, "allow_uncommitted": bool(a.allow_uncommitted)}
    if run_mode == "smoke":
        record["verdict"] = "not evaluable: smoke geometry (no T1 stage before E8.5, so no momentum baseline)"
    else:
        readings = readings_for_decision(cfg, type_rows, pop_rows)
        record.update(dl.decide(cfg["decision"], readings))
        record["readings"] = readings
        if run_mode == "rehearsal":
            record["verdict"] = f"REHEARSAL (code-path check, not a result): {record['verdict']}"
    record_path.write_text(json.dumps(clean(record), indent=1, default=json_default) + "\n")
    print_summary(type_rows, pop_rows)
    print(f"\n[{run_mode}] verdict: {record['verdict']}")
    if run_mode != "smoke":
        print("checks:", json.dumps(clean(record["checks"]), default=json_default))
    print("outputs:", out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
