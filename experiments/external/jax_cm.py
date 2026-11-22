"""Row T71 (experiments/harness/TRIES.md): the cardiomyocyte direction of JAX/Qiu 2024 nuclei (register row X003),
recomputed from raw counts with a skeletal-muscle guard class.

The kill test of row T71 (experiments/t1_m1/cmjax_read.py) blends this direction into the CM group's entry of the
recomputed Feng directions. This script makes the direction and the counts its target-free pre-gate reads, and nothing
else, with the fixed settings of experiments/external/jax_cm.yaml.

  0. The gene universe is the gene table's genes (row T68's rule, `jax_blood.gene_universe`, imported): a count-matrix
     row id that meta/df_gene.csv does not list is removed when each matrix is loaded, before anything else. The
     record keeps the dropped count and, per embryo and per end, the share of the listed nuclei's UMIs they held
     (reported only).
  1. The groups are recomputed by row N54's marker-score machinery (feng_relabel's log1p_cp10k, unique_index,
     control_sets through jax_blood.controls, set_score, assign; none of it copied) with N54's sets as in
     feng_relabel.yaml and one competing class skm (Tnnt3, Tnni2, Mylpf, Myl1, Tnnc2; a set gene is never a control).
     The controls' gene means are taken per end and averaged with equal weight (`jax_blood.end_means`).
  2. Each r_cm nucleus with at least umi_target table-gene UMIs is downsampled to exactly umi_target
     (`jax_blood.downsample` through `jax_blood.embryo_mean`, default_rng([downsample_seed, embryo]) over the embryo's
     nuclei sorted by name). An embryo enters with at least min_embryo_nuclei such nuclei, an end with at least
     min_end_embryos such embryos (`end_direction`). Mean CP10k per embryo, equal-weight mean per end.
  3. log_ratio[g] = log((m_late + eps) / (m_early + eps)); direction = log_ratio x bracket_days / Delta, Delta the mean
     day label of the late end's contributing embryos minus the early end's (windows.parse_stage).
  Recorded per embryo for the pre-gate: the r_cm nuclei, those whose skm score is at least purity_score, the entering
  nuclei and the day label.

What it reads, and nothing else: the count matrices of the mode's embryos (each checked against the JAX folder's
MANIFEST.json); from meta/df_cell.csv the columns of CELL_WHITELIST (cell_id, embryo_id, day) through jax_blood's
reader; from meta/df_gene.csv gene_id and gene_short_name. JAX's own labels and stats/*.npz are never read.

Modes (ROW_EMBRYOS): practice, early end embryos 1-12 and late end 59-63; test, early end 13-22 and late end 59-63, run
only with --captain-ruling naming a docs/PLAN.md decision-log entry that names row T71. Any other list refuses, and
embryos 23-58 always refuse. Every nucleus read must carry one day label per embryo, outside T1's protected window,
and in practice mode outside the practice tasks' window (E8.5, E11.5]; in test mode the early end at or before E9.5 and
the late end after E13.5 (`check_days`, before any count is read).

It runs once per mode (`preflight`, all checked before the output folder is made and before any count is read): row T71
in TRIES.md; X003's line holds AMENDMENT_PHRASE (and in test mode TEST_PHRASE) and the three cell columns in
backticks; X010's line holds SKM_PHRASE and at least two of the five skm genes in backticks; the committed yaml with
sha256 YAML_SHA256 and status frozen, equal to its version at HEAD; a clean tree at --commit; the output folder absent
and inside the checkout's own data/; the inputs' format (`check_inputs`). It writes STARTED.json before the first count
matrix is opened, then the direction file, the per-embryo table and MANIFEST.json. It never removes anything; a folder
left by a stopped run blocks a second run of that mode.

Run (repo root, through experiments/t1_m1/cmjax.sh's `reduce` step; rdata from the external dependency group):
  uv run --frozen --no-sync python -I -B experiments/external/jax_cm.py --mode practice --commit <row commit>
"""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import sys
from collections.abc import Iterable, Mapping
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import feng_relabel as fr  # noqa: E402  (row N54's marker-score machinery and gate helpers)
import filterio as fio  # noqa: E402
import jax_blood as jb  # noqa: E402  (row T68's reducer: the rds loader, readers, gene rule, downsampling, means)
import windows  # noqa: E402

from vec import provenance  # noqa: E402

SOURCE = "jax_cm"
ROW = "T71"
CONFIG = fio.HERE / "jax_cm.yaml"
YAML_SHA256 = "3b618b06b129c6889ddbaa5a2178933ea7985e58028e576431a53cbba3e7e411"
TRIES = fio.REPO / "experiments" / "harness" / "TRIES.md"
REGISTER = fio.REPO / "EXTERNAL_DATA.md"
PLAN = fio.REPO / "docs" / "PLAN.md"
REGISTER_ROW = "X003"
SKM_REGISTER_ROW = "X010"
AMENDMENT_PHRASE = "T1 CM direction input, row T71"
TEST_PHRASE = "T1 CM test early end, embryos 13-22"
SKM_PHRASE = "skeletal-muscle competing class of row T71"
CODE = [Path(__file__), Path(jb.__file__), Path(fr.__file__), Path(fio.__file__), Path(windows.__file__),
        fr.PANEL_SOURCE]
STARTED = "STARTED.json"
MODES = ("practice", "test")
ENDS = jb.ENDS  # ("early", "late")
ROW_EMBRYOS = {"practice": {"early": tuple(range(1, 13)), "late": (59, 60, 61, 62, 63)},
               "test": {"early": tuple(range(13, 23)), "late": (59, 60, 61, 62, 63)}}
NEVER = tuple(range(23, 59))  # embryos no mode may list
CELL_WHITELIST = ("cell_id", "embryo_id", "day")
GENE_COLUMNS = jb.GENE_COLUMNS
AMENDMENT_COLUMNS = CELL_WHITELIST
GENE_UNIVERSE = jb.GENE_UNIVERSE
SKM_GENES = ("Tnnt3", "Tnni2", "Mylpf", "Myl1", "Tnnc2")
SKM_SOURCE = "X010_T71"
MIN_SKM_NAMED = 2
TASK = "T1"
PRACTICE_WINDOW = (8.5, 11.5)  # the practice tasks' window (E8.5, E11.5]: E8.5 exactly is the bracket's own end
TEST_EARLY_MAX, TEST_LATE_MIN = 9.5, 13.5  # test: early at or before E9.5, late after E13.5


class Refused(jb.Refused):
    """A rule of jax_cm.yaml or of this module stops the run (a jax_blood.Refused, so one except catches both)."""


# ------------------------------------------------------------------------------------------------ config


def load_config(path: str | Path = CONFIG) -> dict:
    cfg = yaml.safe_load(Path(path).read_text())
    check_config(cfg)
    return cfg


def check_config(cfg: Mapping) -> None:
    """Structural checks that need no data: the row, the gene rule, both modes' lists, the sets and the settings."""
    if cfg.get("row") != ROW:
        raise Refused(f"the yaml is not row {ROW}'s")
    if cfg.get("gene_universe") != GENE_UNIVERSE:
        raise Refused(f"the yaml's gene_universe is {cfg.get('gene_universe')!r}, not {GENE_UNIVERSE!r}")
    if set(cfg.get("embryos") or {}) != set(MODES):
        raise Refused(f"the yaml's embryo lists must name exactly the modes {MODES}")
    for mode in MODES:
        check_embryos(cfg["embryos"][mode], mode)
    sources = set(cfg.get("sources") or [])
    sets = cfg["sets"]
    for name, s in sets.items():
        if s.get("role") not in fr.ROLES:
            raise Refused(f"set {name}: role {s.get('role')!r} is not one of {fr.ROLES}")
        if s["role"] != "split" and not str(s.get("label", "")).startswith(fr.LABEL_PREFIX):
            raise Refused(f"set {name}: its label must start with {fr.LABEL_PREFIX!r}")
        bad = [g for g, src in s["genes"].items() if src not in sources]
        if bad:
            raise Refused(f"set {name}: genes {bad} have no listed source")
        if len(s["genes"]) < int(cfg["score"]["min_genes_found"]):
            raise Refused(f"set {name} has fewer genes than min_genes_found")
    skm = sets.get(cfg.get("skm_set"), {})
    if (skm.get("role") != "competing" or tuple(skm.get("genes") or ()) != SKM_GENES
            or set(skm["genes"].values()) != {SKM_SOURCE}):
        raise Refused(f"the skm set must be the competing class of {SKM_GENES}, sourced to {SKM_SOURCE}")
    if sets.get("cm", {}).get("label") != cfg.get("cm_label") or sets["cm"].get("role") != "target":
        raise Refused("the CM label is not the cm target set's label")
    for k in ("umi_target", "min_embryo_nuclei", "min_end_embryos"):
        if int(cfg[k]) < 1:
            raise Refused(f"{k} must be at least 1")
    for k in ("eps", "bracket_days"):
        if float(cfg[k]) <= 0:
            raise Refused(f"{k} must be positive")


def check_embryos(embryos: Mapping[str, Iterable[int]], mode: str) -> dict[str, tuple[int, ...]]:
    """The lists must be the mode's lists exactly (ROW_EMBRYOS); an embryo of 23-58 refuses whatever the mode."""
    if mode not in MODES:
        raise Refused(f"mode {mode!r} is not one of {MODES}")
    got = {e: tuple(sorted(int(x) for x in embryos.get(e, ()))) for e in ENDS}
    if any(n in NEVER for e in ENDS for n in got[e]):
        raise Refused("embryos 23-58 are never read by this row")
    if set(embryos) != set(ENDS) or got != ROW_EMBRYOS[mode]:
        raise Refused(f"the {mode} lists {dict(embryos)} are not row {ROW}'s {ROW_EMBRYOS[mode]}")
    return got


def end_of(n: int, mode: str) -> str:
    for e in ENDS:
        if n in ROW_EMBRYOS[mode][e]:
            return e
    raise Refused(f"embryo {n} is outside row {ROW}'s {mode} lists")


def scoring_config(cfg: Mapping) -> dict:
    return jb.scoring_config(cfg)


def output_dir(cfg: Mapping, mode: str) -> Path:
    return fio.REPO / cfg["output"][mode]


# ------------------------------------------------------------------------------------------------ the gate


def _register_line(row: str, path: Path | None) -> str | None:
    for line in Path(path or REGISTER).read_text(encoding="utf-8").splitlines():
        if line.startswith(f"| {row} |"):
            return " ".join(line.split())
    return None


def register_amended(mode: str = "practice", path: Path | None = None) -> bool:
    """X003's row names this use (AMENDMENT_PHRASE; in test mode also TEST_PHRASE) and the three cell columns read,
    each in backticks."""
    text = _register_line(REGISTER_ROW, path)
    if text is None:
        return False
    need = [AMENDMENT_PHRASE] + ([TEST_PHRASE] if mode == "test" else [])
    return all(p in text for p in need) and all(f"`{c}`" in text for c in AMENDMENT_COLUMNS)


def skm_registered(path: Path | None = None) -> bool:
    """X010's row holds the skm amendment (SKM_PHRASE) and names at least MIN_SKM_NAMED of SKM_GENES in backticks."""
    text = _register_line(SKM_REGISTER_ROW, path)
    return text is not None and SKM_PHRASE in text and sum(f"`{g}`" in text for g in SKM_GENES) >= MIN_SKM_NAMED


def input_paths(cfg: Mapping, src: Path, mode: str) -> dict[str, Path]:
    embryos = check_embryos(cfg["embryos"][mode], mode)
    inp = cfg["inputs"]
    wanted = sorted(n for e in ENDS for n in embryos[e])
    return {"cells": src / inp["cells"], "genes": src / inp["genes"],
            **{f"embryo_{n}": src / inp["counts"].format(n=n) for n in wanted}}


def check_inputs(cfg: Mapping, src: Path, mode: str, importer=importlib.import_module, verify=None) -> dict:
    """The inputs' format, reading no count and no cell row: the rds reader imports; every input of the mode's embryos
    is in place and matches the JAX folder's MANIFEST.json; the two csv header lines hold every column read."""
    verify = fio.verify_manifest_output if verify is None else verify
    try:
        importer(jb.RDS_READER)
    except ImportError as e:
        raise Refused(f"the rds reader {jb.RDS_READER} cannot be imported ({e})") from None
    paths = input_paths(cfg, Path(src), mode)
    for k, p in paths.items():
        if not p.is_file():
            raise Refused(f"{k}: {_rel(p)} is not in place")
        try:
            verify(p)
        except (ValueError, FileNotFoundError) as e:
            raise Refused(f"{k}: {e}") from None
    for k, need in (("cells", CELL_WHITELIST), ("genes", GENE_COLUMNS)):
        missing = [c for c in need if c not in jb.header_of(paths[k])]
        if missing:
            raise Refused(f"{paths[k].name}: columns {missing} are not in its header")
    return {"inputs_checked": sorted(paths)}


def preflight(cfg_path: Path, commit: str, dst: Path, mode: str = "practice", ruling: str | None = None,
              src: Path | None = None, importer=importlib.import_module, verify=None) -> dict:
    """Everything the run needs, before the output folder is made and before any count is read."""
    if mode not in MODES:
        raise Refused(f"mode {mode!r} is not one of {MODES}")
    if not fr.row_in_tries(ROW, TRIES):
        raise Refused(f"{ROW} has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if not register_amended(mode):
        raise Refused(f"{REGISTER_ROW}'s row in EXTERNAL_DATA.md does not name this use ({AMENDMENT_PHRASE!r}"
                      + (f", {TEST_PHRASE!r}" if mode == "test" else "") + ") and the columns read")
    if not skm_registered():
        raise Refused(f"{SKM_REGISTER_ROW}'s row does not hold the skm amendment ({SKM_PHRASE!r} and at least "
                      f"{MIN_SKM_NAMED} of {list(SKM_GENES)} in backticks)")
    if mode == "test" and not fr.ruling_logged(ruling, ROW, PLAN):
        raise Refused(f"test mode needs the captain's ruling: docs/PLAN.md has no decision-log entry {ruling!r} "
                      f"naming {ROW}")
    if Path(cfg_path).resolve() != CONFIG.resolve():
        raise Refused(f"only the committed {CONFIG} is run, not {cfg_path}")
    raw = Path(cfg_path).read_bytes()
    if fio.sha256_file(cfg_path) != YAML_SHA256:
        raise Refused(f"the yaml's sha256 is not the row's {YAML_SHA256}")
    cfg = yaml.safe_load(raw)
    if cfg.get("status") != "frozen":
        raise Refused(f"the yaml's status is {cfg.get('status')!r}: it is run only when 'frozen'")
    check_config(cfg)
    if fr.committed_at_head(CONFIG) != raw:
        raise Refused("the yaml differs from its committed version at HEAD")
    head, unclean = provenance.git_state(fio.REPO)
    if unclean is None or unclean:
        raise Refused(f"the tree is not clean ({unclean}): the reducer runs once, on a committed tree")
    if head != commit:
        raise Refused(f"HEAD is {head}, not the run commit {commit}")
    if Path(dst).exists():
        raise Refused(f"{dst} exists: the reducer never writes over a folder")
    check_inputs(cfg, fio.REPO / cfg["inputs"]["dir"] if src is None else Path(src), mode, importer=importer,
                 verify=verify)
    return {"row": ROW, "mode": mode, "ruling": ruling, "commit": head, "clean": True, "yaml_sha256": YAML_SHA256}


# ------------------------------------------------------------------------------------------------ the day labels


def embryo_days(col: Mapping[str, np.ndarray], mode: str) -> dict[int, str]:
    """Each listed embryo's one day label (an embryo with none or with two refuses), checked before any count: outside
    T1's window; practice: outside the practice tasks' window (E8.5, E11.5]; test: early at or before E9.5, late after
    E13.5."""
    out = {}
    for n in sorted(n for e in ENDS for n in ROW_EMBRYOS[mode][e]):
        labs = sorted(set(col["day"][np.char.equal(col["embryo_id"], f"embryo_{n}")].tolist()))
        if len(labs) != 1:
            raise Refused(f"embryo {n} carries {len(labs)} day labels, not one")
        day = windows.parse_stage(labs[0])
        if not windows.allowed(day, [TASK]):
            raise Refused(f"embryo {n}: day label {labs[0]!r} lies in T1's protected window")
        e = end_of(n, mode)
        if mode == "practice" and PRACTICE_WINDOW[0] < day <= PRACTICE_WINDOW[1]:
            raise Refused(f"embryo {n}: day label {labs[0]!r} lies in the practice tasks' window")
        if mode == "test" and ((e == "early" and day > TEST_EARLY_MAX) or (e == "late" and day <= TEST_LATE_MIN)):
            raise Refused(f"embryo {n}: day label {labs[0]!r} is not on the {e} side of the test's bracket")
        out[n] = labs[0]
    return out


# ------------------------------------------------------------------------------------------------ the arithmetic


def end_direction(per_embryo: Mapping[int, tuple[np.ndarray | None, int]], days: Mapping[int, str], mode: str,
                  min_nuclei: int, min_embryos: int, eps: float, bracket_days: float) -> dict:
    """Per end, the equal-weight mean over its embryos with at least min_nuclei entering nuclei, when at least
    min_embryos such embryos exist (else no direction); log_ratio = log((m_late + eps) / (m_early + eps)); Delta = the
    mean day of the late end's contributing embryos less the early end's; direction = log_ratio x bracket_days / Delta.
    """
    m, used = {}, {}
    for e in ENDS:
        used[e] = [n for n, (v, k) in sorted(per_embryo.items())
                   if end_of(n, mode) == e and v is not None and k >= min_nuclei]
        m[e] = np.mean([per_embryo[n][0] for n in used[e]], axis=0) if len(used[e]) >= min_embryos else None
    out = {"m_early": m["early"], "m_late": m["late"], "embryos_used": used, "log_ratio": None, "direction": None,
           "delta_days": None}
    if m["early"] is None or m["late"] is None:
        return out
    delta = (float(np.mean([windows.parse_stage(days[n]) for n in used["late"]]))
             - float(np.mean([windows.parse_stage(days[n]) for n in used["early"]])))
    if not delta > 0:
        raise Refused(f"the late end is not later than the early end (Delta {delta})")
    lr = np.log((m["late"] + eps) / (m["early"] + eps))
    return {**out, "log_ratio": lr, "direction": lr * (float(bracket_days) / delta), "delta_days": delta}


def scores_and_labels(L: sparse.csr_matrix, cfg: Mapping, set_idx: Mapping, ctrl: Mapping) -> tuple[dict, np.ndarray]:
    """Every nucleus's set scores (feng_relabel.set_score) and recomputed label (feng_relabel.assign)."""
    scores = {name: fr.set_score(L, set_idx[name], ctrl[name]) for name in cfg["sets"]}
    return scores, fr.assign(scores, scoring_config(cfg))


def _rel(p: Path) -> str:
    try:
        return str(Path(p).resolve().relative_to(fio.REPO))
    except ValueError:
        return str(p)


def reduce(cfg: Mapping, src: Path, gate: Mapping, loader=jb.load_rds, verify=fio.verify_manifest_output) -> dict:
    """The two passes over the mode's embryos (module docstring); returns the arrays and the record."""
    if not gate or gate.get("row") != ROW or gate.get("clean") is not True or gate.get("mode") not in MODES:
        raise Refused(f"no gate of row {ROW} on a clean tree: the reducer runs only through main(), after preflight")
    mode = gate["mode"]
    embryos = check_embryos(cfg["embryos"][mode], mode)
    wanted = sorted(n for e in ENDS for n in embryos[e])
    paths = input_paths(cfg, src, mode)
    inputs = {k: {"file": _rel(p), "sha256": verify(p).get("sha256")} for k, p in paths.items()}
    cells = jb.read_cells(paths["cells"], wanted, CELL_WHITELIST)
    col = cells["columns"]
    days = embryo_days(col, mode)  # before any count is read
    gene_ids, symbols = jb.read_genes(paths["genes"])
    table = jb.unique_table(gene_ids)
    cfg_s = scoring_config(cfg)
    store, info = {}, {}
    sums, counts = {}, {e: 0 for e in ENDS}
    umis = {e: [0.0, 0.0] for e in ENDS}  # per end: the listed nuclei's UMIs on all matrix genes, on dropped genes
    genes_first, universe = None, None
    for n in wanted:
        C, genes, names = loader(paths[f"embryo_{n}"])
        if genes_first is None:
            genes_first = list(genes)
            universe = jb.gene_universe(genes_first, table)
        elif list(genes) != genes_first:
            raise Refused(f"embryo {n}: its genes differ from the first embryo's")
        mine = np.flatnonzero(np.char.equal(col["embryo_id"], f"embryo_{n}"))
        listed = {str(c) for c in col["cell_id"][mine]}
        pos = [j for j, c in enumerate(names) if c in listed]
        C = sparse.csr_matrix(C)[pos]
        total_all = float(C.sum())
        C = sparse.csr_matrix(C[:, universe["kept"]])  # the dropped genes leave here, before anything else
        e = end_of(n, mode)
        umis[e][0] += total_all
        umis[e][1] += total_all - float(C.sum())
        nm = [names[j] for j in pos]
        order = np.argsort(np.asarray(nm, dtype=str), kind="stable")
        C, nm = C[order], [nm[j] for j in order]
        empty = np.asarray(C.sum(1)).ravel() <= 0
        C = C[~empty]
        L = fr.log1p_cp10k(C)
        sums[e] = sums.get(e, 0.0) + np.asarray(L.sum(0)).ravel()
        counts[e] += C.shape[0]
        store[n] = C
        info[n] = {"end": e, "day": days[n], "nuclei_in_matrix": len(names), "nuclei_listed": int(mine.size),
                   "nuclei_read": int(C.shape[0]), "empty_dropped": int(empty.sum()),
                   "left_out_umi_share": (total_all - float(C.sum())) / total_all if total_all > 0 else None}
    kept_ids = [genes_first[j] for j in universe["kept"]]
    jsym = [symbols[table[g]] for g in kept_ids]
    found, set_idx, uidx = jb.set_indices(cfg, jsym)
    means = jb.end_means(sums, counts)
    ctrl = jb.controls(cfg, jsym, means, set_idx, uidx)
    if any(jsym[i] in SKM_GENES for c in ctrl.values() for i in c):
        raise Refused("an skm gene was drawn as a control gene")
    per_embryo, labels_count = {}, {}
    for n in wanted:
        C = store.pop(n)
        scores, lab = scores_and_labels(fr.log1p_cp10k(C), cfg, set_idx, ctrl)
        cm = lab == cfg["cm_label"]
        high = cm & (scores[cfg["skm_set"]] >= float(cfg["purity_score"]))
        labels_count[n] = {x: int(np.sum(lab == x)) for x in fr.all_labels(cfg_s)}
        v, k = jb.embryo_mean(C[np.flatnonzero(cm)], int(cfg["umi_target"]), [int(cfg["downsample_seed"]), n])
        per_embryo[n] = (v, k)
        info[n].update(cm_assigned=int(cm.sum()), cm_skm_high=int(high.sum()), cm_entering=k)
    res = end_direction(per_embryo, days, mode, int(cfg["min_embryo_nuclei"]), int(cfg["min_end_embryos"]),
                        float(cfg["eps"]), float(cfg["bracket_days"]))
    gnames = np.asarray(jsym, dtype=str)
    nan = np.full(len(gnames), np.nan)
    record = {**universe["record"], "mode": mode, "embryos_read": wanted,
              "embryo_lists": {e: list(embryos[e]) for e in ENDS}, "columns_read": list(CELL_WHITELIST),
              "gene_columns_read": list(GENE_COLUMNS), "authors_labels_read": False, "stats_npz_read": False,
              "inputs": inputs, "rows_other_embryos_dropped": cells["rows_other_embryos"], "per_embryo": info,
              "end_left_out_umi_share": {e: (u[1] / u[0] if u[0] > 0 else None) for e, u in umis.items()},
              "labels_per_embryo": labels_count, "embryos_used": res["embryos_used"],
              "delta_days": res["delta_days"], "bracket_days": float(cfg["bracket_days"]),
              "sets": {k: {"found": found[k], "controls": [str(gnames[i]) for i in ctrl[k]]} for k in cfg["sets"]},
              "has_direction": res["direction"] is not None}
    arrays = {"gene_id": np.asarray(kept_ids, dtype=str), "symbol": gnames,
              "m_early": nan if res["m_early"] is None else res["m_early"],
              "m_late": nan if res["m_late"] is None else res["m_late"],
              "log_ratio": nan if res["log_ratio"] is None else res["log_ratio"],
              "direction": nan if res["direction"] is None else res["direction"],
              "embryo": np.array(sorted(per_embryo), dtype=np.int64),
              "embryo_mean_cp10k": np.vstack([per_embryo[n][0] if per_embryo[n][0] is not None else nan
                                              for n in sorted(per_embryo)]),
              "embryo_nuclei": np.array([per_embryo[n][1] for n in sorted(per_embryo)], dtype=np.int64)}
    return {"arrays": arrays, "record": record}


def write_outputs(dst: Path, cfg: Mapping, cfg_path: Path, gate: Mapping, res: Mapping) -> dict:
    out = cfg["output"]
    fio.save_npz_atomic(dst / out["direction"], **res["arrays"])
    tmp = dst / f".{out['embryos']}.partial"
    with open(tmp, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["embryo", "end", "day", "nuclei_read", "cm_assigned", "cm_skm_high", "cm_entering"])
        for n, r in sorted(res["record"]["per_embryo"].items()):
            w.writerow([n, r["end"], r["day"], r["nuclei_read"], r["cm_assigned"], r["cm_skm_high"],
                        r["cm_entering"]])
    tmp.replace(dst / out["embryos"])
    manifest = {"source": SOURCE, "row": ROW, "register_row": REGISTER_ROW,
                "config": {"file": _rel(cfg_path), "sha256": fio.sha256_file(cfg_path), "status": cfg.get("status")},
                "gate": dict(gate), **res["record"],
                "settings": {k: cfg[k] for k in ("umi_target", "downsample_seed", "min_embryo_nuclei",
                                                 "min_end_embryos", "eps", "bracket_days", "cm_label", "skm_set",
                                                 "purity_score")},
                "outputs": {out["direction"]: {"sha256": fio.sha256_file(dst / out["direction"])},
                            out["embryos"]: {"sha256": fio.sha256_file(dst / out["embryos"])}},
                "provenance": fio.provenance(CODE)}
    fio.write_manifest(dst, manifest)
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--commit", required=True, help="the run commit (full hash): HEAD must be it")
    ap.add_argument("--mode", required=True, choices=list(MODES))
    ap.add_argument("--captain-ruling", default=None, help="test mode only: the decision-log entry naming the row")
    a = ap.parse_args(argv)
    try:  # no count and no cell row of JAX's is read in this block (check_inputs reads two header lines)
        cfg = load_config(CONFIG)
        dst = fr.output_dir_ok(output_dir(cfg, a.mode))
        gate = preflight(CONFIG, a.commit, dst, a.mode, a.captain_ruling)
    except (Refused, fr.Refused, jb.Refused) as e:
        raise SystemExit(f"refused: {e}") from None
    dst.mkdir(parents=True)
    (dst / STARTED).write_text(json.dumps({"gate": gate, "config_sha256": YAML_SHA256}, indent=1) + "\n")
    res = reduce(cfg, fio.REPO / cfg["inputs"]["dir"], gate)
    write_outputs(dst, cfg, CONFIG, gate, res)
    print(f"written {_rel(dst)}")  # no count
    return 0


if __name__ == "__main__":
    sys.exit(main())
