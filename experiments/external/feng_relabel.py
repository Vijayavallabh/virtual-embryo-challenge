"""Feng et al. 2022 (GEO GSE193346, register row X002): cell groups recomputed from permitted cells alone.

An organiser asked on 2026-11-16 that annotations learned using protected cells be left out and recomputed on permitted
data only. Feng's own labels (the metadata column cell_type_spec) were made with every stage pooled, the in-window
stages included, so this script labels the permitted cells again from their raw counts, with the fixed marker scores
of experiments/external/feng_relabel.yaml (option (c) of the design of 2026-11-19), and writes a second statistics
file in the clean reduction's layout with the recomputed label in place of the authors'.

What it reads, and nothing else:
  - cd1_legal_counts.npz and cd1_legal_genes.txt of the clean reduction (experiments/external/feng2022_geo.py), each
    checked against that reduction's MANIFEST.json;
  - from cd1_legal_obs.csv, the cell names (first column) and the columns of OBS_WHITELIST (stage, zone, strain), by a
    reader that refuses any other column and refuses a whitelist naming an annotation, cluster or embedding column. The
    authors' cell_type_spec is never read; nothing here can read it (tests/test_feng_relabel.py).
What it writes (the yaml's output.dir, data/derived/feng_relabel_c: a folder inside the checkout's own data/
directory, which is a real directory in every checkout; never under data/external or data/raw, which a campaign
checkout on the DGX links into the reference checkout; the gate refuses an output folder that resolves anywhere else),
with its own MANIFEST.json: cd1_stats_relabel.npz (keys 'stage|zone|<recomputed label>|half', the fields of
feng2022_geo.group_stats, the same seeded halves), cd1_relabel_cells.csv (per cell: name, stage, zone, label, every
score) and the manifest (inputs' and outputs' sha256, the yaml's sha256, the code's commit and sha256, the per-stage
unassigned shares, the cap's verdict, the groups' cells at each end, the fallbacks the rule fires, the set genes found
and their controls).

It runs once, and only in order (`preflight`, all checked before any Feng file is read): the committed
experiments/external/feng_relabel.yaml and no other config; its status "frozen"; no gene UNSOURCED and no choice
PENDING in it; the yaml as committed at HEAD; row N55's dry run committed at HEAD, with the yaml's min_score and
min_margin the ones N55's pre-stated rule gave and every other value as N55 read it (`check_n55`); a clean tree; ledger
row N54 in experiments/harness/TRIES.md; the captain's ruling, a docs/PLAN.md decision-log entry named by --ruling,
naming N54; and no relabelling recorded or started before (neither its folder under data/ nor the committed copy of
its manifest, RECORD, at HEAD or in the tree). The manifest records that gate, and every reader of the file refuses one
written without it (feng_dirs.check_recomputed). After the run it writes a copy of the manifest to RECORD, which the
main session commits at once: N54's listing and the builder refuse until the copy is committed and equal to the
manifest, so a second run cannot go unseen.

A run that stops after preflight has passed (a set short of score.min_genes_found in Feng's unique symbols, a set that
draws no control gene, a cell with no counts, a strain or shape mismatch, any other error, an interrupt) takes the
yaml's terminal.on_stop route, "none": main() writes a stop record (the manifest and RECORD, with the reason, usable
false and route none), prints the none command and exits with STOP_EXIT. Nothing in the yaml is edited and nothing is
rerun under row N54. A marker (STARTED.json) is written in the output folder before the first Feng file is opened, so a
killed run also blocks a second one. Run (DGX campaign checkout, repo root; minutes):
  uv run --frozen --no-sync python -I -B experiments/external/feng_relabel.py --ruling "<decision-log id>"
"""

from __future__ import annotations

import argparse
import ast
import copy
import csv
import hashlib
import json
import re
import subprocess
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import feng2022_geo  # noqa: E402
import filterio as fio  # noqa: E402

SOURCE = "feng_relabel"
CONFIG = fio.HERE / "feng_relabel.yaml"
DIRECTION_CONFIG = fio.REPO / "experiments" / "t1_direction" / "config.yaml"
# The team's released-data marker panel (experiments/t1_scope/label_evidence.py MARKERS, committed 2026-10-19 at
# b875016): its every gene is kept out of every control set (the yaml's score rule); read by literal parse, not import.
PANEL_SOURCE = fio.REPO / "experiments" / "t1_scope" / "label_evidence.py"
CODE = [Path(__file__), Path(feng2022_geo.__file__), Path(fio.__file__), PANEL_SOURCE]
TRIES = fio.REPO / "experiments" / "harness" / "TRIES.md"
PLAN = fio.REPO / "docs" / "PLAN.md"
ROW = "N54"
# The committed copy of the relabelling's manifest: written after the run (or after its stop), committed by the main
# session at once.
RECORD = fio.HERE / "results" / "feng_relabel" / "MANIFEST.json"
# Row N55's dry run on the released E9.5 (feng_relabel_dryrun.py), committed before the freeze; the gate reads it.
N55_OUT = fio.HERE / "results" / "feng_relabel_dryrun" / "n55.json"
STARTED = "STARTED.json"  # written in the output folder before the first Feng file is opened
STOP_EXIT = 3  # main()'s exit code when a run stopped after preflight and took route none
PENDING = "PENDING"  # a choice the captain's ruling fixes before the yaml is frozen
ROW_ID = re.compile(r"[A-Z][0-9]+[a-z]?")  # a TRIES row's ID (n12_hygiene_size.ROW_ID)
CP10K = 1e4
UNSOURCED = "UNSOURCED"
# What a literature source under the yaml's `sources` must state (besides a quoted sentence per gene under `names`)
LITERATURE_KEYS = ("citation", "doi", "samples", "licence", "names")
# The only metadata columns this module may read (besides the cell names, the first column). Fixed in code, not only
# in the yaml, so that no config can widen it.
OBS_WHITELIST = ("stage", "zone", "strain")
# A column whose name matches this is an annotation, cluster, embedding or fitted value: it can never be whitelisted.
ANNOTATION_PATTERN = re.compile(r"cell_?type|celltype|type_spec|cluster|annot|label|ident|umap|tsne|t_sne|pca|"
                                r"harmony|embed|seurat|snn|leiden|louvain|^res\.|lineage|subtype|class|score|phase",
                                re.IGNORECASE)
LABEL_PREFIX = "r_"  # every recomputed label starts with this; no author label does (readers check it)
ROLES = ("target", "competing", "split")
ROUTES = ("recomputed", "none")  # the two routes of the yaml's `terminal` rule (build_test.py t1 --feng-groups)
# The one base the terminal rule allows (the main session's decision (1) of 2026-11-19, under its delegation):
# library 3's cells at library 3's own mix, the first upload's design. The mid base and its 0.60 are in no later file
# (docs/S2_RULES.md, the amendment of 2026-11-19 for T1), so they are not offered.
BASE_ARGS = {"lib3": ["--base", "lib3"]}
# The first upload's command, the one design either route builds (equal to build_test.FENG_NONE_COMMAND; tested)
FIRST_UPLOAD_COMMAND = "t1 --g-mom --base lib3"


class Refused(RuntimeError):
    """A rule of feng_relabel.yaml or of this module stops the run."""


# ------------------------------------------------------------------------------------------------ config


def file_sha256(path: str | Path) -> str:
    return fio.sha256_file(path)


def load_config(path: str | Path = CONFIG) -> dict:
    cfg = yaml.safe_load(Path(path).read_text())
    check_config(cfg)
    return cfg


def check_config(cfg: Mapping) -> None:
    """Structural checks that need no data: whitelist, roles, split, labels, groups and fallbacks."""
    check_whitelist(cfg["obs_columns"])
    sets = cfg["sets"]
    for name, s in sets.items():
        if s.get("role") not in ROLES:
            raise Refused(f"set {name}: role {s.get('role')!r} is not one of {ROLES}")
        if not s.get("genes"):
            raise Refused(f"set {name} has no genes")
        if s["role"] != "split" and not str(s.get("label", "")).startswith(LABEL_PREFIX):
            raise Refused(f"set {name}: its label must start with {LABEL_PREFIX!r}")
        for g, src in s["genes"].items():
            if src == UNSOURCED:
                if s["role"] != "competing":
                    raise Refused(f"set {name}, gene {g}: only a competing class may hold an UNSOURCED gene")
                continue
            for sid in gene_sources(src):
                if sid == UNSOURCED or sid not in cfg["sources"]:
                    raise Refused(f"set {name}, gene {g}: source {sid!r} is not under `sources`")
                ref = cfg["sources"][sid]
                if isinstance(ref, Mapping):  # a literature source quotes the sentence that names each of its genes
                    missing = [k for k in LITERATURE_KEYS if not ref.get(k)]
                    if missing:
                        raise Refused(f"source {sid} lacks {missing}")
                    if not str((ref.get("names") or {}).get(g, "")).strip():
                        raise Refused(f"set {name}, gene {g}: source {sid} quotes no sentence naming it")
        both = set(s["genes"]) & set(s.get("excluded") or {})
        if both:
            raise Refused(f"set {name}: genes {sorted(both)} are both members and excluded")
    dropped = set(cfg.get("dropped") or {})
    kept = sorted(dropped & {g for s in sets.values() for g in s["genes"]})
    if kept:
        raise Refused(f"genes {kept} are both dropped and members of a set")
    sp =cfg["assign"]["ec_split"]
    if sets.get(sp["class"], {}).get("role") != "target":
        raise Refused(f"the split class {sp['class']!r} is not a target set")
    for k in ("plus", "minus"):
        if sets.get(sp[k], {}).get("role") != "split":
            raise Refused(f"the split's {k} set {sp[k]!r} is not a split set")
    labels = set(all_labels(cfg))
    for group, labs in cfg["groups"].items():
        bad = [lab for lab in labs if lab not in labels or lab == cfg["assign"]["unassigned"]]
        if bad:
            raise Refused(f"group {group}: labels {bad} are not recomputed class labels")
        if any(sets[n]["role"] == "competing" and sets[n]["label"] in labs for n in sets if sets[n]["role"] != "split"):
            raise Refused(f"group {group} holds a competing class")
    for group, fb in cfg["fallback"].items():
        if group not in cfg["groups"] or (fb != "none" and fb not in cfg["groups"]):
            raise Refused(f"fallback {group} -> {fb}: both must be groups (or 'none')")
    if set(cfg["fallback"]) != set(cfg["groups"]):
        raise Refused("every group needs a fallback (or 'none')")
    if not 0.0 < float(cfg["unassigned_cap"]) < 1.0:
        raise Refused("unassigned_cap must lie in (0, 1)")
    t = cfg["terminal"]
    if t["otherwise"] != "none":
        raise Refused("terminal.otherwise must be 'none' (the route when the recomputed groups fail)")
    if not set(t["requires_direction"]) <= set(cfg["groups"]):
        raise Refused("terminal.requires_direction must name groups of the map")
    if t["none_with_g_mom"] not in (True, PENDING):  # main's none file is built with g_mom on (decision (2))
        raise Refused(f"terminal.none_with_g_mom must be true (main's none file has g_mom on) or {PENDING}")
    if t["base"] not in (*BASE_ARGS, PENDING):
        raise Refused(f"terminal.base must be one of {sorted(BASE_ARGS)} or {PENDING}")
    if t.get("on_stop") != "none":
        raise Refused("terminal.on_stop must be 'none' (a run that stops after preflight takes route none)")
    check_n55_rule(cfg)


def check_n55_rule(cfg: Mapping) -> None:
    """Row N55's pre-stated rule (the yaml's `n55_rule`): an ordered list of (min_score, min_margin) settings, and the
    chamber labels it checks, each with the released labels it pools and the recomputed labels that count as its own."""
    r = cfg.get("n55_rule")
    if not isinstance(r, Mapping) or not r.get("order") or not r.get("chamber_labels"):
        raise Refused("n55_rule needs an ordered list of settings and the chamber labels it checks")
    for s in r["order"]:
        if set(s) != {"min_score", "min_margin"}:
            raise Refused(f"n55_rule.order: a setting names min_score and min_margin only, not {sorted(s)}")
    labels = set(all_labels(cfg))
    for name, spec in r["chamber_labels"].items():
        if not spec.get("released") or not spec.get("own") or not set(spec["own"]) <= labels:
            raise Refused(f"n55_rule.chamber_labels.{name}: released labels and own recomputed labels are needed")


def gene_sources(src) -> list[str]:
    """A gene's source ids: the yaml gives one id or a non-empty list of ids."""
    if isinstance(src, str):
        return [src]
    if isinstance(src, list) and src and all(isinstance(s, str) for s in src):
        return list(src)
    raise Refused(f"a gene's source must be an id or a list of ids, not {src!r}")


def all_labels(cfg: Mapping) -> list[str]:
    """Every label the assignment can give, in a fixed order."""
    a = cfg["assign"]
    sp = a["ec_split"]
    out = []
    for name, s in cfg["sets"].items():
        if s["role"] == "split":
            continue
        if name == sp["class"]:
            out += [sp["plus_label"], sp["minus_label"], sp["ambiguous_label"]]
        else:
            out.append(s["label"])
    return [*out, a["unassigned"]]


def classes(cfg: Mapping) -> list[str]:
    """The sets a cell can be assigned to (targets and competing classes), in the yaml's order."""
    return [n for n, s in cfg["sets"].items() if s["role"] in ("target", "competing")]


def unsourced_genes(cfg: Mapping) -> dict[str, list[str]]:
    return {n: [g for g, src in s["genes"].items() if src == UNSOURCED] for n, s in cfg["sets"].items()
            if any(src == UNSOURCED for src in s["genes"].values())}


def pending_choices(obj, where: str = "") -> list[str]:
    """Every place in the yaml that still holds PENDING (a choice the captain's ruling fixes before the freeze)."""
    if isinstance(obj, Mapping):
        return [p for k, v in obj.items() for p in pending_choices(v, f"{where}.{k}" if where else str(k))]
    if isinstance(obj, list):
        return [p for i, v in enumerate(obj) for p in pending_choices(v, f"{where}[{i}]")]
    return [where] if obj == PENDING else []


def without_unsourced(cfg: Mapping) -> dict:
    """The yaml's pre-set deletion rule (`unsourced_rule`), as the main session applies it by hand before the freeze to
    a gene still UNSOURCED: the gene is deleted, and a competing class left with fewer than min_genes_found genes is
    deleted whole. The tests use it to show the result is a valid configuration; the run never applies it."""
    out = copy.deepcopy(dict(cfg))
    keep = {}
    for name, s in out["sets"].items():
        genes = {g: src for g, src in s["genes"].items() if src != UNSOURCED}
        if len(genes) >= int(out["score"]["min_genes_found"]):
            keep[name] = {**s, "genes": genes}
        elif s["role"] != "competing":
            raise Refused(f"set {name} would fall under min_genes_found: a target or split set cannot be deleted")
    out["sets"] = keep
    check_config(out)
    return out


def check_min_cells(cfg: Mapping, direction_cfg: Mapping) -> None:
    """The yaml's minimum must be the readout's (experiments/t1_direction/config.yaml), which feng_dirs applies."""
    if int(cfg["min_cells_ext"]) != int(direction_cfg["metrics"]["min_cells_ext"]):
        raise Refused(f"min_cells_ext {cfg['min_cells_ext']} differs from the readout's "
                      f"{direction_cfg['metrics']['min_cells_ext']}")


def panel_genes(path: str | Path = PANEL_SOURCE) -> set[str]:
    """Every gene of the released-data panel (label_evidence.MARKERS), read by literal parse of its assignment."""
    tree = ast.parse(Path(path).read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "MARKERS" for t in node.targets):
            return {str(g) for genes in ast.literal_eval(node.value).values() for g in genes}
    raise Refused(f"{path} has no MARKERS assignment")


def named_genes(cfg: Mapping) -> set[str]:
    """Every gene the yaml names: set members (UNSOURCED included), the genes a set leaves out by name, and the genes
    taken out of the sets (`dropped`)."""
    return ({g for s in cfg["sets"].values() for g in [*s["genes"], *(s.get("excluded") or {})]}
            | set(cfg.get("dropped") or {}))


def route(usable: bool, fallbacks: Mapping[str, Mapping], cfg: Mapping) -> str:
    """The yaml's terminal rule, fixed before any run: "recomputed" when the file passed the cap and every group of
    terminal.requires_direction has a direction after the fallback; otherwise terminal.otherwise ("none": no Feng
    direction on any group). No other outcome exists, and none is chosen after the result."""
    need = cfg["terminal"]["requires_direction"]
    if usable and all(fallbacks.get(g, {}).get("takes") is not None for g in need):
        return "recomputed"
    return cfg["terminal"]["otherwise"]


def build_command(route_: str, cfg: Mapping) -> list[str]:
    """The one builder command the terminal rule fixes for `route_` (arguments after build_test.py, without --out):
    the first upload's command (FIRST_UPLOAD_COMMAND, t1 --g-mom --base lib3) with --feng-groups and the route. Route
    "none" is main's none file as built (g_mom on); route "recomputed" is the same design on the recomputed groups.
    Refused while a choice it needs is PENDING."""
    t = cfg["terminal"]
    if route_ not in ROUTES:
        raise Refused(f"route must be one of {ROUTES}")
    if t["base"] == PENDING or t["none_with_g_mom"] == PENDING:
        raise Refused(f"the terminal rule's choices are still {PENDING}: they are fixed before the freeze")
    return ["t1", "--g-mom", *BASE_ARGS[t["base"]], "--feng-groups", route_]


# ------------------------------------------------------------------------------------------------ row N55's rule


def chamber_check(released: np.ndarray, labels: np.ndarray, chamber_labels: Mapping) -> dict:
    """Per chamber label of n55_rule (its released labels pooled): its cells, the share given its own recomputed labels
    (pooled), and the largest share given any one other recomputed label (r_unassigned included). It passes when its
    cells go most often to its own labels: the own share is above every other label's share. No cells: it fails."""
    released, labels = np.asarray(released, dtype=str), np.asarray(labels, dtype=str)
    out = {}
    for name, spec in chamber_labels.items():
        m = np.isin(released, list(spec["released"]))
        n = int(m.sum())
        if n == 0:
            out[name] = {"n": 0, "own_share": None, "largest_other": None, "passed": False}
            continue
        lab = labels[m]
        own = float(np.mean(np.isin(lab, list(spec["own"]))))
        others = [float(np.mean(lab == r)) for r in sorted(set(lab.tolist()) - set(spec["own"]))]
        largest = max(others, default=0.0)
        out[name] = {"n": n, "own_share": own, "largest_other": largest, "passed": bool(own > largest)}
    return out


def n55_outcome(released: np.ndarray, scores: Mapping[str, np.ndarray], cfg: Mapping) -> dict:
    """Row N55's pre-stated rule on the released E9.5, applied by code: for each setting of n55_rule.order, in order,
    the labels, the unassigned share, the cap's verdict and the chamber checks; a setting passes when the share is at
    most unassigned_cap and every chamber label passes. (i) The first setting passes: kept. (ii) Otherwise the first
    later setting that passes: changed to it. (iii) None passes: the first is kept, and the guard failed on the
    released dry run (N54's first line says so). The outcome's setting is what the yaml's min_score and min_margin
    must be at the freeze (check_n55); nothing else can move by N55."""
    rule = cfg["n55_rule"]
    cap, un = float(cfg["unassigned_cap"]), cfg["assign"]["unassigned"]
    rows = []
    for s in rule["order"]:
        lab = assign(scores, cfg, s["min_score"], s["min_margin"])
        share = float(np.mean(lab == un))
        chambers = chamber_check(released, lab, rule["chamber_labels"])
        passed = bool(share <= cap and all(c["passed"] for c in chambers.values()))
        rows.append({"setting": dict(s), "unassigned_share": share, "cap_passed": share <= cap,
                     "chambers": chambers, "passed": passed})
    if rows[0]["passed"]:
        outcome, setting = "kept", rows[0]["setting"]
    else:
        first = next((r for r in rows[1:] if r["passed"]), None)
        outcome, setting = ("changed", first["setting"]) if first else ("guard_failed", rows[0]["setting"])
    return {"settings": rows, "outcome": outcome, "setting": setting}


# ------------------------------------------------------------------------------------------------ the gate


def row_in_tries(row: str, path: Path | None = None) -> bool:
    """`row` is a row ID of experiments/harness/TRIES.md's tables (each table row's first word; n12's tries_rows)."""
    for ln in Path(path or TRIES).read_text(encoding="utf-8").splitlines():
        if ln[:2] == "| ":
            w = ln[2:].split("|", 1)[0].split()
            if w and ROW_ID.fullmatch(w[0]) and w[0] == row:
                return True
    return False


def ruling_logged(ruling: str | None, row: str = ROW, plan: Path | None = None) -> bool:
    """The docs/PLAN.md decision-log entry `ruling` (e.g. '2026-11-20 (a)') exists and names `row` anywhere in it: the
    entry runs from its '- <ruling>' line (at the margin) to the next entry or heading."""
    if ruling is None or not re.fullmatch(r"\d{4}-\d{2}-\d{2} \([a-z]+\)", ruling):
        return False
    block: list[str] | None = None
    for line in Path(plan or PLAN).read_text(encoding="utf-8").splitlines():
        if block is not None:
            if line.startswith("- ") or line.startswith("#"):
                break
            block.append(line)
        elif line.startswith(f"- {ruling}"):
            block = [line]
    return block is not None and any(re.search(rf"\b{re.escape(row)}\b", ln) for ln in block)


def committed_at_head(path: Path) -> bytes | None:
    """The bytes of `path` as committed at HEAD, or None if HEAD has no such file."""
    rel = Path(path).resolve().relative_to(fio.REPO).as_posix()
    r = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=fio.REPO, capture_output=True)
    return r.stdout if r.returncode == 0 else None


def yaml_at_commit(commit: str, path: Path | None = None) -> bytes | None:
    """The bytes of the yaml (default CONFIG) as committed at `commit`, or None."""
    rel = Path(CONFIG if path is None else path).resolve().relative_to(fio.REPO).as_posix()
    r = subprocess.run(["git", "show", f"{commit}:{rel}"], cwd=fio.REPO, capture_output=True)
    return r.stdout if r.returncode == 0 else None


# What N55 may move in the yaml (by its rule), and what the freeze sets: nothing else may differ from the yaml N55 read
N55_MOVABLE = (("status",), ("assign", "min_score"), ("assign", "min_margin"))


def _without(cfg: Mapping, paths) -> dict:
    out = copy.deepcopy(dict(cfg))
    for p in paths:
        d = out
        for k in p[:-1]:
            d = d.get(k, {})
        if isinstance(d, dict):
            d.pop(p[-1], None)
    return out


def check_n55(cfg: Mapping) -> dict:
    """Row N55 before the run (decision (4) of 2026-11-19): its output committed at HEAD; its pre-stated rule's
    setting equal to the yaml's min_score and min_margin; and the yaml N55 read (at N55's stamp commit, with the sha256
    N55 recorded) equal to this one in every value but the status and those two. Returns N55's outcome."""
    raw = committed_at_head(N55_OUT)
    if raw is None:
        raise Refused(f"row N55's dry run is not committed at HEAD ({_rel(N55_OUT)}): it runs, and "
                      "its rule is applied, before the freeze")
    rep = json.loads(raw)
    rule = rep.get("rule") or {}
    want = rule.get("setting") or {}
    have = {"min_score": cfg["assign"]["min_score"], "min_margin": cfg["assign"]["min_margin"]}
    if {k: float(v) for k, v in want.items()} != {k: float(v) for k, v in have.items()}:
        raise Refused(f"the yaml's thresholds {have} are not the setting row N55's rule gave ({want}, "
                      f"{rule.get('outcome')})")
    read = yaml_at_commit(str(rep.get("stamp", {}).get("commit")))
    if read is None or hashlib.sha256(read).hexdigest() != rep.get("inputs", {}).get("config_sha256"):
        raise Refused("the yaml row N55 read cannot be found at N55's stamp commit with the sha256 it recorded")
    if _without(yaml.safe_load(read), N55_MOVABLE) != _without(cfg, N55_MOVABLE):
        raise Refused("the yaml differs from the one row N55 read in more than its status, min_score and min_margin")
    return {"outcome": rule.get("outcome"), "setting": want, "commit": rep["stamp"]["commit"]}


def output_dir_ok(dst: Path, repo: Path | None = None) -> Path:
    """The output folder resolves inside the checkout's own data/ directory, which must be a real directory (not a
    link): so a campaign checkout whose data/external and data/raw link into the reference checkout never writes
    there. Returns the resolved folder."""
    data = Path(fio.REPO if repo is None else repo) / "data"
    if data.is_symlink():
        raise Refused(f"{data} is a link: the output folder must be inside the checkout's own data directory")
    real = data.resolve()
    out = Path(dst).resolve()
    if real not in out.parents:
        raise Refused(f"the output folder must be inside {real} (the checkout's own data directory); got {out}")
    return out


def preflight(cfg_path: Path, ruling: str | None, dst: Path) -> dict:
    """Everything the run needs, checked before any Feng file is opened (the module docstring's order). Returns the
    gate the manifest records."""
    from vec import provenance

    if Path(cfg_path).resolve() != CONFIG.resolve():
        raise Refused(f"only the committed {CONFIG} is run, not {cfg_path}")
    raw = Path(cfg_path).read_bytes()
    cfg = yaml.safe_load(raw)
    if cfg.get("status") != "frozen":
        raise Refused(f"the yaml's status is {cfg.get('status')!r}: it is run only when 'frozen'")
    if unsourced_genes(cfg):
        raise Refused(f"the yaml holds UNSOURCED genes {unsourced_genes(cfg)}")
    if pending_choices(cfg):
        raise Refused(f"the yaml holds choices still {PENDING}: {pending_choices(cfg)}")
    n55 = check_n55(cfg)
    if not row_in_tries(ROW):
        raise Refused(f"{ROW} has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if not ruling_logged(ruling):
        raise Refused(f"docs/PLAN.md has no decision-log entry {ruling!r} naming {ROW}: the run follows the ruling")
    if committed_at_head(CONFIG) != raw:
        raise Refused("the yaml differs from its committed version at HEAD")
    commit, unclean = provenance.git_state(fio.REPO)
    if unclean is None or unclean:
        raise Refused(f"the tree is not clean ({unclean}): the relabelling runs once, on a committed tree")
    if Path(dst).exists() or RECORD.exists() or committed_at_head(RECORD) is not None:
        raise Refused("a relabelling is already recorded or was started: a rerun is a new ledger row")
    return {"row": ROW, "ruling": ruling, "commit": commit, "clean": True, "config_at_head": True, "n55": n55}


# ------------------------------------------------------------------------------------------------ metadata (whitelist)


def check_whitelist(columns: Iterable[str]) -> tuple[str, ...]:
    """The columns to read must lie inside OBS_WHITELIST and none may look like an annotation, cluster or embedding."""
    cols = tuple(columns)
    bad = [c for c in cols if c not in OBS_WHITELIST or ANNOTATION_PATTERN.search(c)]
    if bad:
        raise Refused(f"columns {bad} may not be read: only {list(OBS_WHITELIST)} (and the cell names) are read")
    return cols


def read_obs(path: str | Path, columns: Sequence[str] = OBS_WHITELIST) -> dict[str, np.ndarray]:
    """The cell names (first column) and the whitelisted `columns` of a metadata table, and nothing else: each row's
    other fields are dropped as the row is parsed, and a requested column outside the whitelist refuses before the file
    is opened. Returns {"cell": names, column: values}."""
    cols = check_whitelist(columns)
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        missing = [c for c in cols if c not in header]
        if missing:
            raise Refused(f"{Path(path).name}: columns {missing} are not in its header")
        idx = [0, *(header.index(c) for c in cols)]
        if len(set(idx)) != len(idx):
            raise Refused("the cell-name column is also a requested column")
        rows = []
        for i, row in enumerate(reader):
            if len(row) != len(header):
                raise Refused(f"{Path(path).name}: data row {i} has {len(row)} fields, header {len(header)}")
            rows.append([row[j] for j in idx])
    arr = np.array(rows, dtype=str).reshape(len(rows), len(idx))
    return {"cell": arr[:, 0], **{c: arr[:, k + 1] for k, c in enumerate(cols)}}


# ------------------------------------------------------------------------------------------------ scores


def log1p_cp10k(C: sparse.spmatrix) -> sparse.csr_matrix:
    """log1p of counts per 10,000 of each cell's own library (no other cell enters)."""
    C = sparse.csr_matrix(C, dtype=np.float64)
    lib = np.asarray(C.sum(1)).ravel()
    if (lib <= 0).any():
        raise Refused("a cell has no counts")
    L = sparse.diags(CP10K / lib) @ C
    L = sparse.csr_matrix(L)
    L.data = np.log1p(L.data)
    return L


def unique_index(genes: Sequence[str]) -> dict[str, int]:
    """Symbol -> column, for symbols Feng's list holds once (a duplicated symbol is not used, as dirlib.sym_index)."""
    seen: dict[str, int] = {}
    dup: set[str] = set()
    for i, g in enumerate(genes):
        if g in seen:
            dup.add(g)
        seen[g] = i
    return {g: i for g, i in seen.items() if g not in dup}


def gene_means(L: sparse.csr_matrix, stage: np.ndarray, stages: Sequence[str]) -> np.ndarray:
    """Each gene's mean log1p CP10k per stage, averaged over `stages` with equal weight."""
    per = []
    for s in stages:
        m = stage == s
        if not m.any():
            raise Refused(f"no cell at stage {s}")
        per.append(np.asarray(L[np.flatnonzero(m)].mean(0)).ravel())
    return np.mean(per, axis=0)


def control_sets(means: np.ndarray, genes: Sequence[str], set_genes: Mapping[str, Sequence[int]],
                 exclude: Iterable[int], n_bins: int, per_gene: int, seed: int) -> dict[str, np.ndarray]:
    """Per set, its fixed control genes: genes with a non-zero mean are ranked by (mean, symbol) and cut into n_bins
    equal-count bins; for each of the set's genes in its bin, per_gene genes are drawn without replacement from that
    bin, leaving out `exclude` (every set gene, every duplicated symbol) and genes already drawn for the set; the
    generator is default_rng([seed, k]) for the k-th set in `set_genes`' order. A set gene with a zero mean has no bin
    and draws nothing."""
    genes = np.asarray(genes, dtype=str)
    pool = np.flatnonzero(means > 0)
    order = pool[np.lexsort((genes[pool], means[pool]))]
    bins = np.array_split(order, n_bins)
    bin_of = {int(g): b for b, members in enumerate(bins) for g in members}
    excl = {int(i) for i in exclude}
    out = {}
    for k, (name, idx) in enumerate(set_genes.items()):
        rng = np.random.default_rng([int(seed), k])
        chosen: list[int] = []
        for g in idx:
            if int(g) not in bin_of:
                continue
            cand = [int(c) for c in bins[bin_of[int(g)]] if int(c) not in excl and int(c) not in chosen]
            take = min(per_gene, len(cand))
            if take:
                chosen += [cand[j] for j in rng.choice(len(cand), size=take, replace=False)]
        out[name] = np.array(sorted(chosen), dtype=np.int64)
    return out


def set_score(L: sparse.csr_matrix, idx: Sequence[int], ctrl: Sequence[int]) -> np.ndarray:
    """Per cell: mean log1p CP10k over the set's genes minus the mean over its control genes (0 when there are none)."""
    a = np.asarray(L[:, list(idx)].mean(1)).ravel()
    b = np.asarray(L[:, list(ctrl)].mean(1)).ravel() if len(ctrl) else np.zeros(L.shape[0])
    return a - b


def assign(scores: Mapping[str, np.ndarray], cfg: Mapping, min_score: float | None = None,
           min_margin: float | None = None) -> np.ndarray:
    """The recomputed label of every cell (assign rule of the yaml; min_score and min_margin override it only for the
    sensitivity rows, which decide nothing)."""
    a = cfg["assign"]
    t = float(a["min_score"] if min_score is None else min_score)
    m = float(a["min_margin"] if min_margin is None else min_margin)
    names = classes(cfg)
    S = np.vstack([np.asarray(scores[n], dtype=np.float64) for n in names])
    if S.shape[0] < 2:
        raise Refused("the assignment needs at least two classes")
    order = np.argsort(-S, axis=0, kind="stable")
    top, second = order[0], order[1]
    cols = np.arange(S.shape[1])
    best, runner = S[top, cols], S[second, cols]
    ok = (best >= t) & (best - runner >= m)
    sp = a["ec_split"]
    lab_of = [cfg["sets"][n]["label"] for n in names]
    out = np.array([lab_of[i] for i in top], dtype=object)
    out[~ok] = a["unassigned"]
    ec = ok & (np.array(names, dtype=object)[top] == sp["class"])
    d = np.asarray(scores[sp["plus"]], dtype=np.float64) - np.asarray(scores[sp["minus"]], dtype=np.float64)
    split_m = float(a["split_margin"])
    out[ec & (d >= split_m)] = sp["plus_label"]
    out[ec & (d <= -split_m)] = sp["minus_label"]
    out[ec & (d > -split_m) & (d < split_m)] = sp["ambiguous_label"]
    if np.any(out == cfg["sets"][sp["class"]]["label"]):  # the split covers every ec cell; a leftover is a bug
        raise Refused("an ec cell was left unsplit")
    return out.astype(str)


def unassigned_shares(labels: np.ndarray, stage: np.ndarray, stages: Sequence[str], unassigned: str) -> dict:
    out = {}
    for s in stages:
        m = stage == s
        out[s] = float(np.mean(labels[m] == unassigned)) if m.any() else None
    return out


def cap_check(shares: Mapping[str, float | None], cap: float) -> dict:
    """The per-stage cap: passed only if every stage has cells and its unassigned share is at most the cap."""
    over = sorted(s for s, v in shares.items() if v is None or v > cap)
    return {"cap": float(cap), "shares": dict(shares), "stages_over": over, "passed": not over}


def group_cells(labels: np.ndarray, stage: np.ndarray, stages: Sequence[str], groups: Mapping[str, Sequence[str]]
                ) -> dict[str, dict[str, int]]:
    return {g: {s: int(np.sum((stage == s) & np.isin(labels, list(labs)))) for s in stages}
            for g, labs in groups.items()}


def fallbacks(cells: Mapping[str, Mapping[str, int]], min_cells: int, fallback: Mapping[str, str]) -> dict:
    """Per group: whether it has its own direction (at least min_cells at both ends), and if not, what it takes: the
    fallback group's direction when that group has its own, else none. One level only: a fallback is never chained."""
    own = {g: min(c.values()) >= min_cells for g, c in cells.items()}
    out = {}
    for g in cells:
        if own[g]:
            out[g] = {"own": True, "fired": False, "takes": g}
        else:
            fb = fallback.get(g, "none")
            out[g] = {"own": False, "fired": True, "takes": fb if fb != "none" and own.get(fb, False) else None}
    return out


def relabel(C: sparse.spmatrix, genes: Sequence[str], stage: np.ndarray, cfg: Mapping) -> dict:
    """The labelling itself, on counts and stage labels only: scores, labels and the rule checks. C holds the cells to
    label (rows), stage their stage labels (each one of the yaml's stages)."""
    stage = np.asarray(stage, dtype=str)
    stages = list(cfg["stages"])
    if not np.isin(stage, stages).all():
        raise Refused("a cell outside the yaml's stages was handed to the labeller")
    L = log1p_cp10k(C)
    uidx = unique_index(list(genes))
    sc = cfg["score"]
    found, missing, set_idx = {}, {}, {}
    for name, s in cfg["sets"].items():
        f = [g for g in s["genes"] if g in uidx]
        found[name], missing[name] = f, [g for g in s["genes"] if g not in uidx]
        if len(f) < int(sc["min_genes_found"]):
            raise Refused(f"set {name}: {len(f)} of its genes are in Feng's symbols (unique); the yaml needs "
                          f"{sc['min_genes_found']}")
        set_idx[name] = [uidx[g] for g in f]
    dup = set(range(len(genes))) - set(uidx.values())
    kept_out = named_genes(cfg) | panel_genes()  # the yaml's score rule: no named or panel gene is ever a control
    exclude = {i for idx in set_idx.values() for i in idx} | dup | {uidx[g] for g in kept_out if g in uidx}
    means = gene_means(L, stage, stages)
    ctrl = control_sets(means, genes, set_idx, exclude, int(sc["control_bins"]), int(sc["control_per_gene"]),
                        int(sc["control_seed"]))
    empty = [n for n, c in ctrl.items() if len(c) == 0]
    if empty:
        raise Refused(f"sets {empty} drew no control gene: their scores would not be control-adjusted")
    scores = {name: set_score(L, set_idx[name], ctrl[name]) for name in cfg["sets"]}
    labels = assign(scores, cfg)
    low = np.vstack([scores[n] for n in classes(cfg)]).max(0) < float(cfg["assign"]["min_score"])
    shares = unassigned_shares(labels, stage, stages, cfg["assign"]["unassigned"])
    cap = cap_check(shares, float(cfg["unassigned_cap"]))
    cells = group_cells(labels, stage, stages, cfg["groups"])
    fb = fallbacks(cells, int(cfg["min_cells_ext"]), cfg["fallback"])
    unassigned_kinds = {s: {"no_class_reached_min_score": int(np.sum((stage == s) & low)),
                            "margin_failed": int(np.sum((stage == s) & (labels == cfg["assign"]["unassigned"]) & ~low))}
                        for s in stages}  # reported, deciding nothing: the cap counts both
    gnames = np.asarray(genes, dtype=str)
    return {"labels": labels, "scores": scores, "unassigned_share": shares, "cap": cap, "group_cells": cells,
            "fallbacks": fb, "route": route(bool(cap["passed"]), fb, cfg), "unassigned_kinds": unassigned_kinds,
            "sets": {n: {"found": found[n], "missing": missing[n], "controls": [str(gnames[i]) for i in ctrl[n]]}
                     for n in cfg["sets"]},
            "label_counts": {s: {lab: int(np.sum((stage == s) & (labels == lab))) for lab in all_labels(cfg)}
                             for s in stages}}


# ------------------------------------------------------------------------------------------------ run


def _write_cells(path: Path, names, stage, zone, labels, scores: Mapping[str, np.ndarray]) -> None:
    tmp = path.with_name(f".{path.name}.partial")
    with open(tmp, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["cell", "stage", "zone", "label", *(f"score:{n}" for n in scores)])
        for i in range(len(names)):
            w.writerow([names[i], stage[i], zone[i], labels[i], *(f"{scores[n][i]:.17g}" for n in scores)])
    tmp.replace(path)


def _rel(p: Path) -> str:
    try:
        return str(Path(p).resolve().relative_to(fio.REPO))
    except ValueError:
        return str(p)


def run(cfg: Mapping, cfg_path: Path, src: Path, dst: Path, gate: Mapping | None) -> dict:
    """Read the clean reduction's permitted cells (manifest-checked), label them, write the statistics, the cells
    table and the manifest. `gate` is what `preflight` returned (main passes it); the manifest records it, and the
    readers refuse a file whose gate is missing or not clean. Returns the manifest."""
    if not gate:
        raise Refused("no gate: the relabelling runs only through main(), after preflight")
    # the gate's three conditions again, for a caller that hands run() a gate of its own (defence in depth)
    if cfg.get("status") != "frozen":
        raise Refused(f"the yaml's status is {cfg.get('status')!r}: it is run only when 'frozen'")
    if gate.get("row") != ROW or gate.get("clean") is not True:
        raise Refused(f"the gate does not show row {ROW} and a clean tree: {dict(gate)}")
    unsourced = unsourced_genes(cfg)
    if unsourced:
        raise Refused(f"the yaml holds UNSOURCED genes {unsourced}: source each one to a named reference that is not "
                      "the Feng paper, or delete it, before any run")
    check_min_cells(cfg, yaml.safe_load(DIRECTION_CONFIG.read_text()))
    inp = cfg["inputs"]
    paths = {k: src / inp[k] for k in ("counts", "genes", "obs")}
    entries = {k: fio.verify_manifest_output(p) for k, p in paths.items()}
    obs = read_obs(paths["obs"], cfg["obs_columns"])
    if not np.all(obs["strain"] == cfg["expect_strain"]):
        raise Refused(f"a cell's strain is not {cfg['expect_strain']}")
    genes = paths["genes"].read_text().splitlines()
    C = sparse.load_npz(paths["counts"]).tocsr()
    if C.shape != (len(obs["cell"]), len(genes)):
        raise Refused(f"counts {C.shape} do not match the obs table ({len(obs['cell'])}) and genes ({len(genes)})")
    keep = np.flatnonzero(np.isin(obs["stage"], list(cfg["stages"])))
    names, stage, zone = obs["cell"][keep], obs["stage"][keep], obs["zone"][keep]
    C = C[keep]
    res = relabel(C, genes, stage, cfg)
    st = feng2022_geo.group_stats(C, stage, zone, res["labels"], names, genes)
    dst.mkdir(parents=True, exist_ok=True)
    out = cfg["output"]
    fio.save_npz_atomic(dst / out["stats"], **st)
    _write_cells(dst / out["cells"], names, stage, zone, res["labels"], res["scores"])
    usable = bool(res["cap"]["passed"])
    outputs = {out["stats"]: {"sha256": fio.sha256_file(dst / out["stats"]), "usable": usable,
                              "groups": int(len(st["keys"])), "half_seed": feng2022_geo.HALF_SEED},
               out["cells"]: {"sha256": fio.sha256_file(dst / out["cells"])}}
    manifest = {"source": SOURCE, "method": cfg["method"], "register_row": "X002",
                "config": {"file": _rel(cfg_path), "sha256": file_sha256(cfg_path), "status": cfg.get("status")},
                "inputs": {k: {"file": _rel(p), "sha256": entries[k]["sha256"]}
                           for k, p in paths.items()},
                "obs_columns_read": ["<cell names>", *cfg["obs_columns"]], "authors_labels_read": False,
                "stages": list(cfg["stages"]), "unassigned_share": res["unassigned_share"], "cap": res["cap"],
                "usable": usable, "group_cells": res["group_cells"], "fallbacks": res["fallbacks"],
                "route": res["route"], "unassigned_kinds": res["unassigned_kinds"],
                "label_counts": res["label_counts"], "sets": res["sets"], "unsourced": unsourced,
                "gate": dict(gate), "outputs": outputs, "provenance": fio.provenance(CODE)}
    fio.write_manifest(dst, manifest)
    return manifest


def write_record(manifest_path: Path, record: Path | None = None) -> Path:
    """The committed copy of the manifest (byte for byte), which the main session commits right after the run."""
    record = RECORD if record is None else record
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_bytes(Path(manifest_path).read_bytes())
    return record


def write_stop(dst: Path, cfg: Mapping, cfg_path: Path, gate: Mapping, reason: str) -> dict:
    """The record of a run that stopped after preflight (terminal.on_stop): the manifest in the output folder and its
    copy at RECORD, with the reason, usable false and the route "none". It lists no output, so no reader can take a
    statistics file from this run (filterio.verify_manifest_output refuses an unlisted file)."""
    manifest = {"source": SOURCE, "method": cfg.get("method"), "register_row": "X002",
                "config": {"file": _rel(cfg_path), "sha256": file_sha256(cfg_path), "status": cfg.get("status")},
                "gate": dict(gate), "stopped": True, "reason": str(reason)[:500], "usable": False,
                "route": cfg["terminal"]["on_stop"], "outputs": {}, "provenance": fio.provenance(CODE)}
    Path(dst).mkdir(parents=True, exist_ok=True)
    fio.write_manifest(dst, manifest)
    write_record(Path(dst) / "MANIFEST.json")
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ruling", required=True, help="the docs/PLAN.md decision-log entry of the captain's ruling")
    a = ap.parse_args(argv)
    cfg_path = CONFIG
    try:  # nothing of Feng's is opened in this block: a refusal here records nothing and changes no route
        cfg = load_config(cfg_path)
        src = (fio.REPO / cfg["inputs"]["dir"]).resolve()
        dst = output_dir_ok(fio.REPO / cfg["output"]["dir"])
        gate = preflight(cfg_path, a.ruling, dst)
    except Refused as e:
        raise SystemExit(f"refused: {e}") from None
    dst.mkdir(parents=True)
    (dst / STARTED).write_text(json.dumps({"gate": gate, "config_sha256": file_sha256(cfg_path)}, indent=1) + "\n")
    try:
        man = run(cfg, cfg_path, src, dst, gate)
    except BaseException as e:  # terminal.on_stop: any stop after preflight takes route none, recorded, never rerun
        stop = write_stop(dst, cfg, cfg_path, gate, f"{type(e).__name__}: {e}")
        print(f"stopped after preflight ({stop['reason']}): route {stop['route']}, by terminal.on_stop; build "
              f"{' '.join(build_command(stop['route'], cfg))}\ncommit {RECORD.relative_to(fio.REPO)} now; nothing is "
              f"rerun under row {ROW}")
        if not isinstance(e, Exception):
            raise
        return STOP_EXIT
    rec = write_record(dst / "MANIFEST.json")
    fb = ", ".join(f"{g}: {v['takes']}" for g, v in man["fallbacks"].items())
    print(f"written {dst}; usable {man['usable']}; route {man['route']}; unassigned share per stage "
          f"{man['unassigned_share']}; fallbacks {{{fb}}}; build {' '.join(build_command(man['route'], cfg))}\n"
          f"commit {rec.relative_to(fio.REPO)} now: N54's listing and the builder refuse until it is committed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
