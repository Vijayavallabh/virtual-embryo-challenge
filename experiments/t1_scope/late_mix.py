#!/usr/bin/env python
"""Row N64 (experiments/harness/TRIES.md): step 0 of row R54's late-mix line, target-free and not a try. It writes
`shares.json`, the chamber shares row R54's file T1_s2rcMSB_LM16 draws library 3's cells at (the T1 CLOSED design of
the step-change round as its cold challenge fixes it: the map by chamber-classes rule 3, library 3 whole, an
independent recomputation to 1e-9).

Inputs, all target-free:
- library 3's own chamber mix by the builder's own call: `harness.t1_chamber_shares` on E10.5's library-3 cells (the
  cell names ending in _3), as build_test.build_t1 takes it for --base lib3 (the unlocked E10.5 file is the test
  prediction's input, A20); no other E10.5 cell and no truth of any task is read;
- the per-label counts of Feng's permitted CD1 cells at E15.5, E16.5 and E17.5 in row N57's committed manifest
  (experiments/external/results/feng_relabel_late/MANIFEST.json: the frozen relabeller's recomputed groups, no
  author label, no expression value);
- for the change-transfer reading only (reported, deciding nothing, entering no file): the E9.5 and E14.5 counts in
  row N54's committed manifest (experiments/external/results/feng_relabel/MANIFEST.json).

The map (configs/t1_chamber_classes.yaml rule 3; the cold challenge's fix of the design's rule-2 map):
  heart    r_cm, r_endocardial_ec, r_vascular_ec, r_ec_ambiguous, r_lymphatic, r_epicardial
  blood    r_blood, r_immune (T1's blood category holds Macrophage)
  out      r_mesenchyme and r_mural (names that do not say they are chamber tissue, with split evidence: rule 3's
           default, surrounding), r_neural_crest, r_hepatic, r_epithelial
  dropped  r_unassigned
The design's map (r_mesenchyme and r_mural in heart) is reported beside it, deciding nothing. A label of the manifest
that the map does not place refuses the step.

The rule, every constant fixed here, for s = cardiomyocytes of heart cells and s = blood of chamber cells:
  logit s(12.5) = logit s_lib3 + 2 / (t_A - 10.5) x (logit s_A - logit s_lib3)
s_lib3 is library 3's own share; s_A the anchor's pooled share (the three stages' counts summed); t_A the anchor's
cell-weighted mean day over its chamber cells under the map (heart and blood cells, days 15.5, 16.5 and 17.5). Then
the moved mix is `t1_scope.with_share(own, "cm", s_cm)` and `with_share(that, "blood", s_blood, within=list(that))`,
the order the builder's --cm-share and --blood-share take, the other categories in proportion.

Independent recomputation (a second code path for every number the builder takes): library 3's mix counted from the
chamber config's yaml directly, the anchor's shares and day from arrays, the logit rule in odds form with math, the
moved mix in closed form. Any value more than 1e-9 from the main path stops the step (exit 1) and nothing is written.

Outputs, in experiments/t1_scope/results/n64_late_mix/ (never overwritten; the main session commits them):
  shares.json  the file the builder reads by its sha256 (row N64 names it as `shares.json sha256 <hex>`): the row,
               the anchor, library 3's own mix, the two moved shares and the moved mix, the map, the formula and the
               inputs' sha256; a record of shares only, no count of any board
  n64.json     the record: shares.json's fields, the design's map's shares, the change-transfer reading, the
               recomputation's largest difference, and row T72's bar (c) threshold (the logit midpoint of the own and
               moved cardiomyocyte shares plus 0.03)
  n64.txt      the same in words
Refused before anything loads: no single row N64 whose first cell names this script as a whole path; HEAD not the
commit given; a tree that is not clean (tracked files); an output that exists; the X002 amendment of 2026-11-25 not in
EXTERNAL_DATA.md's row X002; a manifest that is not the frozen relabeller's record (config sha256, stages, gate row,
authors_labels_read false). Prints the shares file's sha256 and the record's path only.

Usage (repo root, through experiments/t1_scope/n64_step0.sh):
  uv run --frozen --no-sync python -I -B experiments/t1_scope/late_mix.py --commit <full hash of row N64's commit>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import yaml

from vec import harness, t1_scope
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TRIES = ROOT / "experiments" / "harness" / "TRIES.md"
REGISTER = ROOT / "EXTERNAL_DATA.md"
OUT = HERE / "results" / "n64_late_mix"
OUT_FILES = ("shares.json", "n64.json", "n64.txt")
ROW = "N64"
SCRIPT_PATH = "experiments/t1_scope/late_mix.py"
LATE_MIX = "e15_17"  # build_test.py --late-mix's value (row R54)
FILE = "T1_s2rcMSB_LM16"
E105 = "T1/E10.5_RNA.h5ad"
SOURCE = "E10.5_RNA"
DISSECTED = "3"  # build_test.T1_DISSECTED_LIBRARY: the library whose cell names end in _3
LATE_MANIFEST = ROOT / "experiments" / "external" / "results" / "feng_relabel_late" / "MANIFEST.json"  # row N57
EARLY_MANIFEST = ROOT / "experiments" / "external" / "results" / "feng_relabel" / "MANIFEST.json"  # row N54
RELABEL_YAML_SHA256 = "883e74c2c1ef5f01d5244555f8b3dcae24255568c731a6d53ee2788b7b677b49"  # frozen (row N54)
ANCHOR_DAYS = {"E15.5": 15.5, "E16.5": 16.5, "E17.5": 17.5}
EARLY_DAYS = {"E9.5": 9.5, "E14.5": 14.5}
REFERENCE_DAY = 10.5  # E10.5, the day of library 3's cells
TEST_DAY = 12.5  # E12.5, T1's test stage
TRANSFER_FRACTION = (TEST_DAY - REFERENCE_DAY) / (EARLY_DAYS["E14.5"] - EARLY_DAYS["E9.5"])  # 2/5, reported only
HEART = ("r_cm", "r_endocardial_ec", "r_vascular_ec", "r_ec_ambiguous", "r_lymphatic", "r_epicardial")
BLOOD = ("r_blood", "r_immune")
OUT_LABELS = ("r_mesenchyme", "r_mural", "r_neural_crest", "r_hepatic", "r_epithelial")
DROPPED = ("r_unassigned",)
CM = "r_cm"
DESIGN_HEART = (*HEART, "r_mural", "r_mesenchyme")  # the design's rule-2 map, reported only
DESIGN_OUT = tuple(x for x in OUT_LABELS if x not in ("r_mural", "r_mesenchyme"))
TOL = 1e-9  # the independent recomputation
BAR_C_MARGIN = 0.03  # row T72's bar (c): the logit midpoint plus this
AMENDMENT_PHRASE = "Amendment of 2026-11-25 (rows N64, T72 and R54)"
MAP_RULE = ("configs/t1_chamber_classes.yaml rule 3: r_mesenchyme and r_mural carry names that do not say they are "
            "chamber tissue and their evidence is split (the relabeller defines r_mesenchyme by aSHF, pSHF, aPHM, "
            "pPHM, ST and Pericardium markers, outside T1's chamber, and through its late genes E14.5 fibroblasts), "
            "so rule 3's default, surrounding, holds; r_blood and r_immune are T1's blood category; r_unassigned is "
            "dropped")
FORMULA = ("logit s(12.5) = logit s_lib3 + 2 / (t_A - 10.5) x (logit s_A - logit s_lib3), for s = cardiomyocytes of "
           "heart cells and s = blood of chamber cells; s_A the anchor's pooled share, t_A its cell-weighted mean day "
           "over its chamber cells; then t1_scope.with_share(own, 'cm', s_cm) and with_share(that, 'blood', s_blood, "
           "within=list(that))")


# ---------------------------------------------------------------------------------------------------------------------
# Refusals


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def names_path(cell: str, path: str) -> bool:
    """`cell` names `path` as a whole path, not inside a longer one."""
    return re.search(r"(^|[^A-Za-z0-9_./-])" + re.escape(path) + r"([^A-Za-z0-9_./-]|$)", cell) is not None


def row_first_cell(row: str, tries: Path = TRIES) -> str | None:
    """The first cell of the one line opening `| <row> ` of the ledger, or None unless there is exactly one."""
    lines = [ln for ln in tries.read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {row} ")]
    return lines[0][2:].split("|", 1)[0] if len(lines) == 1 else None


def amendment_in_register(register: Path = REGISTER) -> bool:
    """EXTERNAL_DATA.md's row X002 holds the amendment of 2026-11-25 (docs/proposed_X002_amendment.md)."""
    rows = [ln for ln in register.read_text(encoding="utf-8").splitlines() if ln.startswith("| X002 |")]
    return len(rows) == 1 and AMENDMENT_PHRASE in rows[0]


def manifest_reasons(m: Mapping, stages: Sequence[str], gate_row: str) -> list[str]:
    """Why a committed relabelling manifest is not the frozen relabeller's record for `stages` (empty: it is)."""
    bad = []
    if (m.get("config") or {}).get("sha256") != RELABEL_YAML_SHA256:
        bad.append("its settings file is not the frozen one")
    if list(m.get("stages") or []) != list(stages):
        bad.append(f"its stages are {m.get('stages')}, not {list(stages)}")
    if (m.get("gate") or {}).get("row") != gate_row:
        bad.append(f"its gate row is {(m.get('gate') or {}).get('row')}, not {gate_row}")
    if m.get("authors_labels_read") is not False:
        bad.append("it does not say authors_labels_read false")
    if m.get("usable") is not True:
        bad.append("it is not usable")
    counts = m.get("label_counts") or {}
    if sorted(counts) != sorted(stages):
        bad.append("its label counts do not cover its stages")
    placed = set(HEART) | set(BLOOD) | set(OUT_LABELS) | set(DROPPED)
    for st, c in counts.items():
        unplaced = sorted(set(c) - placed)
        if unplaced:
            bad.append(f"{st}: labels the map does not place: {unplaced}")
        if any(not isinstance(v, int) or v < 0 for v in c.values()):
            bad.append(f"{st}: a count is not a whole number")
    return bad


def _rel(path: Path) -> str:
    try:
        return Path(path).resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def refuse(commit: str | None, out: Path | None = None, tries: Path | None = None,
           register: Path | None = None) -> dict:
    """Every refusal before anything loads (module docstring). Returns the start state (commit, clean)."""
    out = OUT if out is None else Path(out)
    tries = TRIES if tries is None else tries
    register = REGISTER if register is None else register
    first = row_first_cell(ROW, tries)
    if first is None or not first.startswith(f"{ROW} (") or not names_path(first, SCRIPT_PATH):
        sys.exit(f"{ROW} is refused: experiments/harness/TRIES.md has no single row {ROW} whose first cell starts "
                 f"'{ROW} (' and names {SCRIPT_PATH}")
    head = harness._git("rev-parse", "HEAD")
    if not commit or head != commit:
        sys.exit(f"{ROW} is refused: HEAD is not the commit given (the full hash of row {ROW}'s commit)")
    if harness._git("status", "--porcelain", "--untracked-files=no") != "":
        sys.exit(f"{ROW} is refused: the tracked tree is not clean")
    there = [f for f in OUT_FILES if (out / f).exists()]
    if there:
        sys.exit(f"{ROW} is refused: {there} exist in {out} (a step's output is never overwritten)")
    if not amendment_in_register(register):
        sys.exit(f"{ROW} is refused: EXTERNAL_DATA.md's row X002 does not hold \"{AMENDMENT_PHRASE}\"")
    for path, stages, gate in ((LATE_MANIFEST, tuple(ANCHOR_DAYS), "N57"), (EARLY_MANIFEST, tuple(EARLY_DAYS), "N54")):
        bad = manifest_reasons(json.loads(path.read_text()), stages, gate)
        if bad:
            sys.exit(f"{ROW} is refused: {_rel(path)}: {'; '.join(bad)}")
    return {"commit": head, "tree_clean": True}


# ---------------------------------------------------------------------------------------------------------------------
# The main path


def library3_pool(ref):
    """E10.5's library-3 cells, as build_test.build_t1 takes them for --base lib3 (cell names ending in _3)."""
    names = [str(x) for x in ref.obs_names]
    library = np.array([x.rsplit("_", 1)[-1] if "_" in x else "" for x in names])
    kinds = set(library.tolist())
    if not all(k.isdigit() for k in kinds) or DISSECTED not in kinds or len(kinds) < 2:
        raise SystemExit(f"{ROW}: the reference's cell names do not separate library 3 from the others")
    return ref[np.flatnonzero(library == DISSECTED)].copy()


def own_shares(ref) -> dict[str, float]:
    """Library 3's own chamber mix by the builder's call (`harness.t1_chamber_shares` on the pool, file stem E10.5)."""
    return harness.t1_chamber_shares(harness._shallow(library3_pool(ref)), SOURCE)


def pooled(counts: Mapping[str, Mapping[str, int]], stages: Sequence[str]) -> dict[str, int]:
    """The labels' counts summed over `stages`."""
    out: dict[str, int] = {}
    for st in stages:
        for lab, n in counts[st].items():
            out[lab] = out.get(lab, 0) + int(n)
    return out


def anchor_shares(counts: Mapping[str, Mapping[str, int]], days: Mapping[str, float],
                  heart: Sequence[str] = HEART, blood: Sequence[str] = BLOOD) -> dict:
    """The anchor's pooled shares under a map (cardiomyocytes of heart cells, blood of chamber cells) and its
    cell-weighted mean day over its chamber cells (heart and blood)."""
    tot = pooled(counts, list(days))
    h = sum(tot.get(x, 0) for x in heart)
    b = sum(tot.get(x, 0) for x in blood)
    chamber = {st: sum(counts[st].get(x, 0) for x in (*heart, *blood)) for st in days}
    mean_day = sum(days[st] * chamber[st] for st in days) / sum(chamber.values())
    return {"cm_share_of_heart": tot.get(CM, 0) / h, "blood_share_of_chamber": b / (h + b), "mean_day": mean_day,
            "heart_cells": h, "blood_cells": b, "chamber_cells_by_stage": chamber}


def logit(p: float) -> float:
    return float(np.log(p / (1.0 - p)))


def expit(x: float) -> float:
    return float(1.0 / (1.0 + np.exp(-x)))


def toward(own: float, anchor: float, weight: float) -> float:
    """logit s = logit own + weight x (logit anchor - logit own)."""
    return expit(logit(own) + weight * (logit(anchor) - logit(own)))


def moved_mix(own: Mapping[str, float], cm: float, blood: float) -> dict[str, float]:
    """The builder's order: cardiomyocytes of heart cells, then blood of all the mix's categories (of chamber cells)."""
    m = t1_scope.with_share(dict(own), "cm", cm)
    return t1_scope.with_share(m, "blood", blood, within=list(m))


def compute(own: Mapping[str, float], late_counts: Mapping, early_counts: Mapping) -> dict:
    """Every number of the step on the main path (module docstring), both maps' anchors and the change transfer."""
    heart_cats = list(t1_scope.load_spec().heart)
    own_cm = t1_scope.share_of(own, "cm", heart_cats)
    own_blood = float(own.get("blood", 0.0))
    anchor = anchor_shares(late_counts, ANCHOR_DAYS)
    weight = (TEST_DAY - REFERENCE_DAY) / (anchor["mean_day"] - REFERENCE_DAY)
    cm = toward(own_cm, anchor["cm_share_of_heart"], weight)
    blood = toward(own_blood, anchor["blood_share_of_chamber"], weight)
    design = anchor_shares(late_counts, ANCHOR_DAYS, DESIGN_HEART)
    w_design = (TEST_DAY - REFERENCE_DAY) / (design["mean_day"] - REFERENCE_DAY)
    transfer = {}
    for name, heart in (("rule3_map", HEART), ("design_map", DESIGN_HEART)):
        e9 = anchor_shares({"E9.5": early_counts["E9.5"]}, {"E9.5": 9.5}, heart)
        e14 = anchor_shares({"E14.5": early_counts["E14.5"]}, {"E14.5": 14.5}, heart)
        transfer[name] = {
            "feng_cm_share_of_heart": {"E9.5": e9["cm_share_of_heart"], "E14.5": e14["cm_share_of_heart"]},
            "feng_blood_share_of_chamber": {"E9.5": e9["blood_share_of_chamber"],
                                            "E14.5": e14["blood_share_of_chamber"]},
            "cm_share_of_heart": expit(logit(own_cm) + TRANSFER_FRACTION * (logit(e14["cm_share_of_heart"])
                                                                             - logit(e9["cm_share_of_heart"]))),
            "blood_share_of_chamber": expit(logit(own_blood) + TRANSFER_FRACTION * (
                logit(e14["blood_share_of_chamber"]) - logit(e9["blood_share_of_chamber"])))}
    mid = expit((logit(own_cm) + logit(cm)) / 2.0)
    return {"own_shares": dict(own), "own_cm_share_of_heart": own_cm, "own_blood_share_of_chamber": own_blood,
            "anchor": anchor, "weight": weight, "cm_share_of_heart": cm, "blood_share_of_chamber": blood,
            "shares": moved_mix(own, cm, blood),
            "design_map": {"anchor": design, "weight": w_design,
                           "cm_share_of_heart": toward(own_cm, design["cm_share_of_heart"], w_design),
                           "blood_share_of_chamber": toward(own_blood, design["blood_share_of_chamber"], w_design)},
            "change_transfer": transfer,
            "t72_bar_c": {"own_cm_share_of_heart": own_cm, "moved_cm_share_of_heart": cm, "logit_midpoint": mid,
                          "margin": BAR_C_MARGIN, "threshold": mid + BAR_C_MARGIN}}


# ---------------------------------------------------------------------------------------------------------------------
# The independent recomputation


def own_shares_direct(ref, config: Path = t1_scope.CHAMBER_CLASSES) -> dict[str, float]:
    """Library 3's chamber mix counted from the chamber config's yaml directly: each label's category (its per-file
    override for E10.5_RNA where one is written), the categories whose scope is chamber, cells counted with numpy."""
    cfg = yaml.safe_load(Path(config).read_text())
    labels = dict(cfg["labels"])
    labels.update((cfg.get("file_overrides") or {}).get(SOURCE) or {})
    chamber = [c for c, v in cfg["categories"].items() if v["scope"] == "chamber"]
    names = np.asarray([str(x) for x in ref.obs_names])
    lib3 = np.char.rpartition(names.astype(str), "_")[:, 2] == DISSECTED
    lab, n = np.unique(ref.obs["celltype"].astype(str).to_numpy()[lib3], return_counts=True)
    counts = {c: 0 for c in chamber}
    for x, k in zip(lab, n, strict=True):
        cat = labels[str(x)]
        if cat in counts:
            counts[cat] += int(k)
    total = sum(counts.values())
    return {c: v / total for c, v in counts.items() if v > 0}


def recompute(own_direct: Mapping[str, float], late_counts: Mapping) -> dict:
    """The builder's numbers by a second path: arrays for the anchor, the odds form with math for the rule, the
    closed form for the moved mix (cm: s_cm (1 - s_b); another heart category h: own_h / (H - own_cm) (1 - s_cm)
    (1 - s_b), H the heart categories' joint own share; blood: s_b; every other category 0)."""
    stages = list(ANCHOR_DAYS)
    mat = np.array([[late_counts[st].get(x, 0) for x in (*HEART, *BLOOD)] for st in stages], dtype=np.float64)
    heart_n, blood_n = mat[:, :len(HEART)].sum(), mat[:, len(HEART):].sum()
    cm_n = mat[:, HEART.index(CM)].sum()
    s_a_cm, s_a_blood = cm_n / heart_n, blood_n / (heart_n + blood_n)
    t_a = float(np.average([ANCHOR_DAYS[st] for st in stages], weights=mat.sum(axis=1)))
    w = (TEST_DAY - REFERENCE_DAY) / (t_a - REFERENCE_DAY)
    heart_cats = ("cm", "endocardium", "epicardium", "vascular_endothelium")
    big_h = sum(own_direct.get(c, 0.0) for c in heart_cats)
    own_cm = own_direct.get("cm", 0.0) / big_h
    own_blood = own_direct.get("blood", 0.0)

    def odds_toward(p: float, q: float) -> float:
        o = (p / (1 - p)) * math.pow((q / (1 - q)) / (p / (1 - p)), w)
        return o / (1 + o)

    s_cm, s_b = odds_toward(own_cm, s_a_cm), odds_toward(own_blood, s_a_blood)
    rest = big_h - own_direct.get("cm", 0.0)
    cats = t1_scope.load_spec().categories
    shares = {c: 0.0 for c in cats}
    for c in heart_cats:
        shares[c] = s_cm * (1 - s_b) if c == "cm" else own_direct.get(c, 0.0) / rest * (1 - s_cm) * (1 - s_b)
    shares["blood"] = s_b
    return {"own_shares": dict(own_direct), "anchor_cm": s_a_cm, "anchor_blood": s_a_blood, "mean_day": t_a,
            "weight": w, "cm_share_of_heart": s_cm, "blood_share_of_chamber": s_b, "shares": shares}


def largest_difference(main: Mapping, again: Mapping) -> float:
    """The largest absolute difference between the main path and the recomputation over every number the builder
    takes (library 3's mix, the anchor's two shares and day, the weight, the two moved shares, the moved mix)."""
    pairs = [(main["anchor"]["cm_share_of_heart"], again["anchor_cm"]),
             (main["anchor"]["blood_share_of_chamber"], again["anchor_blood"]),
             (main["anchor"]["mean_day"], again["mean_day"]), (main["weight"], again["weight"]),
             (main["cm_share_of_heart"], again["cm_share_of_heart"]),
             (main["blood_share_of_chamber"], again["blood_share_of_chamber"])]
    if set(main["own_shares"]) != set(again["own_shares"]) or set(main["shares"]) != set(again["shares"]):
        return float("inf")
    pairs += [(main["own_shares"][c], again["own_shares"][c]) for c in main["own_shares"]]
    pairs += [(main["shares"][c], again["shares"][c]) for c in main["shares"]]
    return float(max(abs(float(a) - float(b)) for a, b in pairs))


# ---------------------------------------------------------------------------------------------------------------------
# The outputs


def shares_record(res: Mapping, inputs: Mapping, state: Mapping) -> dict:
    """shares.json: what the builder reads (row, anchor, own mix, the two moved shares, the moved mix) and how they
    were made; shares and hashes only."""
    return {"row": ROW, "late_mix": LATE_MIX, "file": FILE, "own_shares": res["own_shares"],
            "own_cm_share_of_heart": res["own_cm_share_of_heart"],
            "own_blood_share_of_chamber": res["own_blood_share_of_chamber"],
            "cm_share_of_heart": res["cm_share_of_heart"], "blood_share_of_chamber": res["blood_share_of_chamber"],
            "shares": res["shares"],
            "anchor": {"stages": list(ANCHOR_DAYS), "cm_share_of_heart": res["anchor"]["cm_share_of_heart"],
                       "blood_share_of_chamber": res["anchor"]["blood_share_of_chamber"],
                       "mean_day": res["anchor"]["mean_day"]},
            "weight": res["weight"], "reference_day": REFERENCE_DAY, "test_day": TEST_DAY, "formula": FORMULA,
            "map": {"heart": list(HEART), "blood": list(BLOOD), "out": list(OUT_LABELS), "dropped": list(DROPPED),
                    "rule": MAP_RULE},
            "inputs": dict(inputs), "commit": state["commit"], "tree_clean": state["tree_clean"]}


def record_text(rec: Mapping) -> str:
    s, d = rec["shares_json"], rec["design_map"]
    lines = [f"Row {ROW} (step 0 of row R54's late-mix line; target-free, not a try), commit {rec['commit']}",
             f"shares.json sha256 {rec['shares_sha256']}  (row {ROW} names it as `shares.json sha256 <hex>`)",
             f"library 3's own mix: cardiomyocytes {s['own_cm_share_of_heart']:.6f} of heart cells, blood "
             f"{s['own_blood_share_of_chamber']:.6f} of chamber cells",
             f"anchor (Feng permitted CD1, E15.5 to E17.5, rule-3 map): cardiomyocytes "
             f"{s['anchor']['cm_share_of_heart']:.6f}, blood {s['anchor']['blood_share_of_chamber']:.6f}, mean day "
             f"{s['anchor']['mean_day']:.4f}; weight {s['weight']:.6f}",
             f"moved (the file's): cardiomyocytes {s['cm_share_of_heart']:.6f} of heart cells, blood "
             f"{s['blood_share_of_chamber']:.6f} of chamber cells",
             "moved mix: " + ", ".join(f"{c} {v:.6f}" for c, v in s["shares"].items() if v > 0),
             f"reported, deciding nothing: the design's map (r_mesenchyme and r_mural in heart) gives cardiomyocytes "
             f"{d['cm_share_of_heart']:.6f} and blood {d['blood_share_of_chamber']:.6f}"]
    for name, t in rec["change_transfer"].items():
        lines.append(f"reported, deciding nothing: the change transfer (Feng E9.5 -> E14.5 at 2/5, {name}) gives "
                     f"cardiomyocytes {t['cm_share_of_heart']:.6f} and blood {t['blood_share_of_chamber']:.6f}")
    b = rec["t72_bar_c"]
    lines.append(f"row T72's bar (c): the logit midpoint of {b['own_cm_share_of_heart']:.6f} and "
                 f"{b['moved_cm_share_of_heart']:.6f} is {b['logit_midpoint']:.6f}; plus {b['margin']:g}: "
                 f"{b['threshold']:.6f}")
    lines.append(f"independent recomputation: largest difference {rec['recomputation']['largest_difference']:.3g} "
                 f"(tolerance {TOL:g})")
    return "\n".join(lines) + "\n"


def run(commit: str | None, data: Path | None = None, out: Path | None = None) -> dict:
    data = harness.DATA_DIR if data is None else Path(data)
    out = OUT if out is None else Path(out)
    state = refuse(commit, out)
    late = json.loads(LATE_MANIFEST.read_bytes())
    early = json.loads(EARLY_MANIFEST.read_bytes())
    ref = load_stage(Path(data) / E105, lineage=False)
    res = compute(own_shares(ref), late["label_counts"], early["label_counts"])
    again = recompute(own_shares_direct(ref), late["label_counts"])
    worst = largest_difference(res, again)
    if not worst <= TOL:
        sys.exit(f"{ROW}: the independent recomputation differs from the main path by {worst:.3g} (tolerance {TOL:g}); "
                 "nothing is written")
    inputs = {"E10.5": {"file": E105, "sha256": sha256_file(Path(data) / E105)},
              "chamber_config": {"file": "configs/t1_chamber_classes.yaml",
                                 "sha256": sha256_file(t1_scope.CHAMBER_CLASSES)},
              "late_manifest": {"file": _rel(LATE_MANIFEST), "row": "N57",
                                "sha256": sha256_file(LATE_MANIFEST)},
              "early_manifest": {"file": _rel(EARLY_MANIFEST), "row": "N54",
                                 "sha256": sha256_file(EARLY_MANIFEST), "use": "the change transfer, reported only"},
              "script": {"file": SCRIPT_PATH, "sha256": sha256_file(Path(__file__))}}
    shares = shares_record(res, inputs, state)
    raw = (json.dumps(shares, indent=1, sort_keys=True) + "\n").encode()
    rec = {"row": ROW, "commit": state["commit"], "shares_json": shares, "shares_sha256": sha256_bytes(raw),
           "design_map": res["design_map"], "change_transfer": res["change_transfer"], "t72_bar_c": res["t72_bar_c"],
           "recomputation": {"largest_difference": worst, "tolerance": TOL, "values": again}}
    out.mkdir(parents=True, exist_ok=True)
    (out / OUT_FILES[0]).write_bytes(raw)
    (out / OUT_FILES[1]).write_text(json.dumps(rec, indent=1) + "\n")
    (out / OUT_FILES[2]).write_text(record_text(rec))
    print(f"shares.json sha256 {rec['shares_sha256']}\n(the record in {out / OUT_FILES[2]})")
    return rec


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--commit", required=True, help="the full hash of row N64's commit: HEAD must be it")
    ap.add_argument("--data", default=None)
    a = ap.parse_args(argv)
    run(a.commit, a.data)


if __name__ == "__main__":
    main()
