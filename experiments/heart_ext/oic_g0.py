#!/usr/bin/env python
"""Row H38's G0 for heart file OIC125 (row R70; HEART_MOON1_design.md P1 with HEART_MOON1_challenge.md section 5's
fixes): target-free, not a try, CPU only. Its rules are oic_core.py's constants; outside data only through rows X021
(the Xenium pup's decoded transcripts, 12 um grid) and X022 (the 23 MOSTA bin1 GEM files, bin24 = 12 um); X002's late
rows are refused, so OIC125 has two studies and G0(c) one pair.

The references: E10.5 at the test (the heart builder's own context for OTK125: D's organ cut of E10.5, start E9.5),
and for the counted read E9.5 (task T2h/<=E9.5->E10.5) and E8.75 (task T2h/<=E8.75->E10.5), each D's pool on that
stage (h29b.organ_pools) with its start stage. A practice reference is loaded from its input stages alone: no target
file is opened.

Steps (each writes into --out; nothing is ever deleted or overwritten):
  ranks    per reference, every panel gene ranked per category by one-vs-rest AUROC on the reference's labelled cells
           (oic_core.marker_ranks); ranks.json and its sha256, before any outside file is read
  units    every file of the fetch's MANIFEST.json (checked against fetch_oic.py's pins; each sha256 recomputed as it
           is read) binned on the 12 um grid (oic_units); units/<section>.npz and units.json
  regions  per reference and section: unit typing with the study's markers, the heart region (largest cm component,
           dilated 25 um), a plot of total counts with the region's outline and no label (plots/<ref>/<section>.png),
           regions.json and reader_sheet.json (one line per plot, verdict empty: "heart" or "doubt" to fill)
  g0       needs the reader's verdicts committed at experiments/heart_ext/oic_reader.json, each bound to the region's
           sha256 (a region that is not the one the reader saw, or a missing verdict, counts as a doubt: the section
           drops). Per reference: G0(a) typing agreement and counts, (b) shared genes, (c) the residual test, (d) not
           idle; the chords (chords_<ref>.json); g0.json, g0.txt and the statement line. Exit 0 PASS, 3 STOP (the test
           reference fails), 5 VOID (the test passes, a practice reference fails: no read can be made).
The statement, which the main session copies into row H38's outcome cell before the read:
  G0 H38: PASS. chords sha256: E10.5 <sha256>; E9.5 <sha256>; E8.75 <sha256>.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fetch_oic  # noqa: E402
import oic_core as oc  # noqa: E402
import oic_units as ou  # noqa: E402

ROW, RULING_ROW, FILE = "H38", "R70", "OIC125"
REFS = ("E10.5", "E9.5", "E8.75")
TEST_REF = "E10.5"
PRACTICE_TASK = {"E9.5": "T2h/<=E9.5->E10.5", "E8.75": "T2h/<=E8.75->E10.5"}
REF_DAY = {"E10.5": 10.5, "E9.5": 9.5, "E8.75": 8.75}
TARGET_DAY = {"E10.5": 12.5, "E9.5": 10.5, "E8.75": 10.5}
STAGE_DAY = {"E15.5": 15.5, "E16.5": 16.5, "P1": oc.P1_DAY}
STUDIES = ("X021", "X022")
READER = HERE / "experiments" / "heart_ext" / "oic_reader.json"
RESULTS = HERE / "experiments" / "heart_ext" / "results"
PASS, STOP, VOID = "PASS", "STOP", "VOID"
EXIT = {PASS: 0, STOP: 3, VOID: 5}
STATEMENT = "G0 {row}: PASS. chords sha256: E10.5 {E10_5}; E9.5 {E9_5}; E8.75 {E8_75}."


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def sha256_file(path: Path) -> str:
    return ou.file_sha256(Path(path))


def code_sha256() -> dict:
    base = Path(__file__).resolve().parent
    return {n: sha256_file(base / n) for n in ("fetch_oic.py", "oic_units.py", "oic_core.py", "oic_g0.py")}


def provenance() -> dict:
    import subprocess

    def git(*a):
        return subprocess.run(["git", "-C", str(HERE), *a], capture_output=True, text=True).stdout.strip()
    return {"commit": git("rev-parse", "HEAD"),
            "tree_clean": git("status", "--porcelain", "--untracked-files=no") == ""}


def write_new(path: Path, text: str) -> None:
    """Writes a new file; an existing one refuses (nothing is overwritten)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "x", encoding="utf-8") as fh:
        fh.write(text)


# ---- the references -----------------------------------------------------------------------------------------------
class Ref:
    """One reference: D's pool (counts over the panel, chamber categories), the start stage's categories, the panel."""

    def __init__(self, name: str, pool, start, genes: list[str], record: dict):
        from vec import t2_build
        from vec.core import counts_of

        spec = t2_build.load_t2_chamber_spec()
        self.name, self.genes, self.record = name, [str(g) for g in genes], record
        self.C = counts_of(pool, self.genes)
        self.cats = t2_build.t2_chamber_categories(np.asarray(pool.obs[t2_build.T2_LABEL_KEY]).astype(str), spec)
        merged = np.asarray(start.obs["merged"]).astype(str)
        scats = t2_build.t2_chamber_categories(np.asarray(start.obs[t2_build.T2_LABEL_KEY]).astype(str), spec)
        # the start stage's own cells of each category's B0 type (momentum_shift's min_cells applies to these)
        own = {"cm": merged == "CM", "endocardium": (merged == "Endothelial") & (scats == "endocardium"),
               "vascular_endothelium": (merged == "Endothelial") & (scats == "vascular_endothelium"),
               "epicardium": (scats == "epicardium") & (merged == "Untracked")}
        self.start_cells = {k: int(v.sum()) for k, v in own.items()}
        self.ref_cells = {k: int((self.cats == k).sum()) for k in oc.CATS}


def test_reference(data_dir: Path | None) -> Ref:
    """E10.5 as the heart builder sees it for OTK125 (row H35's builder_context: the builder's own refusals, stages,
    organ cut and view, stopped at its form step)."""
    h35 = _load("h35_for_oic", HERE / "experiments" / "harness" / "h35_heart_trend_forms.py")
    ctx, _ = h35.builder_context(data_dir)
    return Ref(TEST_REF, ctx.inputs["ref"], ctx.inputs["start"], list(ctx.board.genes),
               {"source": "the heart builder's context for OTK125 (D's organ cut of E10.5; start E9.5)"})


def practice_reference(ref: str, data_dir: Path) -> Ref:
    """A practice reference from its input stages only (no target file): D's pool by h29b.organ_pools, as row H34."""
    from vec import t2_build
    from vec.boards import load_boards
    from vec.prep import load_stage

    h35 = _load("h35_for_oic", HERE / "experiments" / "harness" / "h35_heart_trend_forms.py")
    t = h35.payoff.heart_task(PRACTICE_TASK[ref])
    genes = load_boards()[t.panel_board].genes
    whole = load_stage(Path(data_dir) / t.inputs["ref"], genes=genes, drop_duplicates=True)
    start = load_stage(Path(data_dir) / t.inputs["start"], genes=genes, drop_duplicates=True)
    PC = np.asarray(whole.obsm[h35.h34.COORD_KEY], dtype=float)[:, :3]
    pm, _ = h35.h29b.organ_pools(PC, whole.obs[t2_build.T2_LABEL_KEY].to_numpy(),
                                 np.asarray(whole.obs["core"], dtype=bool), t2_build.load_t2_chamber_spec())
    return Ref(ref, whole[pm[h35.h34.POOL]].copy(), start, list(genes),
               {"source": f"D's pool on {t.inputs['ref']} (h29b.organ_pools), start {t.inputs['start']}",
                "task": PRACTICE_TASK[ref]})


def load_refs(data_dir: Path | None) -> dict[str, Ref]:
    out = {TEST_REF: test_reference(data_dir)}
    for r in REFS[1:]:
        out[r] = practice_reference(r, data_dir if data_dir is not None else Path(HERE / "data" / "raw"))
    return out


# ---- step ranks ---------------------------------------------------------------------------------------------------
def ranks_doc(refs: dict[str, Ref]) -> dict:
    return {"row": ROW, "rule": "oic_core.marker_ranks", "refs": {r: oc.marker_ranks(x.C, x.cats, x.genes)
                                                                 for r, x in refs.items()}}


def run_ranks(out: Path, data_dir: Path | None) -> int:
    doc = ranks_doc(load_refs(data_dir))
    doc["provenance"] = provenance()
    text = json.dumps(doc, indent=1, sort_keys=True)
    write_new(Path(out) / "ranks.json", text)
    print(f"ranks.json sha256 {hashlib.sha256(text.encode()).hexdigest()}", flush=True)
    return 0


# ---- step units ---------------------------------------------------------------------------------------------------
def _one_units(args) -> dict:
    rec, data_dir, units_dir, panel = args
    path = Path(data_dir) / rec["name"]
    if rec["study"] == "X022":
        u = ou.mosta_units(path, panel, sha256=rec["sha256"])
    else:
        u = ou.xenium_file_units(path, panel, sha256=rec["sha256"])
    npz = Path(units_dir) / f"{rec['section']}.npz"
    if npz.exists():
        raise ou.Refused(f"{npz} exists (nothing is overwritten)")
    sha = ou.save_units(u, npz)
    return {"section": rec["section"], "study": rec["study"], "stage": rec["stage"], "file": rec["name"],
            "source_sha256": u["sha256"], "npz": npz.name, "npz_sha256": sha, "n_units": int(len(u["bx"])),
            "n_held": len(u["held"]), "n_rows": int(u["n_rows"]), "unit_um": float(u["unit_um"])}


def panel_genes() -> list[str]:
    from vec.boards import load_boards

    return [str(g) for g in load_boards()["T2:heart:test_extrap"].genes]


def run_units(out: Path, data_dir: Path, jobs: int = 8) -> int:
    man = ou.manifest_files(data_dir)
    panel = panel_genes()
    units_dir = Path(out) / "units"
    units_dir.mkdir(parents=True, exist_ok=True)
    work = [(r, data_dir, units_dir, panel) for r in sorted(man.values(), key=lambda r: -int(r["size"]))]
    with ProcessPoolExecutor(max_workers=max(1, jobs)) as pool:
        recs = list(pool.map(_one_units, work))
    doc = {"row": ROW, "panel": panel, "constants": {"mosta_bin_dnb": fetch_oic.MOSTA_BIN_DNB,
                                                     "xenium_grid_um": fetch_oic.XENIUM_GRID_UM},
           "code_sha256": code_sha256(), "sections": sorted(recs, key=lambda r: r["section"]),
           "missing": sorted(e["name"] for e in fetch_oic.PINNED if e["name"] not in man),
           "provenance": provenance()}
    write_new(Path(out) / "units.json", json.dumps(doc, indent=1, sort_keys=True))
    print(f"units: {len(recs)} sections; missing from the manifest: {doc['missing']}", flush=True)
    return 0


def load_units_doc(units_root: Path) -> dict:
    """units.json, each npz checked against its sha256, and the code that made them this checkout's."""
    doc = json.loads((Path(units_root) / "units.json").read_text())
    now = code_sha256()
    for n in ("fetch_oic.py", "oic_units.py"):
        if doc["code_sha256"].get(n) != now[n]:
            raise ou.Refused(f"the units were made by another {n} than this checkout's")
    for r in doc["sections"]:
        if sha256_file(Path(units_root) / "units" / r["npz"]) != r["npz_sha256"]:
            raise ou.Refused(f"units of {r['section']} are not the file units.json records")
    return doc


# ---- step regions -------------------------------------------------------------------------------------------------
def study_sections(doc: dict) -> dict[str, list[dict]]:
    out = {s: [] for s in STUDIES}
    for r in doc["sections"]:
        out[r["study"]].append(r)
    return out


def study_genes(units: dict[str, dict], secs: list[dict], ref_genes: list[str]) -> list[str]:
    """G_s: the reference panel's genes that the study holds (any count in any of its sections), in panel order."""
    held = set().union(*[set(units[r["section"]]["held"]) for r in secs]) if secs else set()
    return [g for g in ref_genes if g in held]


def _cols(u: dict, genes: list[str]) -> sparse.csr_matrix:
    idx = {g: i for i, g in enumerate(u["genes"])}
    return u["panel"].tocsc()[:, [idx[g] for g in genes]].tocsr()


def region_of(u: dict, genes: list[str], markers: dict) -> tuple[dict, np.ndarray]:
    sc = oc.unit_scores(_cols(u, genes), genes, markers)
    seed = sc["typed"] & (sc["argmax"] == 0) & (np.nan_to_num(sc["margin"], nan=-1.0) >= oc.MARGIN_REGION)
    return sc, oc.heart_region(u["bx"], u["by"], seed, u["unit_um"])


def region_sha256(u: dict, region: np.ndarray) -> str:
    h = hashlib.sha256()
    h.update(np.asarray(u["bx"], np.int64)[region].tobytes())
    h.update(np.asarray(u["by"], np.int64)[region].tobytes())
    return h.hexdigest()


def plot_region(u: dict, region: np.ndarray, path: Path) -> str:
    """Total counts per unit (log scale, grey) with the region's outline; a second panel zoomed on the region. No
    title, axis, legend or label of any kind."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bx, by = np.asarray(u["bx"]) - np.min(u["bx"]), np.asarray(u["by"]) - np.min(u["by"])
    img = np.zeros((int(by.max()) + 1, int(bx.max()) + 1))
    img[by, bx] = np.log1p(np.asarray(u["total"], float))
    reg = np.zeros(img.shape, bool)
    reg[by[region], bx[region]] = True
    fig, axes = plt.subplots(1, 2, figsize=(16, 8), dpi=110)
    for ax in axes:
        ax.imshow(img, cmap="gray_r", origin="lower", interpolation="nearest")
        if reg.any():
            ax.contour(reg.astype(float), levels=[0.5], colors=["#d62728"], linewidths=0.8, origin="lower")
        ax.set_axis_off()
    if reg.any():
        ys, xs = np.nonzero(reg)
        pad = 40
        axes[1].set_xlim(max(0, xs.min() - pad), min(img.shape[1], xs.max() + pad))
        axes[1].set_ylim(max(0, ys.min() - pad), min(img.shape[0], ys.max() + pad))
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ou.Refused(f"{path} exists (nothing is overwritten)")
    fig.savefig(path)
    plt.close(fig)
    return sha256_file(path)


def regions_doc(refs_ranks: dict, units: dict[str, dict], doc: dict, ref_genes: dict[str, list[str]],
                plots: Path | None = None) -> dict:
    out = {}
    secs = study_sections(doc)
    for ref, rk in refs_ranks.items():
        out[ref] = {}
        for s in STUDIES:
            genes = study_genes(units, secs[s], ref_genes[ref])
            markers = oc.markers_for(rk, genes)
            for r in secs[s]:
                u = units[r["section"]]
                sc, region = region_of(u, genes, markers)
                rec = {"study": s, "stage": r["stage"], "region_sha256": region_sha256(u, region),
                       "region_units": int(region.sum()), "seed_units": int(((sc["argmax"] == 0) & sc["typed"]).sum()),
                       "markers_sha256": oc.sha256_json(markers)}
                if plots is not None:
                    rec["plot"] = f"plots/{ref}/{r['section']}.png"
                    rec["plot_sha256"] = plot_region(u, region, Path(plots) / ref / f"{r['section']}.png")
                out[ref][r["section"]] = rec
    return out


def run_regions(out: Path, units_root: Path) -> int:
    out = Path(out)
    ranks = json.loads((out / "ranks.json").read_text())
    doc = load_units_doc(units_root)
    units = {r["section"]: ou.load_units(Path(units_root) / "units" / r["npz"]) for r in doc["sections"]}
    ref_genes = {r: ranks["refs"][r]["genes"] for r in REFS}
    reg = regions_doc({r: ranks["refs"][r] for r in REFS}, units, doc, ref_genes, plots=out / "plots")
    write_new(out / "regions.json", json.dumps({"row": ROW, "regions": reg, "provenance": provenance()}, indent=1,
                                              sort_keys=True))
    sheet = {"reader": None, "read_at_utc": None, "regions_commit": provenance()["commit"],
             "verdicts": [{"ref": r, "section": s, "plot": v["plot"], "plot_sha256": v["plot_sha256"],
                           "region_sha256": v["region_sha256"], "verdict": None}
                          for r in REFS for s, v in sorted(reg[r].items())],
             "instructions": ("For each plot: 'heart' if the outlined region is the heart alone (no lung, diaphragm "
                              "or other skeletal muscle), else 'doubt'. A doubt drops that section at that reference.")}
    write_new(out / "reader_sheet.json", json.dumps(sheet, indent=1))
    print(f"regions: {sum(len(v) for v in reg.values())} plots under {out / 'plots'}; "
          f"sheet {out / 'reader_sheet.json'}",
          flush=True)
    return 0


# ---- step g0 ------------------------------------------------------------------------------------------------------
def reader_verdicts(path: Path = READER) -> dict:
    """The committed reader file: (ref, section) -> (verdict, region_sha256). Refused unless it is tracked and as
    committed."""
    import subprocess

    rel = str(Path(path).resolve().relative_to(HERE))
    tracked = subprocess.run(["git", "-C", str(HERE), "ls-files", "--error-unmatch", rel], capture_output=True)
    dirty = subprocess.run(["git", "-C", str(HERE), "status", "--porcelain", "--", rel], capture_output=True,
                           text=True).stdout.strip()
    if tracked.returncode != 0 or dirty:
        raise ou.Refused(f"{rel} is not committed as it stands")
    doc = json.loads(Path(path).read_text())
    out = {}
    for v in doc.get("verdicts", []):
        out[(v["ref"], v["section"])] = (v.get("verdict"), v.get("region_sha256"))
    return out


def ref_means(ref: Ref, genes: list[str], cats: list[str]) -> dict[str, np.ndarray]:
    idx = {g: i for i, g in enumerate(ref.genes)}
    C = ref.C.tocsc()[:, [idx[g] for g in genes]].tocsr()
    code = np.array([oc.CATS.index(c) if c in oc.CATS else -1 for c in ref.cats])
    m = oc.category_means(C, code)
    return {k: m[k] for k in cats if k in m}


def source_typing(units: dict[str, dict], secs: list[dict], genes: list[str], markers: dict, ref: Ref,
                  confirmed: set[str]) -> dict:
    """Per confirmed section: the region units typed by the marker rule and by the 15-NN vote; pooled per study."""
    idx = {g: i for i, g in enumerate(ref.genes)}
    Cref = ref.C.tocsc()[:, [idx[g] for g in genes]].tocsr()
    t1s, t2s, per_stage = [], [], {}
    for r in secs:
        if r["section"] not in confirmed:
            continue
        u = units[r["section"]]
        sc, region = region_of(u, genes, markers)
        t1 = oc.region_typing(sc, region)
        inside = region & sc["typed"]
        Cu = _cols(u, genes)
        t2 = np.full(len(t1), -1)
        if inside.any():
            t2[inside] = oc.knn_vote(Cu[inside], Cref, ref.cats, genes)
        t1s.append(t1[inside])
        t2s.append(t2[inside])
        st = per_stage.setdefault(r["stage"], {"C": [], "t": []})
        st["C"].append(Cu[t1 >= 0])
        st["t"].append(t1[t1 >= 0])
    t1a = np.concatenate(t1s) if t1s else np.zeros(0, int)
    t2a = np.concatenate(t2s) if t2s else np.zeros(0, int)
    stages = {s: (sparse.vstack(v["C"]).tocsr(), np.concatenate(v["t"])) for s, v in per_stage.items()}
    return {"agreement": oc.agreement(t1a, t2a), "stages": stages, "sections": sorted(confirmed & {r["section"]
                                                                                                    for r in secs})}


def chords_at(ref: Ref, units: dict, doc: dict, rk: dict, verdicts: dict) -> dict:
    """G0 (a) to (d) and the chords at one reference. Returns the record; 'pass' tells whether a chord file is made."""
    secs = study_sections(doc)
    why: list[str] = []
    src: dict = {}
    for s in STUDIES:
        genes = study_genes(units, secs[s], ref.genes)
        markers = oc.markers_for(rk, genes)
        confirmed = set()
        for r in secs[s]:
            u = units[r["section"]]
            _, region = region_of(u, genes, markers)
            v, sha = verdicts.get((ref.name, r["section"]), (None, None))
            if v == "heart" and sha == region_sha256(u, region):
                confirmed.add(r["section"])
        typ = source_typing(units, secs[s], genes, markers, ref, confirmed) if confirmed else None
        ag = typ["agreement"] if typ else {"share": float("nan"), "jaccard": {}, "typed": {}}
        usable = [k for k in oc.CATS if ag["typed"].get(k, 0) >= oc.MIN_TYPED[k]
                  and ag["jaccard"].get(k, float("nan")) >= oc.AGREE_JACCARD]
        ok = bool(typ and len(genes) >= oc.MIN_SHARED_GENES and ag["share"] >= oc.AGREE_ALL
                  and ag["typed"].get("cm", 0) >= oc.MIN_TYPED["cm"])
        src[s] = {"genes": genes, "markers": markers, "confirmed": sorted(confirmed), "typing": typ, "agreement": ag,
                  "usable": usable, "pass": ok, "n_genes": len(genes)}
        if len(genes) < oc.MIN_SHARED_GENES:
            why.append(f"(b) {s}: {len(genes)} shared genes, under {oc.MIN_SHARED_GENES}")
        if not ok:
            why.append(f"(a) {s} drops (agreement, cm count, shared genes or no confirmed section)")
    live = [s for s in STUDIES if src[s]["pass"]]
    rec = {"ref": ref.name, "ref_record": ref.record, "ref_cells": ref.ref_cells, "start_cells": ref.start_cells,
           "sources": {s: {k: v for k, v in src[s].items() if k not in ("typing", "genes")} for s in STUDIES}}
    if len(live) < 2:
        return {**rec, "pass": False, "why": [*why, "(a) fewer than two sources left"]}
    panel = ref.genes
    pidx = {g: i for i, g in enumerate(panel)}
    r_study, p_study, fits, centres = {}, {}, {}, {}
    level = {}
    for s in live:
        genes = src[s]["genes"]
        centre = [k for k in oc.CATS if ref.ref_cells[k] >= oc.MIN_REF_CELLS and k in src[s]["usable"]]
        centres[s] = centre
        mref = ref_means(ref, genes, centre)
        cref = oc.contrasts(mref, centre) if len(centre) >= 2 else {}
        r_stage, p_stage = {k: [] for k in centre}, {k: [] for k in centre}
        for stage, (C, t) in src[s]["typing"]["stages"].items():
            msrc = oc.category_means(C, t)
            if not all(k in msrc for k in centre) or len(centre) < 2:
                continue
            csrc = oc.contrasts({k: msrc[k] for k in centre}, centre)
            a = oc.a_of(TARGET_DAY[ref.name], REF_DAY[ref.name], STAGE_DAY[stage])
            for k in centre:
                guard = (msrc[k] >= oc.GUARD_CP10K) & (mref[k] >= oc.GUARD_CP10K)
                r, fit = oc.residual(csrc[k], cref[k], guard)
                fits[f"{s}/{stage}/{k}"] = {**fit, "a": a}
                full_r, full_p = np.full(len(panel), np.nan), np.full(len(panel), np.nan)
                gi = [pidx[g] for g in genes]
                full_r[gi], full_p[gi] = r, a * r
                r_stage[k].append(full_r)
                p_stage[k].append(full_p)
        r_study[s] = {k: oc.study_mean(v) for k, v in r_stage.items() if v}
        p_study[s] = {k: oc.study_mean(v) for k, v in p_stage.items() if v}
        lv = np.full(len(panel), np.nan)
        for k in centre:
            if k not in level:
                level[k] = lv.copy()
                level[k][[pidx[g] for g in genes]] = np.log(mref[k] + oc.PSEUDO)
    tests, P, kept = {}, {}, {}
    excluded = np.zeros(len(panel), bool)
    for s in live:
        for g in oc.typing_genes(src[s]["markers"]):
            excluded[pidx[g]] = True
    a, b = live[0], live[1]
    for k in oc.CATS:
        if not (k in r_study[a] and k in r_study[b]):
            tests[k] = {"pass": False, "why": "not held by both studies"}
            continue
        tests[k] = oc.spearman_test(r_study[a][k], r_study[b][k], level[k], key=k)
        if not tests[k]["pass"]:
            continue
        if ref.start_cells[k] < oc.MIN_START_CELLS:
            tests[k]["no_p"] = "fewer own cells at the start stage than momentum_shift's min_cells"
            continue
        P[k], kept[k] = oc.consensus({s: p_study[s][k] for s in live}, excluded)
    if not tests["cm"]["pass"]:
        why.append("(c) cm fails the residual test")
    n_big = int(sum((np.abs(np.nan_to_num(v)) >= oc.P_MIN).sum() for v in P.values()))
    if n_big < oc.N_P_MIN:
        why.append(f"(d) {n_big} entries with |P| >= {oc.P_MIN}, under {oc.N_P_MIN}")
    rec.update(centre=centres, lineage_tests=tests, fits=fits, consensus=kept, n_big=n_big,
               clipped=int(sum(v["clipped"] for v in kept.values())),
               P={k: [None if not np.isfinite(x) else float(x) for x in v] for k, v in P.items()}, genes=panel,
               a={s: {st: oc.a_of(TARGET_DAY[ref.name], REF_DAY[ref.name], STAGE_DAY[st])
                      for st in src[s]["typing"]["stages"]} for s in live})
    return {**rec, "pass": not why, "why": why}


def chord_file(rec: dict, ranks_sha: str, units_doc: dict, reader_sha: str) -> dict:
    return {"row": ROW, "file": FILE, "ruling": RULING_ROW, "ref": rec["ref"],
            "task": PRACTICE_TASK.get(rec["ref"], "test"), "genes": rec["genes"], "P": rec["P"],
            "rule": "oic_core: residual chords (FIX 1), consensus (FIX 2), lineage test (FIX 3), centre set (FIX 5)",
            "rules": oc.RULES, "ranks_sha256": ranks_sha, "reader_sha256": reader_sha,
            "units": {r["section"]: r["npz_sha256"] for r in units_doc["sections"]},
            "sources": {r["file"]: r["source_sha256"] for r in units_doc["sections"]},
            "code_sha256": code_sha256(), "provenance": provenance()}


def load_chords(path: Path, sha256: str | None = None, genes: list[str] | None = None) -> dict[str, np.ndarray]:
    """A chord file as category -> P over `genes` (default its own order; NaN where none). With sha256, the file must
    be of that sha256."""
    raw = Path(path).read_bytes()
    if sha256 is not None and hashlib.sha256(raw).hexdigest() != sha256:
        raise ou.Refused(f"{path} is not the file of sha256 {sha256}")
    doc = json.loads(raw)
    own = doc["genes"]
    genes = own if genes is None else [str(g) for g in genes]
    idx = {g: i for i, g in enumerate(own)}
    out = {}
    for k, v in doc["P"].items():
        arr = np.array([np.nan if x is None else float(x) for x in v])
        out[k] = np.array([arr[idx[g]] if g in idx else np.nan for g in genes])
    return out


def run_g0(out: Path, units_root: Path, data_dir: Path | None, reader: Path = READER,
           results: Path = RESULTS) -> int:
    out = Path(out)
    ranks_raw = (out / "ranks.json").read_text() if (out / "ranks.json").exists() else None
    refs = load_refs(data_dir)
    now = ranks_doc(refs)
    if ranks_raw is None or json.loads(ranks_raw)["refs"] != json.loads(json.dumps(now["refs"])):
        print("G0 H38 is refused: ranks.json is missing or not the ranks this checkout computes", file=sys.stderr)
        return 2
    ranks_sha = hashlib.sha256(ranks_raw.encode()).hexdigest()
    doc = load_units_doc(units_root)
    units = {r["section"]: ou.load_units(Path(units_root) / "units" / r["npz"]) for r in doc["sections"]}
    verdicts = reader_verdicts(reader)
    reader_sha = sha256_file(reader)
    recs = {r: chords_at(refs[r], units, doc, now["refs"][r], verdicts) for r in REFS}
    shas = {}
    for r in REFS:
        if recs[r]["pass"]:
            text = json.dumps(chord_file(recs[r], ranks_sha, doc, reader_sha), sort_keys=True)
            path = Path(results) / f"chords_{r}.json"
            write_new(path, text)
            shas[r] = hashlib.sha256(text.encode()).hexdigest()
    if not recs[TEST_REF]["pass"]:
        verdict = STOP
    elif not all(recs[r]["pass"] for r in REFS[1:]):
        verdict = VOID
    else:
        verdict = PASS
    line = (STATEMENT.format(row=ROW, E10_5=shas["E10.5"], E9_5=shas["E9.5"], E8_75=shas["E8.75"])
            if verdict == PASS else
            f"G0 {ROW}: {verdict}" + (" AT PRACTICE" if verdict == VOID else "") + ": "
            + "; ".join(f"{r} {w}" for r in REFS for w in recs[r]["why"]) + ".")
    g0 = {"row": ROW, "axis": "h38-oic:g0", "verdict": verdict, "line": line, "chords_sha256": shas,
          "ranks_sha256": ranks_sha, "reader_sha256": reader_sha, "refs": recs, "rules": oc.RULES,
          "provenance": provenance(), "commit": provenance()["commit"]}
    write_new(Path(results) / "g0.json", json.dumps(g0, indent=1, sort_keys=True, default=float))
    write_new(Path(results) / "g0.txt", line + "\n")
    print(line, flush=True)
    return EXIT[verdict]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=("ranks", "units", "regions", "g0"))
    ap.add_argument("--out", type=Path, required=True, help="this G0's working folder (ranks, units, plots)")
    ap.add_argument("--data", type=Path, help="the competition data folder (default data/raw of this checkout)")
    ap.add_argument("--fetched", type=Path, help="the fetch's folder (MANIFEST.json), for step units")
    ap.add_argument("--units", type=Path, help="the folder that holds units.json (default --out)")
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args(argv)
    data = a.data or (HERE / "data" / "raw")
    try:
        if a.step == "ranks":
            return run_ranks(a.out, data)
        if a.step == "units":
            if a.fetched is None:
                print("units needs --fetched", file=sys.stderr)
                return 2
            return run_units(a.out, a.fetched, a.jobs)
        if a.step == "regions":
            return run_regions(a.out, a.units or a.out)
        return run_g0(a.out, a.units or a.out, data)
    except (ou.Refused, FileExistsError) as e:
        print(f"G0 {ROW} is refused: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
