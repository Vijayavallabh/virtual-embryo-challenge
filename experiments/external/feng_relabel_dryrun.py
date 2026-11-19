"""Row N55 of experiments/harness/TRIES.md (proposed; T1; a diagnostic, not a try; target-free; no Feng file): the dry
run of feng_relabel.yaml's scorer and guards on the challenge's released T1 training stage E9.5, before the yaml is
frozen for row N54 (rule 9, J7: a guard is dry-run on committed data other than the data its row reads).

It labels the released E9.5 cells (data/raw/T1/E9.5_RNA.h5ad, a training stage, not a protected window) with exactly
feng_relabel.relabel, the yaml's sets, score, thresholds and split, on that one stage (the controls are ranked on its
own gene means), and reports:
  - per released label (the organisers' own annotation of released data): its cells and the share given each
    recomputed label;
  - the unassigned share, its two kinds (no class reached min_score; the margin failed) and whether the cap of 0.40
    would have failed at this stage;
  - the four groups' cells at E9.5 and whether any would fall under min_cells_ext (reported only: they follow the
    released file's size, not Feng's);
  - the same counts, shares per released label and cap verdict under each of the yaml's sensitivity settings;
  - `rule`: the yaml's n55_rule applied by code (feng_relabel.n55_outcome): kept, changed (to which setting) or
    guard_failed. Its setting is what the yaml's min_score and min_margin must be at the freeze; nothing else moves
    by N55, and feng_relabel.check_n55 refuses the relabelling otherwise.
It is optimistic: the released-data genes were selected on this same file (sources.released_e95). It is also a
different assay (single nuclei, SoupX) from Feng's whole cells. No Feng value enters, and E14.5's late classes (mural,
immune, lymphatic, definitive blood) are not exercised: the released file has no E14.5.

It reads only the released E9.5 file's counts, gene names and celltype column. It refuses unless row N55 is in
TRIES.md, a docs/PLAN.md decision-log entry named by --ruling names N55, the yaml is still draft and holds no
UNSOURCED gene and no PENDING choice, its first n55_rule setting is the one under `assign`, the tree is clean, and its
output exists neither in the tree nor at HEAD (deleting an untracked output cannot allow a silent rerun). Writes
experiments/external/results/feng_relabel_dryrun/n55.json, which the main session commits at once. Run (DGX campaign
checkout; minutes):
  uv run --frozen --no-sync python -I -B experiments/external/feng_relabel_dryrun.py --ruling "<decision-log id>"
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
import time
from collections.abc import Mapping
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import feng_relabel as fr  # noqa: E402
import filterio as fio  # noqa: E402

ROW = "N55"
STAGE = "E9.5"
RELEASED = fio.REPO / "data" / "raw" / "T1" / "E9.5_RNA.h5ad"
LABEL_KEY = "celltype"
OUT = fr.N55_OUT.parent


def one_stage(cfg: Mapping, stage: str = STAGE) -> dict:
    """The yaml with its stages set to `stage` alone (the scorer, sets, thresholds and guards unchanged)."""
    out = copy.deepcopy(dict(cfg))
    out["stages"] = [stage]
    return out


def by_released_label(released: np.ndarray, recomputed: np.ndarray, labels: list[str]) -> dict:
    """Per released label: its cells and the share of them given each recomputed label."""
    out = {}
    for lab in sorted(set(released.tolist())):
        m = released == lab
        out[lab] = {"n": int(m.sum()), "shares": {r: float(np.mean(recomputed[m] == r)) for r in labels
                                                  if np.any(recomputed[m] == r)}}
    return out


def dry_run(C, genes, released: np.ndarray, cfg: Mapping) -> dict:
    """The dry run itself on counts `C` (cells x genes) with the released labels `released`; no file is read."""
    c1 = one_stage(cfg)
    stage = np.full(C.shape[0], STAGE)
    res = fr.relabel(C, genes, stage, c1)
    labels_all = fr.all_labels(c1)
    sens = []
    for s in c1.get("sensitivity", []):
        lab = fr.assign(res["scores"], c1, s["min_score"], s["min_margin"])
        shares = fr.unassigned_shares(lab, stage, c1["stages"], c1["assign"]["unassigned"])
        sens.append({"setting": dict(s), "group_cells": fr.group_cells(lab, stage, c1["stages"], c1["groups"]),
                     "unassigned_share": shares, "cap": fr.cap_check(shares, float(c1["unassigned_cap"])),
                     "by_released_label": by_released_label(released, lab, labels_all)})
    small = {g: v[STAGE] < int(c1["min_cells_ext"]) for g, v in res["group_cells"].items()}
    return {"by_released_label": by_released_label(released, res["labels"], labels_all),
            "unassigned_share": res["unassigned_share"], "unassigned_kinds": res["unassigned_kinds"],
            "cap": res["cap"], "group_cells": res["group_cells"], "under_min_cells": small,
            "sets": res["sets"], "sensitivity": sens, "rule": fr.n55_outcome(released, res["scores"], c1)}


def refusals(cfg: Mapping, ruling: str | None) -> str | None:
    """Why N55 may not run now, or None (all checked before the released file is read)."""
    if not fr.row_in_tries(ROW):
        return f"{ROW} has no row in experiments/harness/TRIES.md: a row is committed before its run"
    if not fr.ruling_logged(ruling, ROW):
        return f"docs/PLAN.md has no decision-log entry {ruling!r} naming {ROW}"
    if cfg.get("status") != "draft":
        return f"the yaml's status is {cfg.get('status')!r}: N55 runs before the freeze, on the draft"
    if fr.unsourced_genes(cfg):
        return f"the yaml holds UNSOURCED genes {fr.unsourced_genes(cfg)}: the dry run checks the yaml to be frozen"
    if fr.pending_choices(cfg):
        return f"the yaml holds choices still {fr.PENDING}: {fr.pending_choices(cfg)}"
    first = cfg["n55_rule"]["order"][0]
    if {k: float(first[k]) for k in first} != {k: float(cfg["assign"][k]) for k in first}:
        return f"n55_rule's first setting {dict(first)} is not the one under assign: N55 reads the yaml as it stands"
    if fr.N55_OUT.exists() or fr.committed_at_head(fr.N55_OUT) is not None:
        return f"{ROW} has already written its output (in the tree or at HEAD); a rerun is a new row"
    return None


def main(argv=None) -> None:
    from vec import provenance
    from vec.core import counts_of
    from vec.prep import load_stage

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ruling", required=True, help="the docs/PLAN.md decision-log entry ordering N55")
    a = ap.parse_args(argv)
    cfg = fr.load_config()
    why = refusals(cfg, a.ruling)
    if why:
        sys.exit(why)
    commit, unclean = provenance.git_state(fio.REPO)
    if unclean is None or unclean:
        sys.exit(f"the tree is not clean ({unclean})")
    ref = load_stage(RELEASED, lineage=False)
    C = counts_of(ref)
    report = {"label": f"{ROW}: feng_relabel.yaml's scorer and guards on the released T1 E9.5 (training stage; no Feng "
                       "file; decides nothing by itself; optimistic, the markers were selected on this file)",
              "ruling": a.ruling, **dry_run(C, [str(g) for g in ref.var_names],
                                            ref.obs[LABEL_KEY].astype(str).to_numpy(), cfg),
              "inputs": {"released": str(RELEASED.relative_to(fio.REPO)), "sha256": fio.sha256_file(RELEASED),
                         "config_sha256": fio.sha256_file(fr.CONFIG)},
              "stamp": {"commit": commit, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    OUT.mkdir(parents=True, exist_ok=True)
    fr.N55_OUT.write_text(json.dumps(report, indent=1, default=float) + "\n")
    print(json.dumps({k: report[k] for k in ("unassigned_share", "cap", "group_cells", "under_min_cells")}, indent=1))
    print(f"rule: {report['rule']['outcome']}; the yaml's min_score and min_margin at the freeze: "
          f"{report['rule']['setting']}. Commit {fr.N55_OUT.relative_to(fio.REPO)} now.")


if __name__ == "__main__":
    main()
