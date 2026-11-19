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
UNSOURCED gene and no PENDING choice, its first n55_rule setting is the one under `assign`, the tree is clean, its
output exists neither in the tree nor at HEAD (deleting an untracked output cannot allow a silent rerun), and no earlier
dry run left its marker (feng_relabel.N55_STARTED, under the checkout's own data/derived, written before the released
file is opened). Writes experiments/external/results/feng_relabel_dryrun/n55.json, which the main session commits at
once.

A dry run that stops after these refusals (a set short of score.min_genes_found in the released file's symbols, a set
that draws no control gene, any other error or interrupt) writes n55.json as a stop record: "stopped": true, the
reason, and "rule": {"outcome": "stopped", "setting": the first setting of n55_rule}; it exits with
feng_relabel.STOP_EXIT (an interrupt is recorded, then raised again). By n55_rule (iv), min_score and min_margin stay,
N55 is not rerun under its row, nothing in the yaml is edited because of it, and N54's first line says J7 was not met.
A dry run killed outright leaves only the marker: `--record-kill` then writes the same stop record (reason "killed"),
and refuses unless the marker exists and n55.json does not (in the tree or at HEAD). Run (DGX campaign checkout;
minutes; never --record-kill while the dry run still runs):
  uv run --frozen --no-sync python -I -B experiments/external/feng_relabel_dryrun.py --ruling "<decision-log id>"
  uv run --frozen --no-sync python -I -B experiments/external/feng_relabel_dryrun.py --record-kill
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
    if fr.N55_STARTED.exists():
        return (f"{ROW} was started before ({fr.N55_STARTED} exists) and wrote no output: record the stop with "
                "--record-kill; a rerun is a new row")
    return None


def stop_report(cfg: Mapping, ruling: str | None, reason: str, config_sha256: str, stamp: Mapping,
                released_sha256: str | None = None) -> dict:
    """n55.json for a dry run that stopped (n55_rule (iv)): the outcome "stopped" with the rule's first setting, so
    min_score and min_margin stay; feng_relabel.check_n55 accepts it and puts the reason in N54's gate."""
    first = dict(cfg["n55_rule"]["order"][0])
    return {"label": f"{ROW}: the dry run on the released T1 E9.5 STOPPED; by n55_rule (iv) the thresholds stay, {ROW} "
                     "is not rerun, and N54 says J7 was not met",
            "ruling": ruling, "stopped": True, "reason": str(reason)[:500],
            "rule": {"outcome": "stopped", "setting": first},
            "inputs": {"released": str(RELEASED.relative_to(fio.REPO)), "sha256": released_sha256,
                       "config_sha256": config_sha256},
            "stamp": dict(stamp)}


def _write(report: Mapping) -> None:
    fr.N55_OUT.parent.mkdir(parents=True, exist_ok=True)
    fr.N55_OUT.write_text(json.dumps(report, indent=1, default=float) + "\n")


def record_kill() -> int:
    """`--record-kill`: n55.json as a stop record for a dry run killed outright (its marker exists, no n55.json in the
    tree or at HEAD), with the ruling, commit and yaml sha256 the marker holds. Returns STOP_EXIT."""
    if not fr.N55_STARTED.exists():
        sys.exit(f"refused: no started dry run ({fr.N55_STARTED} is absent): there is no kill to record")
    if fr.N55_OUT.exists() or fr.committed_at_head(fr.N55_OUT) is not None:
        sys.exit("refused: n55.json exists (in the tree or at HEAD): there is no kill to record")
    info = json.loads(fr.N55_STARTED.read_text())
    if info.get("config_sha256") != fio.sha256_file(fr.CONFIG):
        sys.exit("refused: the yaml differs from the one the killed dry run started with")
    cfg = fr.load_config()
    _write(stop_report(cfg, info.get("ruling"), "killed: the dry run ended without writing n55.json (SIGKILL or the "
                       "OOM killer); recorded afterwards by --record-kill", info["config_sha256"], info["stamp"]))
    print(f"recorded the killed dry run: outcome stopped; min_score and min_margin stay. Commit "
          f"{fr.N55_OUT.relative_to(fio.REPO)} now.")
    return fr.STOP_EXIT


def main(argv=None) -> int:
    from vec import provenance

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    what = ap.add_mutually_exclusive_group(required=True)
    what.add_argument("--ruling", help="the docs/PLAN.md decision-log entry ordering N55")
    what.add_argument("--record-kill", action="store_true",
                      help="record a dry run killed outright (its marker and no n55.json) as a stop")
    a = ap.parse_args(argv)
    if a.record_kill:
        return record_kill()
    cfg = fr.load_config()
    why = refusals(cfg, a.ruling)
    if why:
        sys.exit(why)
    commit, unclean = provenance.git_state(fio.REPO)
    if unclean is None or unclean:
        sys.exit(f"the tree is not clean ({unclean})")
    try:
        fr.output_dir_ok(fr.N55_STARTED.parent)  # inside the checkout's own data/ (never a linked folder)
    except fr.Refused as e:
        sys.exit(f"refused: {e}")
    config_sha = fio.sha256_file(fr.CONFIG)
    stamp = {"commit": commit, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    fr.N55_STARTED.parent.mkdir(parents=True, exist_ok=True)
    fr.N55_STARTED.write_text(json.dumps({"ruling": a.ruling, "config_sha256": config_sha, "stamp": stamp}) + "\n")
    released_sha = None
    try:  # n55_rule (iv): any stop from here on is recorded, with the thresholds kept, and never rerun under N55
        from vec.core import counts_of
        from vec.prep import load_stage

        ref = load_stage(RELEASED, lineage=False)
        C = counts_of(ref)
        released_sha = fio.sha256_file(RELEASED)
        body = dry_run(C, [str(g) for g in ref.var_names], ref.obs[LABEL_KEY].astype(str).to_numpy(), cfg)
    except BaseException as e:
        report = stop_report(cfg, a.ruling, f"{type(e).__name__}: {e}", config_sha, stamp, released_sha)
        _write(report)
        print(f"the dry run stopped ({report['reason']}): outcome stopped, by n55_rule (iv); min_score and min_margin "
              f"stay; {ROW} is not rerun. Commit {fr.N55_OUT.relative_to(fio.REPO)} now.")
        if not isinstance(e, Exception):
            raise
        return fr.STOP_EXIT
    report = {"label": f"{ROW}: feng_relabel.yaml's scorer and guards on the released T1 E9.5 (training stage; no Feng "
                       "file; decides nothing by itself; optimistic, the markers were selected on this file)",
              "ruling": a.ruling, **body,
              "inputs": {"released": str(RELEASED.relative_to(fio.REPO)), "sha256": released_sha,
                         "config_sha256": config_sha},
              "stamp": stamp}
    _write(report)
    print(json.dumps({k: report[k] for k in ("unassigned_share", "cap", "group_cells", "under_min_cells")}, indent=1))
    print(f"rule: {report['rule']['outcome']}; the yaml's min_score and min_margin at the freeze: "
          f"{report['rule']['setting']}. Commit {fr.N55_OUT.relative_to(fio.REPO)} now.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
