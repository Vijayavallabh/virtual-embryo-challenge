#!/usr/bin/env python
"""T1 tries T13-T18 of experiments/harness/TRIES.md: the late-proxy tournament, run and decided. The captain approved
Feng's E14.5 -> E18.5 change as a stand-in direction on T1's own released backtests on 2026-11-10 (docs/PLAN.md decision
log 2026-11-10 (a); EXTERNAL_DATA.md X002 amendment). No counted task can judge a change to Feng's real direction,
E9.5 -> E14.5, whose E9.5 lies inside every T1 pseudo-task's window; (X002, 14.5, 18.5) is legal on all three T1
pseudo-tasks (vec.harness.check_externals) and it is the first counted check of direction changes on T1's own platform
(nuclei, SoupX).

Candidates (vec.harness.CANDIDATES, one try each): lp_proxy (T13, the stand-in itself), lp_cflip and lp_cperm (T14 and
T15, the sign-flipped and gene-permuted controls), lp_hyg (T16, R2's D and A1 held flat), lp_hyg_n (T17, also N1 and N2)
and lp_eb (T18, reliability shrinkage, one factor per expression decile). The base is b0_scope: B0's chamber copy with
no momentum, the same cells at the same seed, so a difference is the direction alone, and M1 keeps its cells without a
Feng group as they are since R1.

Reads, fixed before any result: read H = T1/E8.5->E9.5:cm, the read no counted candidate has run on. H is not a clean
read: the direction readout's smoke read Feng's real direction, E9.5 -> E14.5, within cardiomyocytes on exactly this
transition with the scorer's metric family (the committed experiments/t1_direction/results/smoke_E8.5_E9.5/per_type.csv,
type CM, predictor feng_E14.5, target full, held flat: sign agreement 0.717, Spearman 0.428, de_direction 0.296,
de_score 0.084 against a chance of 0.056). A pass on H is consistent with that known concordance and is no independent
evidence that a late Feng bracket transfers. The within-type E10.5 reads T1/E9.5->E10.5:cm and
T1/E9.5->E10.5:endocardium (configs/t1_late_proxy.yaml: target, reference and inputs cut to one chamber category, so a
gene class is not judged through E9.5's chamber mix) carry their own known feature, E10.5's earlier reads by the
readout and diag_v1. The three whole-population T1 tasks of configs/pseudo_tasks.yaml are reported against b0_scope and
the floor (they carry the composition confound the within-type reads were designed for, so they decide nothing).

The rule (vec.adoption's constants; fixed in the rows before any run; `gate_verdict` and `variant_verdict` apply it):
  gate (stage "gate"; T13 against b0_scope, with T14 and T15): lp_proxy passes the adoption rule's per-task test (paired
  mean gain above 0.3 points, ahead in at least 8 of 10 seeds) on H and on at least one within-type E10.5 read; on every
  read it passes, it also beats lp_cflip and lp_cperm by that test, and its margin over b0_scope sits on the rank
  metrics (`rank_share_ok`: the de_recovery and direction margins together above 0 and at least half the total margin;
  a gain mainly on mmd_u or the variogram is a fit and stops the variants); neither control itself passes the test on
  H; and lp_proxy loses more than 0.3 points on average on none of the three within-type reads. A within-type read that
  cannot be scored fails the whole gate, which is why the preflight checks them before any candidate is built.
  variants (stage "variants", run only if the gate passed; a variant that does not run does not count; T18 also needs
  N20's target-free stop to let it run): the variant, against its reference (lp_hyg and lp_eb against lp_proxy,
  lp_hyg_n against lp_hyg), has a mean gain above 0.3 and at least 9 of 10 seeds ahead on H and on BOTH within-type
  E10.5 reads, its margin sits on the rank metrics on each, and it loses more than 0.3 on average on none of the three
  within-type reads. Seeds check draw noise only (the ledger's opening); the multiplicity of the tries is carried by
  the try count, which every report states, and by asking for the pass on all three reads rather than one.
  What each outcome changes in M1's real direction is written in the tries' rows. No outcome changes an upload by
  itself: a change to M1's real direction lands only by the captain's ruling (no counted backtest can run it), and
  below N20's similarity floor a pass is platform evidence only, with no recommendation for M1.

Usage (DGX, repo root, a clean tree at the commit that registers the rows and carries X002's amendment; N20 has run at
that commit, and nothing is committed until the tournament's last decision is written; BLAS threads capped):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/late_proxy.py preflight
  ... late_proxy.py run --stage gate, ... late_proxy.py decide --stage gate; if it says the gate passed,
  ... late_proxy.py run --stage variants, ... late_proxy.py decide --stage variants; last, ... late_proxy.py register
  --stage gate and (if they ran) --stage variants. results/registry.csv is tracked, so `register` comes after the last
  run: a run on a tree with a change to a tracked file is refused here and at the decision. `run` takes --task KEY
  (repeatable; its shell-safe form too) to split a stage into parallel processes.

Writes under experiments/t1_m1/results/late_proxy_v1/ (nothing about a test board: the boards are pseudo-boards of
released stages): the preflight's verdict per read, one JSON per read and stage, as vec.harness.run_task returns it, and
the decision JSON and text. A result that exists is refused; each is stamped with the commit, vec.provenance's
clean-tree flag, the UTC time and the inputs' sha256. It prints and stores no test-board count, limit or anchor.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from vec import harness, provenance, registry
from vec.adoption import MAX_LOSS, MIN_AHEAD, MIN_GAIN, SEEDS, keep_change
from vec.boards import GROUPS, REPO_ROOT

CONFIG = REPO_ROOT / "configs" / "t1_late_proxy.yaml"
OUT = REPO_ROOT / "experiments" / "t1_m1" / "results" / "late_proxy_v1"
N20 = REPO_ROOT / "experiments" / "t1_m1" / "results" / "n20_proxy_profile" / "n20_proxy_profile.json"
TRIES = REPO_ROOT / "experiments" / "harness" / "TRIES.md"
PREFLIGHT = "preflight.json"
ROWS = {"gate": ("N20", "T13", "T14", "T15"), "variants": ("T16", "T17", "T18")}  # the rows each stage's run needs
BASE = "b0_scope"
PROXY = "lp_proxy"
CFLIP = "lp_cflip"
CPERM = "lp_cperm"
VARIANTS = {"lp_hyg": PROXY, "lp_hyg_n": "lp_hyg", "lp_eb": PROXY}  # each variant's reference (a nested ladder)
CANDIDATES = {"gate": (harness.FLOOR, "b0", BASE, PROXY, CFLIP, CPERM),  # b0: reported against B0 where it runs
              "variants": (harness.FLOOR, BASE, PROXY, *VARIANTS)}  # lp_proxy rerun: each variant's reference
READ_H = "T1/E8.5->E9.5:cm"  # no counted candidate has run on it; not clean (the module docstring's known feature)
WITHIN_E10 = ("T1/E9.5->E10.5:cm", "T1/E9.5->E10.5:endocardium")
WITHIN = (READ_H, *WITHIN_E10)
WHOLE = ("T1/E8.5->E9.5", "T1/E9.5->E10.5:chamber", "T1/E9.5->E10.5")
READS = (*WITHIN, *WHOLE)
VARIANT_AHEAD = 9  # the variants' seed bar, of 10. Seeds check draw noise only; it is not a multiplicity correction
MIN_RANK_SHARE = 0.5  # the rows' fit criterion: the de_recovery and direction margins' share of the total margin
LANDING = {"lp_hyg": "the first counted support for R2 (D, A1), which stands as built",
           "lp_hyg_n": "recommended to the captain: turn on R2's N1 and N2 lines for M1's real direction (his ruling; "
                       "T1's s1 and s2 rebuilt and hashed by 24 Nov)",
           "lp_eb": "recommended to the captain: the same feng_dirs.eb_shrink rule on M1's real direction behind a new "
                    "build option, default off (his ruling; T1's s1 and s2 rebuilt and hashed by 24 Nov)"}


def safe(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9.]+", "_", key).strip("_")


def within_tasks() -> list[harness.PseudoTask]:
    return harness.load_tasks(CONFIG)["T1:late_proxy"]


def whole_tasks() -> list[harness.PseudoTask]:
    return [t for t in harness.load_tasks()["T1:test"] if t.key in WHOLE]


def all_tasks() -> list[harness.PseudoTask]:
    tasks = whole_tasks() + within_tasks()
    if {t.key for t in tasks} != set(READS):
        raise ValueError(f"the reads are not the expected ones: {sorted(t.key for t in tasks)}")
    return tasks


def registered_rows() -> set[str]:
    """The row ids in experiments/harness/TRIES.md's tables (each table row's first word, where it is a row id)."""
    words = (ln[2:].split("|", 1)[0].split() for ln in TRIES.read_text(encoding="utf-8").splitlines()
             if ln[:2] == "| ")
    return {w[0] for w in words if w and re.fullmatch(r"[A-Z][0-9]+[a-z]?", w[0])}


def _prereqs(rows: Sequence[str], n20_path: Path) -> tuple[dict, str]:
    """What every stage needs before anything is read: its rows committed (TRIES.md rule 1), X002's amendment in the
    register (AGENTS.md), a clean tree and N20's yes, stamped clean at this very commit. Returns N20's report and the
    commit."""
    missing = sorted(set(rows) - registered_rows())
    if missing:
        sys.exit(f"rows {missing} are not in experiments/harness/TRIES.md: a row is committed before its run")
    fd = harness._feng_dirs()
    if not fd.x002_amended():
        sys.exit(f"EXTERNAL_DATA.md's X002 row does not carry the amendment (its phrase {fd.X002_AMENDMENT!r}): a use "
                 "of an external dataset is logged before it enters anything")
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    if unclean is None or unclean:
        sys.exit(f"the tree is not clean ({unclean}): a result from a dirty tree would be refused at the decision")
    if not n20_path.exists():
        sys.exit(f"{n20_path} is missing: N20 runs before the gate")
    n20 = json.loads(n20_path.read_text())
    if not n20.get("proceed_to_tries"):
        sys.exit("N20's gate said no: tries T13-T18 are not run and do not count")
    stamp_ = n20.get("stamp", {})
    if stamp_.get("commit") != commit or stamp_.get("tree_clean_at_end") is not True:
        sys.exit(f"N20 must have run clean at this commit ({commit}); its stamp says {stamp_}. Nothing is committed "
                 "between N20 and the tournament's last decision")
    return n20, commit


def check_ready(stage: str, n20_path: Path, out: Path) -> dict:
    """A row goes in before its run (TRIES.md rule 1); the gate runs only after N20's target-free gate said yes and the
    preflight passed at this commit, the variants only after the gate passed (decision_gate.json in `out`). Returns
    N20's report."""
    n20, commit = _prereqs(ROWS[stage], n20_path)
    if stage == "gate":
        pre = out / PREFLIGHT
        if not pre.exists():
            sys.exit(f"{pre} is missing: the preflight runs before the gate")
        p = json.loads(pre.read_text())
        if not p.get("passed"):
            sys.exit(f"the preflight found a read that cannot be scored ({p.get('reads')}): T13-T15 are not run and "
                     "do not count")
        p_stamp = p.get("stamp", {})
        if p_stamp.get("commit") != commit or p_stamp.get("tree_clean") is not True:
            sys.exit(f"the preflight must have passed clean at this commit ({commit}); its stamp says {p_stamp}")
    else:
        gate_file = out / "decision_gate.json"
        if not gate_file.exists() or not json.loads(gate_file.read_text())["verdict"]["passed"]:
            sys.exit(f"the variants run only after the gate passed ({gate_file}): T16-T18 are pruned, not counted")
    return n20


def variants_to_run(n20: Mapping) -> dict[str, str]:
    """T16-T18 and their references, less T18 unless N20's target-free stop let it run (a missing flag does not)."""
    return {v: ref for v, ref in VARIANTS.items() if v != "lp_eb" or n20.get("t18_runs") is True}


def stage_candidates(stage: str, n20: Mapping) -> tuple[str, ...]:
    """The candidates a stage's runs build (T18, lp_eb, leaves the variants when N20's stop says so)."""
    return CANDIDATES["gate"] if stage == "gate" else (harness.FLOOR, BASE, PROXY, *variants_to_run(n20))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 22):
            h.update(chunk)
    return h.hexdigest()


def stamp() -> dict:
    """The commit, vec.provenance's clean-tree flag (results/, runs/ and experiments/**/results/ exempt) and the UTC."""
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    return {"commit": commit, "tree_clean": None if unclean is None else not unclean,
            "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}


def input_hashes() -> dict:
    """The inputs these tries read outside the released stages (which run_task hashes per task): the reads' config,
    Feng's statistics and the module lists they hold flat (R2's lists, by sha256 as N12 records them)."""
    hygiene = harness._t1_hygiene()
    bundle = harness._late_proxy_load(harness.LP_STATS)
    return {"config": CONFIG.relative_to(REPO_ROOT).as_posix(), "config_sha256": sha256_file(CONFIG),
            "feng_stats": harness.LP_STATS.relative_to(REPO_ROOT).as_posix(),
            "feng_stats_sha256": harness._sha256(harness.LP_STATS), "bracket": list(harness.LP_BRACKET),
            "module_lists_sha256": {m: hygiene.list_sha256(syms) for m, syms in sorted(bundle["modules"].items())},
            "held_flat_named_sha256": hygiene.list_sha256(bundle["named"])}


# ---------------------------------------------------------------------------------------------------------------------
# The rule


def compare(results: Mapping[str, dict], x: str, y: str, reads: Sequence[str]) -> dict[str, dict]:
    """Candidate x against y, read by read, with vec.adoption's paired statistics (`keep_change`, over the 10 seeds):
    mean gain, seeds ahead, its test and veto, plus the variants' seed bar; an absent read or candidate is incomplete
    and passes nothing."""
    paired = {}
    for read in reads:
        r = results.get(read)
        if r is None or x not in r["candidates"] or y not in r["candidates"]:
            paired[read] = {}
            continue
        paired[read] = {row["seed"]: (row["scores"][y], row["scores"][x]) for row in r["rows"]
                        if x in row["scores"] and y in row["scores"]}
    out = {}
    for read, t in keep_change(paired)["tasks"].items():
        complete = t.get("complete", False)
        strict = bool(complete and t["mean_gain"] > MIN_GAIN and t["seeds_ahead"] >= VARIANT_AHEAD)
        out[read] = {"complete": complete, "mean_gain": t.get("mean_gain"), "seeds_ahead": t.get("seeds_ahead"),
                     "passes": bool(t.get("passes")), "passes_strict": strict, "vetoes": bool(t.get("vetoes"))}
        if read in results and complete:
            out[read]["group_margins"] = group_margins(results[read], x, y)
    return out


def group_margins(result: dict, x: str, y: str) -> dict[str, float]:
    """Per scored T1 group (vec.boards.GROUPS: de_recovery, direction, distribution, covariance), x's mean margin over
    y in board points (100 x weight x the group skill difference; their sum is the score margin), so a row's "the gain
    sits on the metrics its reason predicts" can be read."""
    rows = [r for r in result["rows"] if x in r.get("groups", {}) and y in r.get("groups", {})]
    return {name: float(np.mean([100.0 * w * (r["groups"][x][name] - r["groups"][y][name]) for r in rows]))
            for name, w, _ in GROUPS["T1"]} if rows else {}


def rank_share_ok(t: Mapping) -> bool:
    """TRIES.md rule 2 and the rows' fit column as a condition of a pass, not a reading after it: on a read that counts,
    the de_recovery and direction margins together are above 0 and at least MIN_RANK_SHARE of the total margin over the
    four scored groups (a gain mainly on mmd_u or the variogram is a fit). An incomplete read counts as false."""
    m = t.get("group_margins")
    if not t.get("complete") or not m:
        return False
    rank = m["de_recovery"] + m["direction"]
    return rank > 0 and rank >= MIN_RANK_SHARE * sum(m.values())


def gate_verdict(results: Mapping[str, dict]) -> dict:
    """The gate: T13 against b0_scope and its controls. Returns the per-read comparisons and whether each condition of
    the rule (see the module docstring) holds."""
    proxy_base = compare(results, PROXY, BASE, READS)
    proxy_flip = compare(results, PROXY, CFLIP, WITHIN)
    proxy_perm = compare(results, PROXY, CPERM, WITHIN)
    perm_base = compare(results, CPERM, BASE, (READ_H,))
    flip_base = compare(results, CFLIP, BASE, (READ_H,))
    passed_w = [r for r in WITHIN_E10 if proxy_base[r]["passes"]]
    checked = ([READ_H] if proxy_base[READ_H]["passes"] else []) + passed_w
    conditions = {
        "read_h": proxy_base[READ_H]["passes"],
        "within_type_e10": bool(passed_w),
        "beats_cflip": bool(checked) and all(proxy_flip[r]["passes"] for r in checked),
        "beats_cperm": bool(checked) and all(proxy_perm[r]["passes"] for r in checked),
        "gain_on_rank_metrics": bool(checked) and all(rank_share_ok(proxy_base[r]) for r in checked),
        "controls_do_not_pass_h": not (perm_base[READ_H]["passes"] or flip_base[READ_H]["passes"]),
        "no_veto": all(proxy_base[r]["complete"] and not proxy_base[r]["vetoes"] for r in WITHIN)}
    return {"passed": all(conditions.values()), "conditions": conditions, "passed_within_e10": passed_w,
            "comparisons": {f"{PROXY} vs {BASE}": proxy_base, f"{PROXY} vs {CFLIP}": proxy_flip,
                            f"{PROXY} vs {CPERM}": proxy_perm, f"{CPERM} vs {BASE}": perm_base,
                            f"{CFLIP} vs {BASE}": flip_base,
                            f"{PROXY} vs b0 (reported)": compare(results, PROXY, "b0", READS),
                            f"{PROXY} vs the floor (reported)": compare(results, PROXY, harness.FLOOR, READS)}}


def variant_verdict(results: Mapping[str, dict], variant: str, reference: str) -> dict:
    """A variant against its reference: paired mean gain above 0.3 and at least 9 of 10 seeds ahead on H and on both
    within-type E10.5 reads, its margin on the rank metrics on each, and no mean loss above 0.3 on any within-type
    read. The whole-population reads are in the table as a report and decide nothing."""
    table = compare(results, variant, reference, READS)
    passing = [r for r in WITHIN if table[r]["passes_strict"]]
    passed_w = [r for r in WITHIN_E10 if table[r]["passes_strict"]]
    conditions = {"read_h": table[READ_H]["passes_strict"], "both_within_type_e10": len(passed_w) == len(WITHIN_E10),
                  "gain_on_rank_metrics": bool(passing) and all(rank_share_ok(table[r]) for r in passing),
                  "no_veto": all(table[r]["complete"] and not table[r]["vetoes"] for r in WITHIN)}
    return {"variant": variant, "reference": reference, "passed": all(conditions.values()), "conditions": conditions,
            "passed_within_e10": passed_w, "comparison": table}


def landing(v: Mapping, close: bool) -> str:
    """What a variant's outcome may say about M1's real direction (its row's landing, fixed before any run): below
    N20's similarity floor a pass is platform evidence only, with no recommendation for M1."""
    if not v["passed"]:
        return "not passed: nothing lands"
    if not close:
        return "passed, but the stand-in is below N20's similarity floor: platform evidence only, no recommendation " \
               "for M1"
    return f"passed: {LANDING[v['variant']]}"


# ---------------------------------------------------------------------------------------------------------------------
# Report


def _fmt(t: Mapping) -> str:
    if not t.get("complete"):
        return "incomplete"
    m = "; ".join(f"{g} {v:+.2f}" for g, v in t.get("group_margins", {}).items())
    return f"{t['mean_gain']:+.2f} ({t['seeds_ahead']}/10){' [' + m + ']' if m else ''}"


def report_text(stage: str, verdict: Mapping, stamp_: Mapping, inputs: Mapping) -> str:
    n20 = verdict["n20"]
    lines = [f"T1 late-proxy tournament, stage {stage} (experiments/harness/TRIES.md T13-T18)",
             f"Commit {stamp_['commit']}, tree clean {stamp_['tree_clean']}, UTC {stamp_['utc']}",
             f"Inputs: {json.dumps(inputs, sort_keys=True)}",
             f"N20: stand-in close to the real direction: {n20['proxy_close_to_real']} ({n20['landing_floor']}); "
             f"T18 allowed to run: {n20['t18_runs']}", ""]
    if stage == "gate":
        lines.append(f"GATE: {'PASSED' if verdict['passed'] else 'FAILED'}; conditions {verdict['conditions']}")
        for name, table in verdict["comparisons"].items():
            lines.append(f"  {name}")
            lines += [f"    {r}: {_fmt(t)}" for r, t in table.items()]
    else:
        for v in verdict["variants"]:
            lines.append(f"{v['variant']} against {v['reference']}: {'PASSED' if v['passed'] else 'not passed'}; "
                         f"conditions {v['conditions']}")
            lines.append(f"  landing: {v['landing']}")
            lines += [f"    {r}: {_fmt(t)}" for r, t in v["comparison"].items()]
        lines += [f"{v} not run: {why}" for v, why in verdict["not_run"].items()]
    return "\n".join(lines) + "\n"


def load_results(folder: Path, reads: Sequence[str]) -> dict[str, dict]:
    """The stage's result JSONs by read key; a read with no file (a failed run) is simply absent and passes nothing."""
    out = {}
    for t in reads:
        p = folder / f"{safe(t)}.json"
        if p.exists():
            out[t] = json.loads(p.read_text())
    return out


def check_consistent(results: Mapping[str, dict], tasks: Sequence[harness.PseudoTask]) -> None:
    """vec.harness's refusals for one decision: one commit, a clean tree, the same config, protocol, seeds and hashes.
    The reads' own inputs (Feng's statistics, the module lists, the reads' config) must also be the same throughout."""
    harness.check_results(list(results.values()), tasks)
    if len({json.dumps(r["late_proxy"]["inputs"], sort_keys=True) for r in results.values()}) != 1:
        raise harness.InconsistentResults("reads with different inputs (Feng's statistics, module lists or config)")


def registry_rows(results: Mapping[str, dict], folder: Path, stage: str) -> list[dict]:
    """One results/registry.csv row per read and candidate (AGENTS.md, "Evidence over claims"), as run.py writes."""
    def rel(p: Path) -> str:
        return str(p.relative_to(REPO_ROOT)) if p.is_relative_to(REPO_ROOT) else str(p)

    rows = []
    for key, r in results.items():
        for n in r["candidates"]:
            s = np.array([row["scores"][n] for row in r["rows"] if n in row["scores"]], float)
            rows.append({"utc_time": r["provenance"]["utc"], "commit": r["provenance"]["commit"],
                         "tree_clean": r["provenance"]["tree_clean"], "source": "experiments/t1_m1/late_proxy.py",
                         "run_id": f"late_proxy_v1/{stage}", "test_board": r["task"]["test_board"], "pseudo_task": key,
                         "scoring_board": r["board"]["key"], "candidate": n,
                         "params": json.dumps(r["params"].get(n, {})),
                         "seeds": ",".join(map(str, r["seeds"])), "score_mean": f"{s.mean():.3f}",
                         "score_sd": f"{s.std(ddof=1) if len(s) > 1 else 0.0:.3f}", "protocol": r["protocol"],
                         "build": r["build"], "anchors": r["anchors"], "gate": "floor_normalised", "frac": r["frac"],
                         "scores_json": rel(folder / f"{safe(key)}.json")})
    return rows


# ---------------------------------------------------------------------------------------------------------------------
# Commands


def cmd_preflight(a) -> None:
    """The three within-type reads' floor and anchors on all ten seeds, before any candidate is built: a read that
    cannot be scored (for example inverted anchors on a small cut) would fail the whole gate and spend T13-T15. It
    builds no candidate and is no try; a failure stops the tournament here, uncounted."""
    _prereqs(ROWS["gate"], Path(a.n20))
    out = Path(a.out)
    path = out / PREFLIGHT
    if path.exists():
        sys.exit(f"{path} exists: the preflight is never overwritten (a rerun is a new row)")
    reads = {}
    for t in within_tasks():
        try:
            harness.run_task(t, [], seeds=tuple(SEEDS), data_dir=Path(a.data), build=harness.BUILD, log=None)
            reads[t.key] = "ok"
        except Exception as e:  # noqa: BLE001 (recorded; the tournament then stops here)
            reads[t.key] = f"{type(e).__name__}: {e}"
    report = {"passed": all(v == "ok" for v in reads.values()), "reads": reads, "stamp": stamp()}
    out.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=1) + "\n")
    print(f"preflight: {'PASSED' if report['passed'] else 'FAILED'}; {reads}")
    if not report["passed"]:
        sys.exit("a within-type read cannot be scored: T13-T15 are not run and do not count")


def cmd_run(a) -> None:
    n20 = check_ready(a.stage, Path(a.n20), Path(a.out))
    tasks = all_tasks()
    chosen = [t for t in tasks if not a.task or t.key in a.task or safe(t.key) in a.task]
    if a.task and len(chosen) != len(set(a.task)):
        sys.exit(f"unknown reads among {a.task}; known: {[t.key for t in tasks]}")
    out = Path(a.out) / a.stage
    out.mkdir(parents=True, exist_ok=True)
    cands = [harness.CANDIDATES[n] for n in stage_candidates(a.stage, n20)]
    inputs = input_hashes()
    failed = []
    for t in chosen:
        path = out / f"{safe(t.key)}.json"
        if path.exists():
            sys.exit(f"{path} exists: a late-proxy result is never overwritten (a rerun is a new row)")
        t0 = time.time()
        try:
            res = harness.run_task(t, cands, seeds=tuple(SEEDS), data_dir=Path(a.data), build=harness.BUILD)
        except Exception as e:  # noqa: BLE001 (recorded; the read then passes nothing in the decision)
            failed.append(t.key)
            (out / f"{safe(t.key)}.failed.txt").write_text(f"{type(e).__name__}: {e}\n")
            print(f"FAILED {t.key}: {type(e).__name__}: {e}", flush=True)
            continue
        res["late_proxy"] = {"stage": a.stage, "stamp": stamp(), "inputs": inputs,
                             "seconds": round(time.time() - t0, 1)}
        partial = path.with_suffix(".json.partial")
        partial.write_text(json.dumps(res, indent=1, default=float))
        partial.replace(path)
        print(f"wrote {path} ({res['late_proxy']['seconds']} s)", flush=True)
    if failed:
        sys.exit(f"{len(failed)} read(s) failed: {failed} (the decision treats them as not passing)")


def cmd_decide(a) -> None:
    out = Path(a.out)
    json_path, text_path = out / f"decision_{a.stage}.json", out / f"decision_{a.stage}.txt"
    if json_path.exists() or text_path.exists():
        sys.exit(f"{json_path} exists: a late-proxy decision is never overwritten")
    tasks = all_tasks()
    results = load_results(Path(a.runs) / a.stage, READS)
    if not results:
        sys.exit(f"no results for stage {a.stage} under {a.runs}")
    if not Path(a.n20).exists():
        sys.exit(f"{a.n20} is missing: the decision records N20's similarity floor and T18's stop")
    n20 = json.loads(Path(a.n20).read_text())
    check_consistent(results, tasks)
    commit = next(iter(results.values()))["provenance"]["commit"]
    if a.stage == "gate":
        verdict = gate_verdict(results)
    else:
        gate_file = out / "decision_gate.json"
        if not gate_file.exists():
            sys.exit("no decision_gate.json: the variants are read only after the gate is decided")
        gate = json.loads(gate_file.read_text())
        if not gate["verdict"]["passed"]:
            sys.exit("the gate did not pass: the variants are not run or read (T16-T18 are pruned, not counted)")
        if gate["stamp"]["reads_commit"] != commit:
            sys.exit(f"the gate's reads ran at {gate['stamp']['reads_commit']} and the variants' at {commit}: "
                     "one commit")
        reference_scores = {r: (results[r]["rows"][0]["scores"].get(PROXY) if r in results else None) for r in READS}
        gate_runs = load_results(Path(a.runs) / "gate", READS)
        for r, v in reference_scores.items():
            if v is None or r not in gate_runs:
                continue
            g = gate_runs[r]["rows"][0]["scores"].get(PROXY)
            if g is not None and abs(float(g) - float(v)) > 1e-9:
                sys.exit(f"{r}: lp_proxy scored {g} in the gate run and {v} in the variants run (seed 0): "
                         "not reproducible")
        active = variants_to_run(n20)
        verdict = {"variants": [variant_verdict(results, v, ref) for v, ref in active.items()],
                   "not_run": {v: "N20's target-free stop: its factors barely change the direction on the reads' "
                                  "groups, so it was not run and does not count" for v in VARIANTS if v not in active}}
        close = n20.get("proxy_close_to_real") is True
        for v in verdict["variants"]:
            v["landing"] = landing(v, close)
    verdict["n20"] = {"proxy_close_to_real": n20.get("proxy_close_to_real") is True,
                      "landing_floor": n20.get("landing_floor"), "t18_runs": n20.get("t18_runs") is True}
    inputs = next(iter(results.values()))["late_proxy"]["inputs"]  # the same for every read (check_consistent)
    stamp_ = {**stamp(), "reads_commit": commit}
    report = {"stage": a.stage, "verdict": verdict, "stamp": stamp_, "inputs": inputs,
              "constants": {"min_gain": MIN_GAIN, "min_ahead": MIN_AHEAD, "variant_ahead": VARIANT_AHEAD,
                            "max_loss": MAX_LOSS, "min_rank_share": MIN_RANK_SHARE, "read_h": READ_H,
                            "within_e10": list(WITHIN_E10), "whole": list(WHOLE), "variants": VARIANTS, "base": BASE}}
    out.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(report, indent=1, default=float) + "\n")
    text = report_text(a.stage, verdict, stamp_, inputs)
    text_path.write_text(text)
    print(text, end="")


def cmd_register(a) -> None:
    """Append a decided stage's rows to the registry, once. results/registry.csv is tracked, so a change to it dirties
    the tree: this runs after the tournament's last run and decision, never between the gate and the variants."""
    runs, out, path = Path(a.runs), Path(a.out), Path(a.registry)
    decision = out / f"decision_{a.stage}.json"
    if not decision.exists():
        sys.exit(f"{decision} is missing: a stage is registered after its decision")
    results = load_results(runs / a.stage, READS)
    if not results:
        sys.exit(f"no results for stage {a.stage} under {runs}")
    check_consistent(results, all_tasks())
    run_id = f"late_proxy_v1/{a.stage}"
    if path.exists():
        with path.open(newline="") as f:
            if any(row.get("run_id") == run_id for row in csv.DictReader(f)):
                sys.exit(f"{path} already holds rows of {run_id}: a stage is registered once")
    n = registry.append(registry_rows(results, runs / a.stage, a.stage), path)
    print(f"registered {n} rows of {run_id} in {path}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("preflight")
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p.add_argument("--out", default=str(OUT))
    p.add_argument("--n20", default=str(N20), help="N20's report (its gate must have said proceed)")
    r = sub.add_parser("run")
    r.add_argument("--stage", required=True, choices=("gate", "variants"))
    r.add_argument("--task", action="append", help="read key (repeatable, or its shell-safe form); default: all six")
    r.add_argument("--data", default=str(harness.DATA_DIR))
    r.add_argument("--out", default=str(OUT))
    r.add_argument("--n20", default=str(N20), help="N20's report (its gate must have said proceed)")
    d = sub.add_parser("decide")
    d.add_argument("--stage", required=True, choices=("gate", "variants"))
    d.add_argument("--runs", default=str(OUT), help="the folder of `run`'s results")
    d.add_argument("--out", default=str(OUT))
    d.add_argument("--n20", default=str(N20), help="N20's report (its similarity floor and T18's stop are recorded)")
    g = sub.add_parser("register")
    g.add_argument("--stage", required=True, choices=("gate", "variants"))
    g.add_argument("--runs", default=str(OUT), help="the folder of `run`'s results")
    g.add_argument("--out", default=str(OUT), help="the folder of the stage's decision")
    g.add_argument("--registry", default=str(registry.REGISTRY))
    a = ap.parse_args(argv)
    {"preflight": cmd_preflight, "run": cmd_run, "decide": cmd_decide, "register": cmd_register}[a.cmd](a)


if __name__ == "__main__":
    main()
