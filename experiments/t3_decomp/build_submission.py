#!/usr/bin/env python
"""Task 3 submission builder: the pre-registered decision file in, a submission file out (docs/PLAN.md section 4.4).

`t3_decompose.py gata4` (release day) writes gata4_decisions.json: for each organiser answer, which components enter
submissions 1 and 2 and with which weights, plus the per-gene factors under <run>/factors/. This script re-decides
nothing. It reads the named submission from the decision file, loads the factor file the decision run wrote (checked
against the sha256 recorded with it) and builds the prediction with the functions the scored variants were built
with: `t3_decompose.build` resamples real WT cells (vec.prep.load_stage(drop_duplicates=True)) to the composition,
applies the factor in count space (binomial thinning below 1, Poisson top-up from the kNN-smoothed profiles above 1)
and the per-label depth ratios, and mixes in stage-pool cells when stage mixing is on; the composition is
`t3_decompose.rule_comp` (scenario states depleted by 2^-log2_magnitude) and `t3_decompose.coverage_comp`. Rows are then
put in random order, each cell keeps its source WT cell's coordinates (the T3 board requires them; they are not
ranked), and the file is written with vec.io.make_submission and checked with vec.io.check_submission against the T3
board. A JSON sidecar next to it records the decision file's sha256, the branch, the submission, the hashes of this
script and of the pre-registered files, and every input.

Branches (organiser follow-up question 3: may the released Gata4 files be used inside the beta-catenin prediction?):
  allowed  submissions.allowed of the decision file, with rules.yaml note 8 (the user, 2026-10-19): when the factor
           table reaches its last row (none of D2, D2u and D10 passes), s1 = the Mab21l2 factor if D11b passes, else
           WT cells, and s2 = the Gata4 factor at exponent 0.5 (see `note8`); then note 6 (lead, 2026-10-19): unless
           D11b passes, s1 and s2 must not both use the Gata4 factor (s2 becomes WT cells with the central composition
           when D7 passes, mild (m/2) when D8 fires (note 7), otherwise the Gata4 factor at exponent 0.5 with s1's
           composition; see `note6`)
  refused  submissions.refused (training on the Gata4 files refused, scoring against them allowed)
  barred   rules.yaml note 3 (lead, 2026-10-19): the Gata4 files may not be used at all, not even to score candidates.
           Both submissions use the factor source mab21l2 (transfer:mab21l2_all) at exponent 1; s1 has the central
           composition at m = D7m, s2 the mild one (m/2).
A decision file without a submissions table (the Mab21l2 smoke run) gets its table from t3_decompose.decide_submissions
applied to its own rule outcomes, and the sidecar says so.

Refusals (exit status 2): an unknown branch or submission; a decision file whose script or rules sha256 differs from
t3_decompose.py / rules.yaml on disk (--allow-hash-mismatch overrides and is recorded); a decision file from a
non-registered protocol (--allow-unregistered-protocol); a stored submissions table that the rules on disk do not
reproduce; a factor file whose sha256 differs from its record; a WT or stage file that differs from the one the run
recorded; a board that is not T3; a cell count outside the board's limits; a build seed below
vec.harness.BUILD_SEED_OFFSET (docs/P3_RUNBOOK.md section 5: never seed 0).

Commands (DGX, repository root; competition data and built files stay there):

  uv run python experiments/t3_decomp/build_submission.py build DECISIONS.json --branch allowed --submission s1 \\
      --wt data/raw/T2/E8.75.h5ad --run-dir runs/t3_decomp/gata4_<STAMP> --board <T3 board key> \\
      --seed 1000001 --out runs/t3_build/<name>.h5ad
  (--n-cells defaults to "max", the board's max_cells read in code and never printed: the approved count)
  ... plan DECISIONS.json [--out FILE]     both submissions of every branch as `build` would make them, with the
                                           notes applied; reads only the decision file
  ... factors --preset mab21l2 --out DIR   per-gene factors written by t3_decompose's own code for a decision file
                                           whose run wrote none (the Mab21l2 smoke test); --preset gata4 writes the
                                           Mab21l2 transfer factor and R on WT E8.75 without any Gata4 file
  ... replay --out FILE                    rebuilds scored Mab21l2 variants through this builder and compares them
                                           with the run's records seed by seed
  ... score FILE --target KO --ref WT --board KEY [--records CSV --variant NAME]
                                           vec.evaluate scores of a built file over seeds, against a recorded variant
"""

from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import time
from contextlib import contextmanager
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import typer
import yaml
from rich.console import Console
from rich.markup import escape
from rich.table import Table

import vec.provenance  # noqa: F401  (records this build's start, before the method code below is loaded)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
import t3_decompose as t3d  # noqa: E402

BUILDER_PATH = Path(__file__).resolve()
T3D_PATH = HERE / "t3_decompose.py"
RULES_PATH = HERE / "rules.yaml"
DATA = REPO / "data" / "raw"
BRANCHES = ("allowed", "refused", "barred")
SUBMISSIONS = ("s1", "s2")
# factor source (decision tables) -> file under <run>/factors/ ({preset}: the decision file's own knockout(s))
FACTOR_FILES = {"gata4": "{preset}_pooled", "gata4_unvalidated": "{preset}_pooled", "reference": "reference_R",
                "mab21l2": "transfer_mab21l2_all"}
GATA4_SOURCES = ("gata4", "gata4_unvalidated")          # factor sources that are the pooled Gata4 factor
COMPOSITION_KINDS = ("none", "central", "mild", "mild_raw")
EXTRAS = ("depth", "stage", "coverage")
# The captain's ruling of 2026-11-13 (docs/PLAN.md decision log 2026-11-13 (d); experiments/harness/TRIES.md row R5):
# T3's second upload, the refused branch's s1, also carries the wild-type-only factor of tries K7b and K8b, applied per
# lineage group at the tried exponent. It is fitted on released wild-type stage files alone (no knockout file), by the
# tried estimator itself (k7_wt_revert). --wt-revert switches it on; without the flag every file builds as before.
WT_REVERT_SLOT = ("refused", "s1")
WT_REVERT_STREAM = "wt_revert"
JITTER_NN = 0.25          # duplicate positions are moved by at most this x the WT median nearest-neighbour distance
COORD_KEY = "spatial_3D"
# scored Mab21l2 variants that `replay` rebuilds: which components each one switches on (rules.yaml components)
REPLAY_VARIANTS = {
    "fac:none": {},
    "fac:factor": {"factor": True},
    "fac:depth+factor": {"factor": True, "depth": True},
    "fac:stage+factor": {"factor": True, "stage": True},
    "fac[top]:cov+factor": {"factor": True, "coverage": True},
}

app = typer.Typer(no_args_is_help=True, add_completion=False)
console = Console(emoji=False, highlight=False)


class Refused(Exception):
    """The builder will not build; the message says why (the command exits with status 2)."""


def say(msg: str, style: str | None = None) -> None:
    console.print(f"[{style}]{escape(msg)}[/]" if style else escape(msg), soft_wrap=True)


def sha256(path) -> str:
    return t3d.sha256(path)


def load_rules() -> dict:
    return yaml.safe_load(RULES_PATH.read_text())


def dump_json(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(t3d.json_safe(obj), indent=1))


def roundtrip(obj):
    """obj as it reads back from a decision file (JSON types)."""
    return json.loads(json.dumps(t3d.json_safe(obj)))


# ---------------------------------------------------------------------------------------------------------------------
# Decision file: hashes, the named submission


def load_decision(path: Path, rules: dict) -> dict:
    try:
        dec = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as e:
        raise Refused(f"cannot read the decision file {path}: {e}") from e
    missing = [k for k in ("preset", "rules_sha256", "script_sha256", "rules") if k not in dec]
    if missing:
        raise Refused(f"{path} is not a t3_decompose decision file (missing {missing})")
    if dec["preset"] not in rules["presets"]:
        raise Refused(f"{path}: unknown preset {dec['preset']!r}")
    return dec


def check_hashes(dec: dict, allow_mismatch: bool) -> dict:
    """The decision file's script and rules sha256 against t3_decompose.py and rules.yaml on disk."""
    disk = {"script": sha256(T3D_PATH), "rules": sha256(RULES_PATH)}
    rec = {"t3_decompose_on_disk": disk["script"], "rules_on_disk": disk["rules"],
           "decision_script_sha256": dec["script_sha256"], "decision_rules_sha256": dec["rules_sha256"],
           "mismatch": [k for k in ("script", "rules") if dec[f"{k}_sha256"] != disk[k]],
           "allow_hash_mismatch": bool(allow_mismatch)}
    if rec["mismatch"] and not allow_mismatch:
        raise Refused(f"the decision file's {' and '.join(rec['mismatch'])} sha256 differs from the file on disk "
                      "(it was made by other code or rules); pass --allow-hash-mismatch to build anyway (recorded)")
    return rec


def composition_spec(kind: str, m: float, rules: dict, source: str) -> dict:
    """A composition in decide_submissions' form: the scenario states depleted by 2^-(m x weight)."""
    if kind == "none":
        return {"kind": "none"}
    scen = rules["presets"]["gata4"]["scenario"]
    l2 = float((m or 0.0) * rules["weights"]["composition"][kind])
    return {"kind": kind, "depleted_states": scen["depleted_states"], "spared_states": scen["spared_states"],
            "log2_magnitude": l2, "depletion_factor": float(2.0 ** -l2), "magnitude_source": source}


def barred_table(dec: dict, rules: dict) -> dict:
    """rules.yaml note 3 (lead, 2026-10-19): both submissions use the factor source mab21l2 at exponent 1; s1 central
    composition at m = D7m, s2 mild (m/2). If m is 0 the two files would be the same, so s2 takes half of s1's
    exponent (submissions.identical_fallback)."""
    d7m = dec["rules"].get("D7m") or {}
    m = float(d7m.get("m") or 0.0)
    src = "Mab21l2 (D7m)"
    s1 = {"factor": {"source": "mab21l2", "exponent": 1.0}, "composition": composition_spec("central", m, rules, src)}
    s2 = {"factor": {"source": "mab21l2", "exponent": 1.0}, "composition": composition_spec("mild", m, rules, src)}
    note = None
    if m <= 0:
        s2["factor"]["exponent"] = 0.5
        note = "D7m gives m = 0, so both compositions are unchanged WT: s2 takes factor exponent 0.5"
    warn = None if d7m.get("outcome") else "D7m did not pass (no counted state): no composition change"
    return {"factor_row": {"when": "rules.yaml note 3", "s1": "mab21l2", "s2": "mab21l2"}, "s1": s1, "s2": s2,
            "note": note, "warning": warn}


def note6(dec: dict, rules: dict, tab: dict) -> tuple[dict, dict]:
    """rules.yaml note 6 (lead, 2026-10-19, before the Gata4 download), allowed branch: if D11b (the Mab21l2 factor
    scored on Gata4) passes, the allowed tables apply as written. If it fails or could not be evaluated, the two
    submissions must not both use the Gata4 factor: s1 keeps the table's choice; s2 = no factor (WT cells) with the
    central composition at m when D7 passes (mild, m/2, when D8 fires: note 7), otherwise the Gata4 factor at
    exponent 0.5 with s1's composition.
    t3_decompose.decide_submissions predates the note, so it is applied here. Reading: the note changes s2 when the
    table's s1 and s2 both use the Gata4 factor; when one of them does not, the note's condition already holds and the
    table stands. The replacement carries no depth, stage or coverage, and if its factor would equal s1's (the same
    source at the same exponent) it takes half of s1's exponent (note 7, in the direction of the hedge). Returns s2
    and the record (with D11b's per-target scores)."""
    d11b = dec["rules"].get("D11b") or {}
    s1, s2 = tab["s1"], tab["s2"]
    rec = {"D11b_outcome": d11b.get("outcome"), "D11b_per_target": d11b.get("per_target"), "table_s2": s2}
    if d11b.get("outcome") is True:
        return s2, {**rec, "applies": False, "reason": "D11b passed: note 6 changes nothing"}
    if not (s1["factor"]["source"] in GATA4_SOURCES and s2["factor"]["source"] in GATA4_SOURCES):
        return s2, {**rec, "applies": False, "reason": "D11b did not pass, but the table's s1 and s2 do not both use "
                                                       "the Gata4 factor, so the note's condition holds as written"}
    d7, d8 = dec["rules"].get("D7") or {}, dec["rules"].get("D8") or {}
    if d7.get("outcome"):
        kind = "mild" if d8.get("outcome") is True else "central"   # note 7: keep the D8 hedge on the composition
        new = {"factor": {"source": "none", "exponent": 0.0},
               "composition": composition_spec(kind, d7.get("m"), rules, "Gata4 replicates (D7)")}
        how = f"D11b did not pass and D7 passed: s2 = no factor (WT cells) with the {kind} composition" + \
              (" (D8: the composition sign is a live risk)" if kind == "mild" else " at m")
    else:
        new = {"factor": {"source": s1["factor"]["source"], "exponent": 0.5},
               "composition": json.loads(json.dumps(s1["composition"]))}
        how = "D11b and D7 did not pass: s2 = the Gata4 factor at exponent 0.5 with s1's composition"
        if new["factor"] == s1["factor"]:    # note 7: compared on the factor alone, since the replacement has no extras
            new["factor"]["exponent"] = 0.5 * float(s1["factor"]["exponent"])
            how += f"; that is s1's factor, so s2 takes half of s1's exponent ({new['factor']['exponent']})"
    new.update({k: False for k in EXTRAS if k in s1 or k in s2})     # the allowed tables' form; no extras (note 7)
    return new, {**rec, "applies": True, "reason": how}


def rule_outcomes(dec: dict) -> dict:
    """Rule outcomes under the names the decision tables use (as t3_decompose.decide_submissions reads them)."""
    return {k: v.get("outcome") for k, v in dec["rules"].items()
            if isinstance(v, dict) and "outcome" in v and not isinstance(v.get("outcome"), dict)}


def note8(dec: dict, rules: dict, tab: dict) -> tuple[dict, dict]:
    """rules.yaml note 8 (the user's decision, 2026-10-19, before the Gata4 download), allowed branch: when the factor
    table reaches its last row (none of D2, D2u and D10 passes), no check supports the Gata4 factor, so s1 does not bet
    on it: s1 = the Mab21l2 factor (exponent 1) if D11b passes, otherwise no factor (WT cells), with composition_s1;
    s2 = the pooled Gata4 factor at exponent 0.5 without extras, with composition_s2_other_factor. The two never
    coincide, and note 6 then changes nothing. Applied before note 6. Returns the table (replaced, or as it was) and
    the record."""
    b = rules["submissions"]["allowed"]
    rec = {"table_factor_row": tab.get("factor_row")}
    if tab.get("factor_row") != b["factor_table"][-1]:
        return tab, {**rec, "applies": False, "reason": "the factor table did not reach its last row"}
    names = rule_outcomes(dec)
    d7 = dec["rules"].get("D7") or {}

    def comp(table: str) -> dict:
        kind = t3d.first_match(b[table], names)["value"]
        return composition_spec(kind, d7.get("m_raw" if kind == "mild_raw" else "m"), rules, "Gata4 replicates (D7)")

    off = {k: False for k in EXTRAS if k in tab["s1"] or k in tab["s2"]}
    mab = names.get("D11b") is True
    s1 = {"factor": {"source": "mab21l2", "exponent": 1.0} if mab else {"source": "none", "exponent": 0.0},
          "composition": comp("composition_s1"), **off}
    s2 = {"factor": {"source": "gata4_unvalidated", "exponent": 0.5},
          "composition": comp("composition_s2_other_factor"), **off}
    first = "the Mab21l2 factor (D11b passed)" if mab else "no factor (WT cells)"
    how = (f"none of D2, D2u and D10 passed: s1 = {first} with the {s1['composition']['kind']} composition; s2 = the "
           f"Gata4 factor at exponent 0.5 with the {s2['composition']['kind']} composition")
    plain = not mab and s1["composition"]["kind"] == "none"
    new = {**tab, "s1": s1, "s2": s2, "note": None,
           "warning": "no validated component: s1 is plain WT resampling; review before any upload" if plain else None}
    return new, {**rec, "applies": True, "reason": how, "D11b_outcome": names.get("D11b"), "table_s1": tab["s1"],
                 "table_s2": tab["s2"], "table_note": tab.get("note"), "table_warning": tab.get("warning")}


def resolve_spec(dec: dict, rules: dict, branch: str, submission: str, hashes: dict) -> tuple[dict, dict]:
    """The submission the decision file names for this branch (with rules.yaml note 6 on the allowed branch), and
    where it came from."""
    if branch not in BRANCHES:
        raise Refused(f"unknown branch {branch!r}; expected one of {list(BRANCHES)}")
    if submission not in SUBMISSIONS:
        raise Refused(f"unknown submission {submission!r}; expected one of {list(SUBMISSIONS)}")
    if branch == "barred":
        tab = barred_table(dec, rules)
        return tab[submission], {"source": "rules.yaml note 3 (barred branch)", "factor_row": tab["factor_row"],
                                 "note": tab["note"], "warning": tab["warning"]}
    stored = dec.get("submissions")
    params = stored.get("parameters", {}) if isinstance(stored, dict) else {}
    derived = roundtrip(t3d.decide_submissions(rules, dec["rules"], None, {"submission_parameters": params}))
    if isinstance(stored, dict):
        if branch not in stored:
            raise Refused(f"the decision file has no submissions.{branch} table")
        tab = stored[branch]
        keys = ("factor_row", "s1", "s2", "note", "warning")
        same = all(derived[branch].get(k) == tab.get(k) for k in keys)
        if not same and not hashes["mismatch"]:
            raise Refused(f"submissions.{branch} in the decision file is not what rules.yaml gives for the file's rule "
                          "outcomes (t3_decompose.decide_submissions); the file was changed after the run")
        info = {"source": f"decision file submissions.{branch}.{submission}", "recomputed_matches": same}
    else:
        if dec["preset"] == "gata4":
            raise Refused("a gata4 decision file without a submissions table")
        tab = derived[branch]
        info = {"source": "derived: t3_decompose.decide_submissions on the decision file's rule outcomes (this "
                          f"{dec['preset']} decision file carries no submissions table)", "recomputed_matches": None}
    info.update(factor_row=tab.get("factor_row"), note=tab.get("note"), warning=tab.get("warning"))
    if branch == "allowed":
        tab, info["note8"] = note8(dec, rules, tab)
        s2, info["note6"] = note6(dec, rules, tab)
        info["warning"] = tab.get("warning")
        # once a note changed the pair, the table's note describes a pair that is not built (it is kept beside the
        # note's reason), and a replaced submission names the note as its source (at most one note applies)
        n = next((k for k in ("8", "6") if info[f"note{k}"]["applies"]), None)
        if n:
            info["table_note"], info["note"] = info["note"], f"rules.yaml note {n}: {info[f'note{n}']['reason']}"
            if submission == "s2" or n == "8":
                src = info["source"]
                info["table_source"], info["source"] = src, f"rules.yaml note {n} (replacing {src})"
                info["table_factor_row"], info["factor_row"] = info["factor_row"], None
        if submission == "s2":
            return s2, info
    return tab[submission], info


def validate_spec(spec: dict, rules: dict) -> None:
    src = spec.get("factor", {}).get("source")
    if src != "none" and src not in FACTOR_FILES:
        raise Refused(f"unknown factor source {src!r}")
    if not isinstance(spec["factor"].get("exponent"), int | float):
        raise Refused(f"factor exponent {spec['factor'].get('exponent')!r} is not a number")
    comp = spec.get("composition", {})
    if comp.get("kind") not in COMPOSITION_KINDS:
        raise Refused(f"unknown composition {comp.get('kind')!r}")
    if comp["kind"] != "none":
        scen = rules["presets"]["gata4"]["scenario"]
        if list(comp.get("depleted_states", [])) != list(scen["depleted_states"]):
            raise Refused("the composition's depleted states differ from the rules.yaml scenario")
        l2 = comp.get("log2_magnitude")
        if isinstance(l2, bool) or not isinstance(l2, int | float) or not np.isfinite(l2) or l2 < 0:
            raise Refused(f"composition magnitude {l2!r} is not a non-negative number")
    for k in EXTRAS:
        if k in spec and not isinstance(spec[k], bool):
            raise Refused(f"{k}: expected true or false, got {spec[k]!r}")


def factor_name(source: str, preset: str) -> str | None:
    return None if source == "none" else FACTOR_FILES[source].format(preset=preset)


def run_parameters(dec: dict, run_dir: Path | None) -> tuple[dict, str]:
    """Pooled depth ratios, stage alphas and coverage shares for the extras: from the decision file
    (submissions.parameters, written by the gata4 run) or from <run>/factors/factors.json (`factors`)."""
    stored = dec.get("submissions")
    if isinstance(stored, dict) and stored.get("parameters"):
        return stored["parameters"], "decision file submissions.parameters"
    fj = Path(run_dir) / "factors" / "factors.json" if run_dir else None
    if fj is not None and fj.is_file():
        return json.loads(fj.read_text()).get("parameters") or {}, str(fj)
    return {}, "none"


def load_factor(run_dir: Path | None, name: str, genes: list[str], dec: dict) -> tuple[np.ndarray, dict]:
    """<run>/factors/<name>.csv, checked against the sha256 recorded with it, in the reference's gene order."""
    if run_dir is None:
        raise Refused(f"the submission uses the factor {name}: pass --run-dir (the run directory with factors/)")
    path = Path(run_dir) / "factors" / f"{name}.csv"
    if not path.is_file():
        raise Refused(f"factor file not found: {path}")
    h = sha256(path)
    expected, where = None, None
    stored = dec.get("submissions")
    files = stored.get("parameters", {}).get("factor_files", {}) if isinstance(stored, dict) else {}
    if name in files:
        expected, where = files[name]["sha256"], "decision file submissions.parameters.factor_files"
    else:
        fj = Path(run_dir) / "factors" / "factors.json"
        if fj.is_file():
            rec = json.loads(fj.read_text()).get("files", {}).get(name)
            if rec:
                expected, where = rec["sha256"], str(fj)
    if expected is None:
        raise Refused(f"{path}: no recorded sha256 (neither the decision file nor factors/factors.json has it)")
    if expected != h:
        raise Refused(f"{path}: sha256 {h[:16]} differs from the one recorded in {where} ({expected[:16]})")
    df = pd.read_csv(path, float_precision="round_trip")
    if list(df["gene"].astype(str)) != list(genes):
        raise Refused(f"{path}: genes differ from the WT reference's genes (order included)")
    phi = df["factor"].to_numpy(np.float64)
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused(f"{path}: factor values must be finite and positive")
    return phi, {"name": name, "file": str(path), "sha256": h, "recorded_in": where,
                 "log2_sd": float(np.log2(phi).std()), "n_up": int((phi > 1).sum()), "n_down": int((phi < 1).sum())}


def run_records(run_dir: Path | None, preset: str) -> dict:
    """What the run recorded about its inputs: <run>/<preset>_meta.json (t3_decompose) or factors/factors.json."""
    out = {}
    if run_dir is None:
        return out
    meta = Path(run_dir) / f"{preset}_meta.json"
    if meta.is_file():
        m = json.loads(meta.read_text())
        out["meta"] = {"file": str(meta), "preset": m.get("preset"),
                       "wt_sha256": m.get("files", {}).get("wt_pool", {}).get("sha256"),
                       "stage_sha256": {k: v.get("sha256") for k, v in (m.get("stage_pools") or {}).items()}}
    fj = Path(run_dir) / "factors" / "factors.json"
    if fj.is_file():
        f = json.loads(fj.read_text())
        out["factors_json"] = {"file": str(fj), "preset": f.get("preset"),
                               "wt_sha256": f.get("inputs", {}).get("wt", {}).get("sha256"),
                               "stage_sha256": {k.split(":", 1)[1]: v.get("sha256")
                                                for k, v in f.get("inputs", {}).items() if k.startswith("stage:")}}
    return out


def check_inputs_against_run(records: dict, preset: str, wt_sha: str, stage_sha: dict) -> list[str]:
    """Refuse a WT or stage file that differs from the one the run recorded; returns warnings."""
    warns = []
    if not records:
        return ["no run record (<preset>_meta.json or factors/factors.json) to check the WT file against"]
    for r in records.values():
        if r.get("preset") not in (None, preset):
            raise Refused(f"{r['file']} is a {r['preset']} run, the reference preset is {preset}")
        if r.get("wt_sha256") and r["wt_sha256"] != wt_sha:
            raise Refused(f"the WT file differs from the one the run used ({r['file']})")
        for name, h in stage_sha.items():
            rec_h = (r.get("stage_sha256") or {}).get(name)
            if rec_h and rec_h != h:
                raise Refused(f"stage pool {name}: the file differs from the one the run used ({r['file']})")
    return warns


# ---------------------------------------------------------------------------------------------------------------------
# The WT reference in t3_decompose's form


@dataclass
class Reference:
    """The WT cells a submission resamples (deduplicated), as t3_decompose.build_context prepares them."""

    preset: str
    rules: dict
    genes: list
    V: object                       # t3_decompose.Vocab
    pool: object                    # t3_decompose.Pool: counts, label index, kNN-smoothed proportions
    coords: np.ndarray              # positions of the pool cells
    stages: list = field(default_factory=list)   # t3_decompose.StagePool, only with stage mixing
    record: dict = field(default_factory=dict)

    @property
    def ctx(self) -> SimpleNamespace:
        """The fields of t3_decompose.Context that t3_decompose.build reads."""
        return SimpleNamespace(wt=SimpleNamespace(pool=self.pool), V=self.V, rules=self.rules, stages=self.stages)


@contextmanager
def thread_limits(threads: int):
    """BLAS threads as in t3_decompose's parent process; no OpenMP pool. threadpoolctl limits only the runtimes
    loaded when the context is entered, and sklearn loads the OpenMP one, so it is imported first (as
    t3_decompose.run does through provenance()). Otherwise the kNN smoothing would start an OpenMP pool, and forked
    workers deadlock on it (libgomp is not fork-safe)."""
    import sklearn  # noqa: F401
    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=threads, user_api="blas"), threadpool_limits(limits=1, user_api="openmp"):
        yield


def openmp_limited() -> bool:
    """True when every loaded OpenMP runtime is limited to one thread (safe to fork workers)."""
    from threadpoolctl import threadpool_info

    omp = [i for i in threadpool_info() if i["user_api"] == "openmp"]
    return bool(omp) and all(i["num_threads"] == 1 for i in omp)


def load_reference(wt_path: Path, preset: str, rules: dict, stage_paths: dict | None = None,
                   threads: int = 16) -> Reference:
    """WT pool as in t3_decompose.build_context: load_stage(drop_duplicates=True), the preset's label mode, the
    vocabulary, int counts and kNN-smoothed proportions. With stage_paths (stage mixing on) build_context itself runs
    (no knockout target), which also builds the stage pools."""
    from vec.prep import load_stage

    pre = rules["presets"][preset]
    with thread_limits(threads):
        if stage_paths:
            missing = [n for n in pre["stage_pools"] if n not in stage_paths]
            if missing:
                raise Refused(f"stage mixing is on: pass --stage NAME=PATH for {missing}")
            paths = {"wt": str(wt_path), "targets": [],
                     **{f"stage:{n}": str(stage_paths[n]) for n in pre["stage_pools"]}}
            with tempfile.TemporaryDirectory() as tmp:
                ctx, _ = t3d.build_context(preset, paths, rules, n_jobs=threads, threads=threads, out=Path(tmp))
            wt_ad, V, pool, genes, stages = ctx.harm.ref, ctx.V, ctx.wt.pool, ctx.genes, ctx.stages
        else:
            wt_ad = load_stage(wt_path, drop_duplicates=True)
            t3d.add_label_column(wt_ad, pre["label_mode"])
            genes = list(wt_ad.var_names)
            labels = np.array(sorted(set(wt_ad.obs["t3_label"].astype(str))))
            V = t3d.make_vocab(wt_ad, "t3_label", SimpleNamespace(classes=labels), rules)
            pos = {t: i for i, t in enumerate(V.labels)}
            wlab = np.array([pos[t] for t in wt_ad.obs["t3_label"].astype(str)])
            Cw = t3d.dense_counts(wt_ad)
            sm = rules["components"]["smoothing"]
            pool = t3d.Pool(Cw, wlab, t3d.smoothed_props(Cw, sm["n_pcs"], sm["k"]), V.L)
            stages = []
    if COORD_KEY not in wt_ad.obsm:
        raise Refused(f"{wt_path}: obsm['{COORD_KEY}'] missing (the T3 board needs coordinates)")
    coords = np.asarray(wt_ad.obsm[COORD_KEY], dtype=np.float64)[:, :3]
    return Reference(preset=preset, rules=rules, genes=genes, V=V, pool=pool, coords=coords, stages=stages,
                     record={**t3d.file_record(wt_path, wt_ad), "label_mode": pre["label_mode"],
                             "labels": int(V.L), "stage_pools": [s.name for s in stages]})


# ---------------------------------------------------------------------------------------------------------------------
# Prediction: composition, factor, extras -> t3_decompose.build


def composition(V, spec: dict, params: dict) -> tuple[np.ndarray, dict]:
    """Label shares: the WT composition, the scenario states depleted by 2^-log2_magnitude (t3_decompose.rule_comp)
    and, with coverage on, the knockout's strict non-Mesp1 shares (t3_decompose.coverage_comp; every other label keeps
    the scenario's relative share)."""
    f = V.p.copy()
    info: dict = {"kind": spec["composition"]["kind"], "coverage": bool(spec.get("coverage"))}
    comp = spec["composition"]
    if comp["kind"] != "none":
        f, missing = t3d.rule_comp(V, {"deplete": dict.fromkeys(comp["depleted_states"], 0.0)}, {},
                                   float(comp["log2_magnitude"]))
        found = [s for s in comp["depleted_states"] if s not in missing]
        info.update(log2_magnitude=float(comp["log2_magnitude"]), states_found=found, states_missing=missing,
                    depleted_share_wt=float(V.p[V.index(found)].sum()) if found else 0.0,
                    depleted_share=float(f[V.index(found)].sum()) if found else 0.0)
    if spec.get("coverage"):
        shares = params.get("coverage_shares_mean")
        if not shares:
            raise Refused("coverage is on but no coverage shares were recorded by the run")
        fc = V.p.copy()
        for lab, s in shares.items():
            i = V.index([lab])
            if not i.size:
                raise Refused(f"coverage share for {lab!r}, which is not a label of the WT reference")
            fc[i] = float(s)
        f = t3d.coverage_comp(replace(V, p=f), fc)
        info["nonmesp1_share"] = float(f[V.nm].sum())
    info["shares"] = {t: float(x) for t, x in zip(V.labels, f, strict=True)}
    return f, info


def depth_ratios(V, params: dict) -> np.ndarray:
    d = params.get("depth_ratio_pooled")
    if d is None:
        raise Refused("depth is on but no depth ratios were recorded by the run")
    r = np.ones(V.L)
    for lab, x in d.items():
        i = V.index([lab])
        if not i.size:
            raise Refused(f"depth ratio for {lab!r}, which is not a label of the WT reference")
        r[i] = float(x)
    return r


def stage_alphas(ref: Reference, params: dict) -> np.ndarray:
    a = params.get("stage_alphas_mean")
    names = [s.name for s in ref.stages]
    if a is None or set(a) != set(names):
        raise Refused(f"stage mixing is on but the recorded alphas {a} do not match the stage pools {names}")
    return np.array([float(a[n]) for n in names])


def prediction_inputs(ref: Reference, spec: dict, params: dict, phi: np.ndarray | None) -> tuple[dict, dict]:
    """t3_decompose.build keyword arguments for a submission: composition, factor and exponent, depth, stage."""
    f, info = composition(ref.V, spec, params)
    kw = {"f": f, "phi": phi, "w": float(spec["factor"]["exponent"]) if phi is not None else 1.0,
          "ratio": depth_ratios(ref.V, params) if spec.get("depth") else None,
          "alphas": stage_alphas(ref, params) if spec.get("stage") else None}
    return kw, info


def wt_revert_series(wt: Path, given: dict[str, Path] | None = None) -> dict[str, Path]:
    """The other released wild-type stage files the factor is fitted on: every stage of k7_wt_revert.FULL_PANEL but the
    reference's own (`given` replaces the default paths under data/raw/T2; its names must be exactly those stages)."""
    import k7_wt_revert as k7

    want = sorted(n for n in k7.FULL_PANEL if n != Path(wt).stem)
    if Path(wt).stem not in k7.FULL_PANEL:
        raise Refused(f"--wt-revert: the reference {Path(wt).name} is not one of the stages {list(k7.FULL_PANEL)}")
    series = {n: DATA / "T2" / f"{n}.h5ad" for n in want} if not given else {k: Path(v) for k, v in given.items()}
    if sorted(series) != want:
        raise Refused(f"--wt-revert needs exactly the stages {want}, got {sorted(series)}")
    absent = [str(v) for v in series.values() if not v.is_file()]
    if absent:
        raise Refused(f"--wt-revert: wild-type stage files not found: {absent} (pass --wt-series NAME=PATH)")
    return series


def wt_revert_factors(wt: Path, series: dict[str, Path], genes: list[str],
                      pool_counts: np.ndarray | None = None) -> tuple[np.ndarray, dict, np.ndarray, dict]:
    """The wild-type-only factors of tries K7 and K8 for the reference `wt`, fitted exactly as tried
    (k7_wt_revert.group_means, trend_at and revert_factors; the stage days of vec.harness.STAGE_DAYS; the samples in
    sorted name order, as the tries' roles were): (the pooled factor, one factor per lineage group, each reference
    cell's lineage group in the pool's row order, a record without any cell number). pool_counts: the reference pool's
    counts (`Reference.pool.C`); when given, the reference as loaded here must hold exactly those rows and columns."""
    import k7_wt_revert as k7

    from vec.harness import STAGE_DAYS
    from vec.prep import load_stage

    def day(stem: str) -> float:
        return float(STAGE_DAYS[f"T2/{stem}.h5ad"])

    names = sorted(series)
    try:
        ref = load_stage(wt, genes=genes, drop_duplicates=True)
        others = [k7.group_means(load_stage(series[n], genes=genes, drop_duplicates=True))[0] for n in names]
    except KeyError as e:  # a stage file without one of the reference's genes
        raise Refused(f"--wt-revert: {e}") from e
    if pool_counts is not None and not np.array_equal(t3d.dense_counts(ref), pool_counts):
        raise Refused("--wt-revert: the reference as loaded for the factor is not the pool the cells are drawn from")
    ref_means, ref_sizes = k7.group_means(ref)
    trend = k7.trend_at(others, [day(n) for n in names], day(Path(wt).stem))
    phi, per, rec = k7.revert_factors(ref_means, ref_sizes, trend)
    record = {"rule": k7.WT_REVERT_RULE, "exponent": k7.EXPONENT, "applied": "per lineage group (tries K7b, K8b)",
              "estimator": "experiments/t3_decomp/k7_wt_revert.py",
              "estimator_sha256": sha256(HERE / "k7_wt_revert.py"),
              "series": {n: {"path": str(series[n]), "sha256": sha256(series[n])} for n in names},
              "groups_with_a_factor": list(rec["groups"]), "sd_log_factor_pooled": rec["sd_log_phi"],
              "sd_log_factor_per_group": rec["sd_log_phi_per_group"], "genes_at_clip_pooled": rec["genes_at_clip"]}
    return phi, per, np.asarray(ref.obs["group"]).astype(str), record


def apply_wt_revert(ref: Reference, built: Built, phi: np.ndarray, per: dict, group: np.ndarray,
                    rng: np.random.Generator) -> tuple[Built, dict]:
    """The drawn cells with each one's counts moved by its lineage group's factor at the tried exponent, through
    t3_decompose.apply_gene_factor with the pool's smoothed profile (the arm of tries K7b and K8b); a group without a
    factor of its own takes the pooled one. Draws, labels and source rows are unchanged."""
    import k7_wt_revert as k7

    if not np.all(built.source == "wt"):
        raise Refused("--wt-revert expects every drawn cell to come from the reference pool (no stage mixing)")
    if len(group) != ref.pool.C.shape[0]:
        raise Refused("--wt-revert: the lineage groups do not line up with the reference pool's rows")
    C, P = ref.pool.C[built.row].copy(), ref.pool.P[built.row]
    drawn = group[built.row]
    how = {}
    for g in sorted(set(drawn)):
        rows = np.flatnonzero(drawn == g)
        C[rows] = t3d.apply_gene_factor(C[rows], per.get(g, phi) ** k7.EXPONENT, P[rows], rng)
        how[g] = "its own factor" if g in per else "the pooled factor"
    return Built(X=t3d.lognorm(C), lab=built.lab, source=built.source, row=built.row), {"factor_by_group": how}


@dataclass
class Built:
    """t3_decompose.build output and where each cell came from (its source pool and row)."""

    X: np.ndarray        # log1p CP10k, cells x reference genes, in build's row order
    lab: np.ndarray      # vocabulary index
    source: np.ndarray   # "wt" or "stage:<name>"
    row: np.ndarray      # row in the source pool


@contextmanager
def observe_draws(calls: list):
    """Record every t3_decompose.draw_rows call (the rows build takes, in order) without changing any of them."""
    original = t3d.draw_rows

    def spy(pool_idx, k, rng):
        r = original(pool_idx, k, rng)
        calls.append((id(pool_idx), r))
        return r

    t3d.draw_rows = spy
    try:
        yield
    finally:
        t3d.draw_rows = original


def construct(ref: Reference, n: int, rng: np.random.Generator, f, phi=None, w: float = 1.0, ratio=None,
              alphas=None) -> Built:
    """t3_decompose.build on the reference, exactly as the scored variants call it, plus each cell's source row."""
    calls: list = []
    with observe_draws(calls):
        X = t3d.build(ref.ctx, n, rng, f=f, alphas=alphas, phi=phi, w=w, ratio=ratio)
    ident = {id(ix): ("wt", t) for t, ix in enumerate(ref.pool.idx)}
    for sp in ref.stages:
        ident.update({id(ix): (f"stage:{sp.name}", t) for t, ix in enumerate(sp.pool.idx)})
    src, rows, labs = [], [], []
    for key, r in calls:
        if len(r):
            kind, t = ident[key]
            src.append(np.full(len(r), kind, dtype=object))
            rows.append(np.asarray(r, dtype=np.int64))
            labs.append(np.full(len(r), t))
    lab = np.concatenate(labs)
    if len(lab) != X.shape[0]:
        raise RuntimeError(f"tracked {len(lab)} source rows for {X.shape[0]} cells")
    return Built(X=X, lab=lab, source=np.concatenate(src), row=np.concatenate(rows))


def median_nn(coords: np.ndarray) -> float:
    from sklearn.neighbors import NearestNeighbors

    d, _ = NearestNeighbors(n_neighbors=2).fit(coords).kneighbors(coords)
    return float(np.median(d[:, 1]))


def assemble(ref: Reference, built: Built, panel: list[str], rng: np.random.Generator) -> tuple[np.ndarray,
                                                                                                np.ndarray, dict]:
    """Panel columns, rows in random order, coordinates of each cell's source WT cell (a random WT cell of the same
    label for stage-pool cells, whose own positions are in another embryo's frame), and exact duplicate positions
    moved in xy by at most JITTER_NN x the WT median nearest-neighbour distance."""
    pos = {g: i for i, g in enumerate(ref.genes)}
    missing = [g for g in panel if g not in pos]
    if missing:
        raise Refused(f"{len(missing)} board genes are not in the WT reference, e.g. {missing[:5]}")
    n = built.X.shape[0]
    coords = np.empty((n, 3))
    wt = built.source == "wt"
    coords[wt] = ref.coords[built.row[wt]]
    for t in np.unique(built.lab[~wt]):
        i = np.flatnonzero(~wt & (built.lab == t))
        coords[i] = ref.coords[rng.choice(ref.pool.idx[t], len(i), replace=True)]
    perm = rng.permutation(n)
    X = built.X[perm][:, [pos[g] for g in panel]]
    coords = coords[perm]
    _, first = np.unique(coords, axis=0, return_index=True)
    dup = np.ones(n, dtype=bool)
    dup[first] = False
    radius = JITTER_NN * median_nn(ref.coords) if dup.any() else 0.0
    if dup.any():
        r = radius * np.sqrt(rng.random(int(dup.sum())))
        th = rng.uniform(0.0, 2.0 * np.pi, int(dup.sum()))
        coords[dup, 0] += r * np.cos(th)
        coords[dup, 1] += r * np.sin(th)
    sources = pd.Series(built.source).value_counts().to_dict()
    return X, coords.astype(np.float32), {"row_order": "random permutation", "cells_by_source": sources,
                                         "coords": "source WT cell (stage-pool cells: a random WT cell of the label)",
                                         "duplicate_positions_moved": int(dup.sum()), "jitter_radius": radius}


# ---------------------------------------------------------------------------------------------------------------------
# build


def board_spec(key: str, boards_dir: Path | None = None):
    """A board's submission spec from boards/ (vec.boards); its anchors are dropped unread."""
    from vec.boards import BOARDS_DIR, get_board

    try:
        b = get_board(key, boards_dir or BOARDS_DIR)
    except KeyError as e:
        raise Refused(str(e)) from e
    return replace(b, anchors={})


def preflight(decision: Path, branch: str, submission: str, board, n_cells: int, preset: str | None = None,
              allow_hash_mismatch: bool = False, allow_unregistered_protocol: bool = False) -> dict:
    """Every check that needs no data: the branch and submission, the decision file, its hashes and protocol, the
    submission it names, the reference preset and the board."""
    rules = load_rules()
    if branch not in BRANCHES:
        raise Refused(f"unknown branch {branch!r}; expected one of {list(BRANCHES)}")
    if submission not in SUBMISSIONS:
        raise Refused(f"unknown submission {submission!r}; expected one of {list(SUBMISSIONS)}")
    dec = load_decision(decision, rules)
    hashes = check_hashes(dec, allow_hash_mismatch)
    warnings: list[str] = []
    if not dec.get("protocol_as_registered", False):
        if not allow_unregistered_protocol:
            raise Refused("the decision file was not made with the registered protocol (seeds); pass "
                          "--allow-unregistered-protocol to build anyway (recorded)")
        warnings.append("decision file made with a non-registered protocol")
    preset = preset or dec["preset"]
    if preset not in rules["presets"]:
        raise Refused(f"unknown preset {preset!r}")
    if preset != dec["preset"] and branch != "barred":
        raise Refused(f"--preset {preset} differs from the decision file's {dec['preset']}: only the barred branch "
                      "(which reads nothing but D7m) may build on another preset's reference")
    try:
        t3d.check_scenario(rules)
    except SystemExit as e:
        raise Refused(str(e)) from e
    spec, spec_info = resolve_spec(dec, rules, branch, submission, hashes)
    validate_spec(spec, rules)
    if spec_info.get("warning"):
        warnings.append(f"decision file: {spec_info['warning']}")
    if board.task != "T3":
        raise Refused(f"board {board.key} is a {board.task} board")
    if n_cells < board.min_cells or (board.max_cells is not None and n_cells > board.max_cells):
        raise Refused("--n-cells is outside the board's limits (values not shown: test board)" if board.counts_hidden
                      else f"--n-cells {n_cells} is outside the board's limits [{board.min_cells}, {board.max_cells}]")
    if n_cells < 2000:
        warnings.append(("fewer than the 2,000 cells the MMD draws (count not shown: test board)" if board.counts_hidden
                         else f"{n_cells} cells: fewer than the 2,000 the MMD draws") + " (docs/PLAN.md section 3.3)")
    return {"rules": rules, "dec": dec, "hashes": hashes, "spec": spec, "spec_info": spec_info, "preset": preset,
            "warnings": warnings}


def resolve_n_cells(n_cells: str | int, board) -> int:
    """--n-cells: "max" (the default) is the board's max_cells, the approved T3 count, read in code and never printed
    or typed (a test board's limits can encode the held-out cell count); an integer is for tests and dry runs."""
    if str(n_cells).strip().lower() in ("max", "max_cells"):
        if board.max_cells is None:
            raise Refused(f"board {board.key} publishes no max_cells; the approved T3 count is max_cells")
        return int(board.max_cells)
    try:
        return int(n_cells)
    except ValueError:
        raise Refused("--n-cells takes 'max' or an integer") from None


def count_label(n: int, board):
    """A cell count as it may be shown: the number on validation and synthetic boards; on a test board "max_cells"
    when it equals the board's max_cells, else "not max_cells" (values never shown)."""
    if not board.counts_hidden:
        return int(n)
    return "max_cells" if board.max_cells is not None and int(n) == int(board.max_cells) else "not max_cells"


def inputs_for(spec: dict, dec: dict, ref: Reference, run_dir: Path | None) -> tuple[dict, dict]:
    """What t3_decompose.build needs for this submission: the factor file it names (sha256-checked), the run's
    parameters for the extras, then prediction_inputs."""
    fname = factor_name(spec["factor"]["source"], dec["preset"])
    needs_params = any(spec.get(k) for k in EXTRAS)
    params, params_src = run_parameters(dec, run_dir) if needs_params else ({}, "not needed")
    phi, frec = load_factor(run_dir, fname, ref.genes, dec) if fname else (None, None)
    kw, comp_info = prediction_inputs(ref, spec, params, phi)
    return kw, {"factor": frec, "parameters_source": params_src, "composition": comp_info}


def build_submission(decision: Path, branch: str, submission: str, wt: Path, board, n_cells: int, out: Path,
                     run_dir: Path | None = None, preset: str | None = None, stage_paths_in: dict | None = None,
                     seed: int = 1_000_001, allow_hash_mismatch: bool = False,
                     allow_unregistered_protocol: bool = False, overwrite: bool = False, threads: int = 16,
                     wt_revert: dict[str, Path] | None = None) -> dict:
    """Build one submission from the decision file; returns the sidecar (also written next to the file). The seed is
    never below vec.harness.BUILD_SEED_OFFSET (docs/P3_RUNBOOK.md section 5: numpy's default_rng([s, 0]) is
    default_rng(s), and the server draws its working copies with seed 0)."""
    from vec.harness import BUILD_SEED_OFFSET
    from vec.io import check_submission, make_submission, write_submission

    if int(seed) < BUILD_SEED_OFFSET:
        raise Refused(f"build seed {seed} is below {BUILD_SEED_OFFSET} (docs/P3_RUNBOOK.md section 5: never build a "
                      "test file at seed 0)")
    if wt_revert is not None and (branch, submission) != WT_REVERT_SLOT:
        raise Refused("--wt-revert is ruled for T3's second upload only (--branch refused --submission s1; "
                      "docs/S2_RULES.md, ruling of 2026-11-13)")
    t0 = time.time()
    pf = preflight(decision, branch, submission, board, n_cells, preset, allow_hash_mismatch,
                   allow_unregistered_protocol)
    rules, dec, spec, preset, warnings = pf["rules"], pf["dec"], pf["spec"], pf["preset"], pf["warnings"]
    if wt_revert is not None:  # the form tries K8 read: the refused s1 as the decision file states it
        import k8_s2_form as k8

        tried, comp = k8.spec()["composition"], spec.get("composition", {})
        if spec.get("coverage") or any(comp.get(k) != tried[k] for k in ("kind", "depleted_states", "log2_magnitude")):
            raise Refused("--wt-revert was tried on the refused branch's s1 with the central composition of the "
                          "decision file (tries K8); this submission's composition differs")
    out = Path(out)
    side = out.with_suffix(".json")
    if not overwrite and (out.exists() or side.exists()):
        raise Refused(f"{out} or its sidecar exists; pass --overwrite")
    stage_paths = {}
    if spec.get("stage"):
        pools = rules["presets"][preset]["stage_pools"]
        stage_paths = {k: Path(v) for k, v in (stage_paths_in or {}).items()} or {
            n: DATA / "T2" / f"{n}.h5ad" for n in pools}
        absent = [str(p) for p in stage_paths.values() if not p.is_file()]
        if absent:
            raise Refused(f"stage mixing is on: stage pool files not found: {absent} (pass --stage NAME=PATH)")

    ref = load_reference(Path(wt), preset, rules, stage_paths or None, threads)
    stage_sha = {n: sha256(p) for n, p in stage_paths.items()}
    records = run_records(run_dir, preset)
    warnings += check_inputs_against_run(records, preset, ref.record["sha256"], stage_sha)
    kw, info = inputs_for(spec, dec, ref, run_dir)
    if info["composition"].get("states_missing"):
        warnings.append("scenario states not in the WT vocabulary (left unchanged, as in the scored scenario "
                        f"variant): {info['composition']['states_missing']}")
    stream = f"submission:{branch}:{submission}"
    built = construct(ref, n_cells, t3d.stable_rng(seed, stream), **kw)
    revert, streams = None, [stream, f"{stream}:assemble"]
    if wt_revert is not None:
        if kw["phi"] is not None or kw["ratio"] is not None or kw["alphas"] is not None:
            raise Refused("--wt-revert is for a submission with no other factor, depth ratio or stage mixing")
        phi_wr, per_wr, group_wr, revert = wt_revert_factors(Path(wt), wt_revert, ref.genes, ref.pool.C)
        built, applied = apply_wt_revert(ref, built, phi_wr, per_wr, group_wr,
                                         t3d.stable_rng(seed, f"{stream}:{WT_REVERT_STREAM}"))
        revert = {**revert, **applied, "ruling": "the captain, 2026-11-13 (docs/S2_RULES.md)"}
        streams.append(f"{stream}:{WT_REVERT_STREAM}")
    X, coords, ainfo = assemble(ref, built, board.genes, t3d.stable_rng(seed, f"{stream}:assemble"))
    write_submission(make_submission(X, board.genes, coords), out)
    rep = check_submission(out, board)
    sidecar = {
        "builder": "experiments/t3_decomp/build_submission.py",
        "created": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "argv": sys.argv,
        "decision_file": {"path": str(decision), "sha256": sha256(decision), "preset": dec["preset"],
                          "protocol_as_registered": dec.get("protocol_as_registered"), "seeds": dec.get("seeds")},
        "branch": branch,
        "submission": submission,
        "spec": spec,
        "spec_source": pf["spec_info"],
        "hashes": {**pf["hashes"], "build_submission_sha256": sha256(BUILDER_PATH)},
        "inputs": {"wt": ref.record, "reference_preset": preset, "run_dir": str(run_dir) if run_dir else None,
                   "run_records": records, "factor": info["factor"], "parameters_source": info["parameters_source"],
                   "stage_files": {n: {"path": str(p), "sha256": stage_sha[n]} for n, p in stage_paths.items()},
                   "board": {"key": board.key, "task": board.task, "n_genes": board.n_genes,
                             "genes_sha256": hashlib.sha256("\n".join(board.genes).encode()).hexdigest(),
                             **({"limits": "hidden (test board)"} if board.counts_hidden else
                                {"min_cells": board.min_cells, "max_cells": board.max_cells})},
                   "n_cells": count_label(n_cells, board), "seed": int(seed),
                   "rng_streams": streams},
        "construction": {"function": "t3_decompose.build", "wt_revert": revert,
                         "factor_exponent": kw["w"] if kw["phi"] is not None else None,
                         "depth_ratios": None if kw["ratio"] is None else
                         {t: float(r) for t, r in zip(ref.V.labels, kw["ratio"], strict=True) if r != 1.0},
                         "stage_alphas": None if kw["alphas"] is None else
                         dict(zip([s.name for s in ref.stages], map(float, kw["alphas"]), strict=True)),
                         "composition": info["composition"],
                         **{k: (v if k != "cells_by_source" or not board.counts_hidden else
                                {src: round(c / max(1, int(X.shape[0])), 4) for src, c in v.items()})
                            for k, v in ainfo.items()}},
        "output": {"h5ad": str(out), "sha256": sha256(out), "n_cells": count_label(X.shape[0], board),
                   "check": {"ok": rep.ok, "errors": rep.errors, "warnings": rep.warnings}},
        "warnings": warnings,
        "provenance": t3d.provenance({"wt": str(wt), "decision": str(decision)}),
        "seconds": round(time.time() - t0, 1),
    }
    dump_json(sidecar, side)
    return sidecar


def parse_pairs(items: list[str] | None, what: str) -> dict[str, Path]:
    out = {}
    for s in items or []:
        if "=" not in s:
            raise Refused(f"{what} {s!r}: expected NAME=PATH")
        k, v = s.split("=", 1)
        out[k] = Path(v)
    return out


@app.command()
def build(
    decision: Path = typer.Argument(..., help="decision file (gata4_decisions.json; the Mab21l2 one for testing)"),
    branch: str = typer.Option(..., help="allowed | refused | barred (rules.yaml note 3: Gata4 files not usable)"),
    submission: str = typer.Option(..., help="s1 | s2"),
    wt: Path = typer.Option(..., help="WT reference (test: WT E8.75), loaded with drop_duplicates=True"),
    board: str = typer.Option(..., help="T3 board key in boards/ (gene panel, cell limits)"),
    n_cells: str = typer.Option("max", help="'max' (default): the board's max_cells, read in code and never "
                                "printed (the approved count); an integer only for tests and dry runs"),
    out: Path = typer.Option(..., help="output .h5ad (the sidecar is the same name with .json)"),
    run_dir: Path = typer.Option(None, help="run directory holding factors/ (DGX only)"),
    preset: str = typer.Option(None, help="reference preset (default: the decision file's; barred branch only)"),
    stage: list[str] = typer.Option(None, help="NAME=PATH stage pool, needed only when stage mixing is on"),
    seed: int = typer.Option(1_000_001, help="seed of the resampling and row order; never below "
                             "vec.harness.BUILD_SEED_OFFSET (docs/P3_RUNBOOK.md section 5)"),
    boards_dir: Path = typer.Option(None, help="boards directory (default: the repository's boards/)"),
    allow_hash_mismatch: bool = typer.Option(False, help="build although the decision file's script/rules hashes "
                                             "differ from the files on disk (recorded)"),
    allow_unregistered_protocol: bool = typer.Option(False, help="accept a decision file from a non-registered "
                                                     "protocol (recorded)"),
    overwrite: bool = typer.Option(False, help="replace an existing output"),
    threads: int = typer.Option(16, help="BLAS threads while preparing the reference"),
    wt_revert: bool = typer.Option(False, "--wt-revert", help="T3's second upload only (--branch refused --submission "
                                   "s1; the captain's ruling of 2026-11-13): also apply the wild-type-only factor of "
                                   "tries K7b and K8b per lineage group. Off: the file as the frozen rules built it"),
    wt_series: list[str] = typer.Option(None, help="NAME=PATH wild-type stage files for --wt-revert (default: the "
                                        "other full-panel stages under data/raw/T2)"),
):
    """Build submission s1 or s2 of one branch from the decision file."""
    try:
        b = board_spec(board, boards_dir)
        series = wt_revert_series(wt, parse_pairs(wt_series, "--wt-series")) if wt_revert else None
        if wt_series and not wt_revert:
            raise Refused("--wt-series is given without --wt-revert")
        s = build_submission(decision, branch, submission, wt, b, resolve_n_cells(n_cells, b), out, run_dir=run_dir,
                             preset=preset,
                             stage_paths_in=parse_pairs(stage, "--stage"), seed=seed,
                             allow_hash_mismatch=allow_hash_mismatch,
                             allow_unregistered_protocol=allow_unregistered_protocol, overwrite=overwrite,
                             threads=threads, wt_revert=series)
    except Refused as e:
        say(f"refused: {e}", "red")
        raise typer.Exit(2) from None
    sp, c = s["spec"], s["construction"]
    comp = sp["composition"]
    say(f"{branch} {submission} from {s['spec_source']['source']}")
    if (s["spec_source"].get("note8") or {}).get("applies"):
        say(f"  rules.yaml note 8 (the factor table's last row replaced): {s['spec_source']['note8']['reason']}",
            "yellow")
    if s["spec_source"].get("note6"):
        n6 = s["spec_source"]["note6"]
        say(f"  rules.yaml note 6 ({'s2 replaced' if n6['applies'] else 's2 as in the table'}; s1 never changes): "
            f"{n6['reason']}; D11b per target {n6['D11b_per_target']}",
            "yellow" if n6["applies"] and submission == "s2" else None)
    say(f"  factor {sp['factor']['source']} exponent {sp['factor']['exponent']}"
        + (f" ({s['inputs']['factor']['file']})" if s["inputs"]["factor"] else ""))
    if c.get("wt_revert"):
        say(f"  plus the wild-type-only factor per lineage group at exponent {c['wt_revert']['exponent']} "
            "(--wt-revert; the captain's ruling of 2026-11-13)")
    say(f"  composition {comp['kind']}" + (f", log2 magnitude {comp['log2_magnitude']:.3f}"
                                           if comp["kind"] != "none" else "")
        + "".join(f", {k} on" for k in EXTRAS if sp.get(k)))
    if b.counts_hidden:
        say(f"  cells: {s['output']['n_cells']} (counts not shown: test board); sources "
            f"{sorted(c['cells_by_source'])}")
    else:
        say(f"  {s['output']['n_cells']} cells ({c['cells_by_source']}); board max_cells {b.max_cells}")
    for w in s["warnings"]:
        say(f"  warning: {w}", "yellow")
    if s["hashes"]["mismatch"]:
        say(f"  hash mismatch accepted and recorded: {s['hashes']['mismatch']}", "yellow")
    chk = s["output"]["check"]
    for w in chk["warnings"]:
        say(f"  check warning: {w}", "yellow")
    for e in chk["errors"]:
        say(f"  check error: {e}", "red")
    say(f"{'OK' if chk['ok'] else 'FAILED the format check'}: {s['output']['h5ad']} (sidecar "
        f"{Path(s['output']['h5ad']).with_suffix('.json')})", "green" if chk["ok"] else "red")
    raise typer.Exit(0 if chk["ok"] else 1)


# ---------------------------------------------------------------------------------------------------------------------
# factors: per-gene factors by t3_decompose's own code, without any Gata4 file


def preset_paths(preset: str, wt=None, ko=None, mab21l2=None, e875=None, e825=None, e95=None, e85=None) -> dict:
    """Input paths as t3_decompose.main builds them; the gata4 preset gets no knockout replicate (no Gata4 file)."""
    T2, T3 = DATA / "T2", DATA / "T3"
    if preset == "mab21l2":
        ko = ko or T3 / "E9.5_mab21l2_ko.h5ad"
        e875 = e875 or T2 / "E8.75.h5ad"
        e825 = e825 or T2 / "E8.25_late.h5ad"
        return {"wt": str(wt or T2 / "E9.5.h5ad"), "targets": [("mab21l2", str(ko))], "stage:E8.75": str(e875),
                "series:E8.25_late": str(e825), "series:E8.75": str(e875), "transfer_ko": str(ko)}
    if preset == "gata4":
        e95 = e95 or T2 / "E9.5.h5ad"
        e825 = e825 or T2 / "E8.25_late.h5ad"
        e85 = e85 or (T2 / "E8.5.h5ad" if (T2 / "E8.5.h5ad").is_file() else None)
        return {"wt": str(wt or T2 / "E8.75.h5ad"), "targets": [], "stage:E9.5": str(e95),
                "stage:E8.25_late": str(e825), "series:E8.25_late": str(e825), "series:E9.5": str(e95),
                "series:E8.5": None if e85 is None else str(e85),
                "transfer_ko": str(mab21l2 or T3 / "E9.5_mab21l2_ko.h5ad")}
    raise Refused(f"unknown preset {preset!r}")


def write_factors(preset: str, paths: dict, out: Path, threads: int = 16) -> dict:
    """Factor files as t3_decompose.submission_parameters writes them (gene, factor), from released files:
    reference_R, the Mab21l2 transfer factors and, for the mab21l2 preset, the knockout's pooled non-Mesp1 factor
    (mab21l2_pooled) with its depth ratios, stage alphas and coverage shares. factors/factors.json records sha256s."""
    rules = load_rules()
    need = [p for k, p in paths.items() if k != "targets" and p is not None and not Path(p).is_file()]
    need += [p for _, p in paths["targets"] if not Path(p).is_file()]
    if need:
        raise Refused(f"released files not found: {need}")
    out.mkdir(parents=True, exist_ok=True)
    with thread_limits(threads):
        ctx, meta = t3d.build_context(preset, paths, rules, n_jobs=threads, threads=threads, out=out)
        tmeta: dict = {}
        transfer = t3d.transfer_factors(ctx, paths["transfer_ko"], tmeta)
    V = ctx.V
    factors = {"reference_R": ctx.phi_R, **{f"transfer_{k}": v for k, v in transfer.items()}}
    params: dict = {}
    if preset == "mab21l2":
        E = ctx.targets[0].full
        factors[f"{preset}_pooled"] = E.phi_nm
        params = {"depth_ratio_pooled": {V.labels[i]: float(E.ratio[i]) for i in range(V.L) if E.ratio[i] != 1.0},
                  "stage_alphas_mean": dict(zip([s.name for s in ctx.stages], map(float, E.alphas), strict=True)),
                  "coverage_shares_mean": {V.labels[i]: float(E.f["top"]["cov"][i]) for i in np.flatnonzero(V.nm)}}
    fdir = out / "factors"
    fdir.mkdir(exist_ok=True)
    files = {}
    for k, v in factors.items():
        if v is None:
            continue
        p = fdir / f"{k}.csv"
        pd.DataFrame({"gene": ctx.genes, "factor": v}).to_csv(p, index=False)
        files[k] = {"file": str(p), "sha256": sha256(p)}
    inputs = {k: {"path": p, "sha256": sha256(p)} for k, p in paths.items()
              if k not in ("targets",) and p is not None}
    inputs.update({f"target:{n}": {"path": p, "sha256": sha256(p)} for n, p in paths["targets"]})
    rec = {"preset": preset, "created": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "files": files, "parameters": params,
           "inputs": inputs, "wt_pool": meta["files"]["wt_pool"], "R": meta.get("R"),
           "transfer_ko": tmeta.get("transfer_ko"),
           "hashes": {"t3_decompose": sha256(T3D_PATH), "rules": sha256(RULES_PATH),
                      "build_submission": sha256(BUILDER_PATH)},
           "git": t3d.git_state(),
           "note": "written by build_submission.py factors with t3_decompose.build_context, estimate (full pool) and "
                   "transfer_factors; no Gata4 file is read"}
    dump_json(rec, fdir / "factors.json")
    return rec


@app.command()
def plan(
    decision: Path = typer.Argument(..., help="decision file (gata4_decisions.json)"),
    out: Path = typer.Option(None, help="also write the plan as JSON"),
    allow_hash_mismatch: bool = typer.Option(False, help="show the plan although the decision file's script/rules "
                                             "hashes differ from the files on disk (recorded)"),
):
    """Both submissions of every branch as `build` would make them (the decision tables with rules.yaml notes 3, 6, 7
    and 8 applied). Reads only the decision file and builds nothing."""
    try:
        rules = load_rules()
        dec = load_decision(decision, rules)
        hashes = check_hashes(dec, allow_hash_mismatch)
        registered = bool(dec.get("protocol_as_registered", False))
        branches: dict = {}
        for branch in BRANCHES:
            branches[branch] = {}
            for sub in SUBMISSIONS:
                spec, info = resolve_spec(dec, rules, branch, sub, hashes)
                validate_spec(spec, rules)
                branches[branch][sub] = {"spec": spec, "source": info}
    except Refused as e:
        say(f"refused: {e}", "red")
        raise typer.Exit(2) from None
    for branch, subs in branches.items():
        for sub, r in subs.items():
            sp, comp = r["spec"], r["spec"]["composition"]
            say(f"{branch} {sub}: factor {sp['factor']['source']} exponent {sp['factor']['exponent']}, composition "
                f"{comp['kind']}" + (f" (log2 {comp['log2_magnitude']:.3f})" if comp["kind"] != "none" else "")
                + "".join(f", {k} on" for k in EXTRAS if sp.get(k)))
            if r["source"].get("note"):
                say(f"  {r['source']['note']}")
            if r["source"].get("warning"):
                say(f"  warning: {r['source']['warning']}", "yellow")
    if not registered:
        say("warning: the decision file was not made with the registered protocol; build refuses it without "
            "--allow-unregistered-protocol", "yellow")
    if out:
        # the builder's hash too: notes 6-8 live in this file, which the decision file's hashes do not cover
        dump_json({"command": "plan", "created": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
                   "decision_file": {"path": str(decision), "sha256": sha256(decision),
                                     "protocol_as_registered": registered},
                   "hashes": {**hashes, "build_submission": sha256(BUILDER_PATH)}, "branches": branches}, out)
        say(f"written {out}", "green")


@app.command()
def factors(
    preset: str = typer.Option(..., help="mab21l2 (smoke test, vs WT E9.5) | gata4 (WT E8.75, no Gata4 file)"),
    out: Path = typer.Option(..., help="output directory (factors/ goes inside)"),
    wt: Path = typer.Option(None, help="WT reference (default: E9.5 for mab21l2, E8.75 for gata4)"),
    ko: Path = typer.Option(None, help="mab21l2: the Mab21l2 KO (default data/raw/T3/E9.5_mab21l2_ko.h5ad)"),
    mab21l2: Path = typer.Option(None, help="gata4: the Mab21l2 KO for the transfer factor"),
    e875: Path = typer.Option(None, help="mab21l2: WT E8.75 (stage pool and R series)"),
    e825: Path = typer.Option(None, help="WT E8.25_late"),
    e95: Path = typer.Option(None, help="gata4: WT E9.5 (stage pool and R series)"),
    e85: Path = typer.Option(None, help="gata4: WT E8.5 for R (default data/raw/T2/E8.5.h5ad if present)"),
    overwrite: bool = typer.Option(False, help="write into an existing factors/ directory"),
    threads: int = typer.Option(16, help="BLAS threads"),
):
    """Write the per-gene factors a decision names, from released files by t3_decompose's own functions."""
    try:
        if (out / "factors").exists() and not overwrite:
            raise Refused(f"{out / 'factors'} exists; pass --overwrite")
        if preset == "gata4" and ko is not None:
            raise Refused("--preset gata4 reads no knockout replicate (no Gata4 file); --ko is for mab21l2")
        paths = preset_paths(preset, wt, ko, mab21l2, e875, e825, e95, e85)
        rec = write_factors(preset, paths, out, threads)
    except Refused as e:
        say(f"refused: {e}", "red")
        raise typer.Exit(2) from None
    for k, v in rec["files"].items():
        say(f"{k}: {v['file']} ({v['sha256'][:16]})")
    say(f"recorded in {out / 'factors' / 'factors.json'}", "green")


# ---------------------------------------------------------------------------------------------------------------------
# replay: scored Mab21l2 variants rebuilt through this builder


_REPLAY: dict = {}


def same_reference(a: Reference, ctx) -> dict:
    """Is the builder's reference (light path) the run's WT pool? Exact comparisons, plus the largest difference of
    the smoothed proportions."""
    p = ctx.wt.pool
    return {"labels": a.V.labels == ctx.V.labels, "shares": bool(np.array_equal(a.V.p, ctx.V.p)),
            "counts": bool(np.array_equal(a.pool.C, p.C)), "label_index": bool(np.array_equal(a.pool.lab, p.lab)),
            "smoothed_props_max_abs_diff": float(np.abs(a.pool.P.astype(np.float64) - p.P).max()),
            "nonmesp1": bool(np.array_equal(a.V.nm, ctx.V.nm))}


def replay_seed(seed: int) -> list[dict]:
    """One seed as t3_decompose.score_seed builds it: the seed's working copy, the knockout pool without it, the
    per-seed estimate; then each variant through prediction_inputs and construct, scored by vec.evaluate."""
    from vec.evaluate import Stage, score_board

    st = _REPLAY
    ctx, ref, ref_stage, pr = st["ctx"], st["ref"], st["ref_stage"], st["rules"]["protocol"]
    T = ctx.targets[0]
    frac = pr["working_copy_fraction"]
    rows_t = t3d.working_copy_rows(len(T.X_rel), frac, seed)
    n_ref = len(t3d.working_copy_rows(len(ctx.ref_X), frac, seed + 1))
    n_pred = pr["prediction_cells_multiplier"] * max(len(rows_t) // 2, n_ref)
    pool_rows = np.flatnonzero(~np.isin(T.hash_pool, T.hash_rel[rows_t]))
    est = t3d.estimate(ctx, T, pool_rows, t3d.stable_rng(seed, f"estimate:{T.name}"))
    V = ctx.V
    params = {"depth_ratio_pooled": {V.labels[i]: float(est.ratio[i]) for i in range(V.L) if est.ratio[i] != 1.0},
              "stage_alphas_mean": dict(zip([s.name for s in ctx.stages], map(float, est.alphas), strict=True)),
              "coverage_shares_mean": {V.labels[i]: float(est.f["top"]["cov"][i]) for i in np.flatnonzero(V.nm)}}
    out = []
    for name in st["variants"]:
        on = REPLAY_VARIANTS[name]
        spec = {"factor": {"source": "replay" if on.get("factor") else "none", "exponent": 1.0},
                "composition": {"kind": "none"}, **{k: bool(on.get(k)) for k in EXTRAS}}
        r = ref_stage if spec["stage"] else ref
        kw, _ = prediction_inputs(r, spec, params, est.phi_nm if on.get("factor") else None)
        built = construct(r, n_pred, t3d.stable_rng(seed, name), **kw)
        sc = score_board(st["board"], Stage(X=built.X), Stage(X=T.X_rel), Stage(X=ctx.ref_X), seed=seed)
        out.append({"seed": seed, "variant": name, "n_pred": int(n_pred), "t3_replay": sc.score,
                    **{f"replay_{m}": sc.metrics.get(m) for m in ("de_score", "de_direction", "severity_slope",
                                                                  "mmd_u", "variogram")}})
    print(f"[replay seed {seed}] {len(out)} variants, n_pred {n_pred}", flush=True)
    return out


def _replay_task(seed: int) -> list[dict]:
    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=_REPLAY["threads"]):
        return replay_seed(seed)


@app.command()
def replay(
    out: Path = typer.Option(..., help="output JSON (per-seed table next to it as .csv)"),
    records: Path = typer.Option(HERE / "results" / "mab21l2" / "mab21l2_records.csv", help="the run's records"),
    variants: str = typer.Option(",".join(REPLAY_VARIANTS), help="comma-separated scored variants to rebuild"),
    seeds: str = typer.Option(None, help="seeds (default: rules.yaml protocol seeds)"),
    wt: Path = typer.Option(None, help="WT E9.5"),
    ko: Path = typer.Option(None, help="Mab21l2 KO"),
    e875: Path = typer.Option(None, help="WT E8.75"),
    e825: Path = typer.Option(None, help="WT E8.25_late"),
    workers: int = typer.Option(8, help="worker processes (the recorded run used 8 x 2 threads)"),
    threads: int = typer.Option(2, help="BLAS threads per worker"),
):
    """Rebuild scored Mab21l2 variants with this builder's construction and compare with the recorded scores."""
    import multiprocessing as mp

    from vec.boards import Board

    rules = load_rules()
    names = [v.strip() for v in variants.split(",") if v.strip()]
    unknown = [v for v in names if v not in REPLAY_VARIANTS]
    if unknown:
        say(f"refused: cannot replay {unknown}; known: {list(REPLAY_VARIANTS)}", "red")
        raise typer.Exit(2)
    seed_list = t3d.parse_seeds(seeds) if seeds else list(rules["protocol"]["seeds"])
    paths = preset_paths("mab21l2", wt, ko, None, e875, e825)
    t0 = time.time()
    n_par = workers * threads
    rows: list = []
    # as t3_decompose.run: the parent never starts an OpenMP pool before forking (libgomp is not fork-safe)
    with thread_limits(n_par), tempfile.TemporaryDirectory() as tmp:
        ctx, _ = t3d.build_context("mab21l2", paths, rules, n_jobs=n_par, threads=threads, out=Path(tmp))
        ref = load_reference(Path(paths["wt"]), "mab21l2", rules, threads=n_par)
        same = same_reference(ref, ctx)
        say(f"builder reference vs the run's WT pool: {same}")
        ref_stage = Reference(preset="mab21l2", rules=rules, genes=ctx.genes, V=ctx.V, pool=ctx.wt.pool,
                              coords=ref.coords, stages=ctx.stages)
        _REPLAY.update(ctx=ctx, ref=ref, ref_stage=ref_stage, rules=rules, variants=names, threads=threads,
                       board=Board(key="T3:emulated", task="T3", genes=ctx.genes, needs_coords=False))
        if workers > 1:
            if not openmp_limited():
                raise RuntimeError("the OpenMP runtime is not limited to one thread: forked workers would deadlock")
            with mp.get_context("fork").Pool(min(workers, len(seed_list))) as pool:
                for r in pool.imap_unordered(_replay_task, seed_list):
                    rows.extend(r)
        else:
            for s in seed_list:
                rows.extend(_replay_task(s))
    df = pd.DataFrame(rows)
    rec = pd.read_csv(records)
    rec = rec[rec["target"] == ctx.targets[0].name][["seed", "variant", "t3", "n_pred"]]
    df = df.merge(rec.rename(columns={"t3": "t3_recorded", "n_pred": "n_pred_recorded"}), on=["seed", "variant"],
                  how="left").sort_values(["variant", "seed"])
    df["diff"] = df["t3_replay"] - df["t3_recorded"]
    summary = {v: {"seeds": int(len(g)), "t3_recorded_mean": float(g["t3_recorded"].mean()),
                   "t3_replay_mean": float(g["t3_replay"].mean()), "max_abs_diff": float(g["diff"].abs().max()),
                   "exact_1e-9": bool((g["diff"].abs() <= 1e-9).all()), "n_pred_same": bool(
                       (g["n_pred"] == g["n_pred_recorded"]).all())} for v, g in df.groupby("variant")}
    res = {"command": "replay", "records": str(records), "records_sha256": sha256(records), "seeds": seed_list,
           "reference_same_as_run": same, "summary": summary, "hashes": {"t3_decompose": sha256(T3D_PATH),
                                                                          "rules": sha256(RULES_PATH),
                                                                          "build_submission": sha256(BUILDER_PATH)},
           "inputs": {k: v for k, v in paths.items() if k != "targets"}, "minutes": round((time.time() - t0) / 60, 1)}
    dump_json(res, out)
    df.to_csv(out.with_suffix(".csv"), index=False)
    t = Table("variant", "seeds", "recorded mean", "replay mean", "max |diff|", "exact")
    for v, s in summary.items():
        t.add_row(escape(v), str(s["seeds"]), f"{s['t3_recorded_mean']:.4f}", f"{s['t3_replay_mean']:.4f}",
                  f"{s['max_abs_diff']:.2e}", str(s["exact_1e-9"]))
    console.print(t)
    say(f"written {out} and {out.with_suffix('.csv')}", "green")


# ---------------------------------------------------------------------------------------------------------------------
# score: a built file under vec.evaluate, optionally against a recorded variant


def compare_with_records(scores: dict[int, float], records: Path, variant: str, target: str | None) -> dict:
    """Paired per-seed comparison with a recorded variant; within seed noise = |mean difference| <= the recorded
    variant's SD over seeds."""
    rec = pd.read_csv(records)
    rec = rec[rec["variant"] == variant]
    if target is None:
        tg = rec["target"].unique()
        if len(tg) != 1:
            raise Refused(f"records hold targets {list(tg)}: pass --record-target")
        target = str(tg[0])
    s = rec[rec["target"] == target].set_index("seed")["t3"]
    seeds = [k for k in scores if k in s.index]
    if not seeds:
        raise Refused(f"no recorded seeds of {variant} for {target}")
    ours = np.array([scores[k] for k in seeds])
    theirs = s.loc[seeds].to_numpy(np.float64)
    d = ours - theirs
    sd_d = float(d.std(ddof=1)) if len(d) > 1 else 0.0
    rec_sd = float(theirs.std(ddof=1)) if len(theirs) > 1 else 0.0
    return {"variant": variant, "target": target, "seeds": seeds, "recorded_mean": float(theirs.mean()),
            "recorded_sd": rec_sd, "ours_mean": float(ours.mean()),
            "ours_sd": float(ours.std(ddof=1)) if len(ours) > 1 else 0.0, "diff_mean": float(d.mean()),
            "diff_sd": sd_d, "diff_se": sd_d / np.sqrt(len(d)),
            "paired_z": float(d.mean() / (sd_d / np.sqrt(len(d)))) if sd_d > 0 else None,
            "within_seed_noise": bool(abs(d.mean()) <= rec_sd)}


@app.command()
def score(
    submission_file: Path = typer.Argument(..., help="built .h5ad"),
    target: Path = typer.Option(..., help="knockout as released (e.g. the Mab21l2 KO)"),
    ref: Path = typer.Option(..., help="matched WT as released (e.g. WT E9.5)"),
    board: str = typer.Option(..., help="board key (its gene panel; anchors are computed, never read)"),
    seeds: str = typer.Option("0-9", help="seeds (e.g. 0-9)"),
    records: Path = typer.Option(None, help="the run's records CSV, to compare with a recorded variant"),
    variant: str = typer.Option(None, help="recorded variant to compare with (e.g. fac:factor)"),
    record_target: str = typer.Option(None, help="target name in the records (default: the only one)"),
    out: Path = typer.Option(None, help="write the scores as JSON"),
    boards_dir: Path = typer.Option(None, help="boards directory"),
):
    """Score a built file with the local emulator (vec.evaluate) over seeds."""
    import anndata as ad

    from vec.boards import Board
    from vec.evaluate import Stage, score_board
    from vec.io import dense
    from vec.prep import load_stage

    try:
        b = board_spec(board, boards_dir)
        a = ad.read_h5ad(submission_file)
        if [str(g) for g in a.var_names] != b.genes:
            raise Refused("the file's genes are not the board's panel in order")
        X = dense(a.X).astype(np.float32)
        tX, rX = (load_stage(p, genes=b.genes, lineage=False, drop_duplicates=False).X.toarray().astype(np.float32)
                  for p in (target, ref))
        emu = Board(key=f"{b.key}:emulated", task=b.task, genes=b.genes, needs_coords=False)
        per = []
        for s in t3d.parse_seeds(seeds):
            r = score_board(emu, Stage(X=X), Stage(X=tX), Stage(X=rX), seed=s)
            per.append({"seed": s, "t3": r.score, **{f"skill_{k}": v for k, v in r.skills.items()},
                        **{k: r.metrics.get(k) for k in ("de_score", "de_direction", "severity_slope", "mmd_u",
                                                         "variogram")}})
        t3s = np.array([p["t3"] for p in per])
        res = {"command": "score", "file": str(submission_file), "file_sha256": sha256(submission_file),
               "target": str(target), "ref": str(ref), "board": b.key, "cells": int(X.shape[0]),
               "t3_mean": float(t3s.mean()), "t3_sd": float(t3s.std(ddof=1)) if len(t3s) > 1 else 0.0,
               "per_seed": per}
        if records and variant:
            res["comparison"] = compare_with_records({p["seed"]: p["t3"] for p in per}, records, variant,
                                                     record_target)
    except Refused as e:
        say(f"refused: {e}", "red")
        raise typer.Exit(2) from None
    say(f"{submission_file}: T3 {res['t3_mean']:.3f} +/- {res['t3_sd']:.3f} over {len(per)} seeds "
        f"({res['cells']} cells)")
    c = res.get("comparison")
    if c:
        say(f"  recorded {c['variant']}: {c['recorded_mean']:.3f} +/- {c['recorded_sd']:.3f}; paired difference "
            f"{c['diff_mean']:+.3f} (SE {c['diff_se']:.3f}); within seed noise: {c['within_seed_noise']}",
            "green" if c["within_seed_noise"] else "yellow")
    if out:
        dump_json(res, out)


if __name__ == "__main__":
    app()
