#!/usr/bin/env python
"""The barred branch's T3 file with its Mab21l2 factor multiplied by one outside-knockout table: the forms row K22's
reading may name (ledger row K22 and the decision-log entry that fixes its names before the run). A file is built only
in a form the committed K22 reading names, on the base B that decision log 2026-11-20 (q) fixes; nothing here chooses
a form, a table or a base.

A thin builder, in the manner of build_factor_form.py and build_two_knockout.py (which it does not change). It imports
the frozen builder (build_submission.py) and script (t3_decompose.py) unchanged and builds a barred-branch file exactly
as `build_submission.py build DECISIONS --branch barred --submission s1|s2` does, with one change, for this process
only: the factor the frozen loader returns for `transfer_mab21l2_all` (phi_M, sha256-checked by the frozen loader
against the decision file and here against the pin) is replaced by

  identity-B   phi_M ** x                       the base file byte for byte (BASES' sha256), the launcher's gate;
                                                no outside table, folder or reading is needed
  xall-B       phi_M ** x * T_all ** 1.0        outside_forms.product, the function K22's arm x_all applies
  xunrel-B     phi_M ** x * T_unrel ** 1.0      the same, arm x_unrel

for B in s1, s1x1.5, s1x2.0, s2 (x = 1, 1.5, 2.0, 1). T is the pooled table loaded by outside_forms.load_table: by
name, from the folder of the environment variable VEC_T3_EXT_TABLES, at the sha256 of the committed
experiments/t3_ext/TABLES.sha256 and of the pins in outside_forms.py, with COMBINE.json, factors/factors.json and the
manifest agreeing, its genes the reference's (order included). Nothing is rescaled or clipped after the product;
t3_decompose.build receives it at w = 1 and applies it once (apply_gene_factor, after the cells are drawn).

Refused by name: the controls (xall_perm, xunrel_perm), rep_all and anything else. No number or path is taken from the
command line: seed 1000001 and n 'max' are the only values accepted.

Exit status 2 (refused, nothing built): another form; another seed or n; frozen files off their pins; a tree that is
not clean; the frozen decision-file check; a spec factor other than {mab21l2, 1.0}; phi_M off its pin; for an outside
form, pins not filled, TABLES.sha256 untracked or off its pin, the folder variable unset or refused, bytes or records
that disagree, or a committed K22 reading (experiments/t3_decomp/results/k22_outside_forms/reading.json, tracked by
git) that is missing, VOID, does not name the form's arm, or records other tables than the pins. Exit status 4: the
form's own refusal (outside_forms.FormRefused: the product cannot be made, or it changes fewer than 10 genes against
phi_M ** x); the launcher lists it, builds nothing in its place, and goes on.

--table-exponent E (row K23; default 1.0): only 1.0, 2.0 or 4.0 (else exit 2). At 1.0 nothing changes: the code path,
the file and the sidecar are as above. At 2.0 or 4.0 only the kind xall is built (identity and xunrel are refused), as

  xall-B --table-exponent E   phi_M ** x * T_all ** E   outside_forms.product at exponent E, the factor K23's arm x_eE
                                                        applies

and only when the committed K23 reading (experiments/t3_decomp/results/k23_outside_exponents/reading.json, tracked by
git) is not VOID or incomplete, names exactly x_eE, records the pinned tables, and holds E among its pre-gate's
exponents, not dropped; K22's reading is not used then. The form's own refusal (exit 4) is as above. The sidecar's
`outside_build` then records `exponent` E, `k23_reading`, `priced_at` 'x = 1 only' and K23's list entry, with
`k22_reading` None. K23's arms (x_e1, x_e2, x_e4, their controls and twins) are refused by name.

--captain-ruling R22 (row R22 of experiments/harness/TRIES.md, the captain's order of a score-guided search of round
steps; off by default: without it every byte and behaviour above is unchanged): the T3 search line's one file, xall-s2
at --table-exponent 4.0 (phi_M * T_all ** 4 on the barred s2), nothing else changed. Taken only on that form and
exponent, and only while TRIES.md holds exactly one row R22 whose first cell starts with "R22 (the captain's order",
says "ruling, no run" and names `--table-exponent 4.0` and `--captain-ruling R22` (each in backticks); any other row,
form or exponent with the flag is refused (exit 2). K23's reading is not used then (row K23 did not name exponent 4).
Before anything is built, the guard of the line's rule: from released WT E8.75 (load_stage, duplicates dropped; den =
its mean within-cell proportions), the pinned phi_M and the pinned table alone, the den-weighted total of the factor
must lie within GUARD_MAX_DEN_SHIFT (5%) of phi_M's (row K23's pre-gate criterion (iii)); else exit 5, nothing written.
The sidecar's `outside_build` then records `captain_ruling`, `line`, `priced_at` 'x = 1 only' with `priced_note`
(priced on practice at the base's exponent 1 only; row K23 did not name this exponent), `guard` (its figures) and
`k23_reading` None.

The sidecar is the frozen builder's plus one block `outside_build` (this script and its sha256, outside_forms.py's
sha256, the form, the base, x, the table's name, sha256, listing sha256, studies, tiers and register accessions, the
exponent 1.0, the product's figures, the pinned frozen hashes, the K22 reading's sha256, the ledger row). It adds no
cell number. `inputs.factor` stays the frozen loader's record of phi_M; the factor applied is the one `outside_build`
names, as `inputs.factor_note` says.

  VEC_T3_EXT_TABLES=~/v/vec-ext/t3/tables uv run --frozen --no-sync python -I -B \\
      experiments/t3_decomp/build_outside_form.py build DECISIONS.json --form xall-s1 --wt data/raw/T2/E8.75.h5ad \\
      --run-dir runs/t3_decomp/gata4_<STAMP> --board T3:bcatenin --seed 1000001 --out runs/<folder>/T3_xall_s1.h5ad
"""

from __future__ import annotations

import hashlib
import json
import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import typer

import vec.provenance  # noqa: F401  (records this build's start, before the method code below is loaded)

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import build_submission as bs  # noqa: E402
import outside_forms as ext  # noqa: E402  (the tables and the product, shared with k22_outside_forms.py)

from vec import provenance  # noqa: E402

SCRIPT = "experiments/t3_decomp/build_outside_form.py"
SCRIPT_PATH = Path(__file__).resolve()
BRANCH = "barred"
FROZEN_FACTOR = {"source": "mab21l2", "exponent": 1.0}
TABLE_M = "transfer_mab21l2_all"
TABLE_SHA256 = "64b1e0b8572c7193ed1e1f0b7339a26081e601f1e706e6b08193e1718dc5c7e8"  # gata4_decisions.json's record
FROZEN_SHA256 = dict(ext.FROZEN_SHA256)
BASES = {
    "s1": ("s1", 1.0, "3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2"),
    "s1x1.5": ("s1", 1.5, "c060fcaa6c7fdc6d6b1a58c8f850a1895a8295abbc50f297de71a32c97f9c420"),
    "s1x2.0": ("s1", 2.0, "7393810fbe3be2cde67cf23ed54f45d22da738d7d3c9e124780a8f6281075661"),
    "s2": ("s2", 1.0, "1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f"),
}
KINDS = {"identity": None, "xall": ("all", "x_all"), "xunrel": ("unrelated", "x_unrel")}  # kind -> (table, K22 arm)
FORMS = tuple(f"{k}-{b}" for k in KINDS for b in BASES)
REFUSED_KINDS = {"xall_perm": "a control of try K22, never built", "xunrel_perm": "a control of try K22, never built",
                 "rep_all": "the labelled control of the planner's doubt in try K22, never built",
                 "x_all": "use xall-<base>", "x_unrel": "use xunrel-<base>", "barred_form": "use identity-<base>",
                 "x_e1": "the reproduction guard of try K23, never built",
                 "x_e2": "use xall-<base> --table-exponent 2.0", "x_e4": "use xall-<base> --table-exponent 4.0",
                 "x_e2_perm": "a control of try K23, never built", "x_e4_perm": "a control of try K23, never built",
                 "x_e2_size": "a size twin of try K23, never built", "x_e4_size": "a size twin of try K23, never built"}
READING = HERE / "results" / "k22_outside_forms" / "reading.json"
READING_K23 = HERE / "results" / "k23_outside_exponents" / "reading.json"
TABLE_EXPONENTS = (1.0, 2.0, 4.0)  # --table-exponent: 1.0 (K22's forms, the default) or one of K23's exponents
K23_LIST_ENTRY = "experiments/harness/TRIES.md row K23 and decision log 2026-11-21 (a)"
K23_PRICED_AT = "x = 1 only"
# the outside datasets by register accession (experiments/t3_ext/datasets.yaml: study -> series), for the sidecar
STUDY_SERIES = {"GSE167493": ["GSE168768", "GSE186934"], "GSE137337": ["GSE137337"], "GSE204908": ["GSE204906"],
                "GSE205117": ["GSE205117"]}
LIST_ENTRY = "experiments/harness/TRIES.md row K22 and the decision-log entry that fixes its names"
LIST_SEED, LIST_N_CELLS = 1_000_001, "max"
LIST_INPUTS_REFUSAL = "the list fixes seed 1000001 and max_cells"
FACTOR_NOTE = "the frozen loader's record of phi_M only; the factor applied is outside_build.factor"
EXIT_REFUSED, EXIT_FORM_REFUSED = 2, 4
# Row R22's T3 search line (the captain's order of 2026-11-22; the line's rule text in docs/PLAN.md decision log
# 2026-11-22): the outside-table exponent on row K23's file, base B the barred s2, E in 0 (B), 2.0 (the K23 file), 4.0.
TRIES_FILE = HERE.parents[1] / "experiments" / "harness" / "TRIES.md"
R22_ROW = "R22"
R22_EXPONENT = 4.0
R22_FORM = "xall-s2"
R22_OPTION = "`--table-exponent 4.0`"
R22_LINE = ("T3 search line (R22): the outside-table exponent on row K23's file; base B the barred s2 (E = 0), the "
            "K23 file on B (E = 2.0), this file (E = 4.0); the line ends at 4.0, and 8.0 and 1.0 are never tried")
R22_LIST_ENTRY = "experiments/harness/TRIES.md row R22 and docs/PLAN.md decision log 2026-11-22 (the T3 search line)"
R22_PRICED_NOTE = ("priced on practice at the base's exponent 1 only (row K23's arms, x = 1); row K23 did not name "
                   "exponent 4")
GUARD_MAX_DEN_SHIFT = 0.05  # row K23's pre-gate criterion (iii), the line's guard
GUARD_TABLE = "all"
EXIT_GUARD_REFUSED = 5

Refused = bs.Refused


class FormRefused(Refused):
    """The form's own refusal (exit status 4)."""


class GuardRefused(Refused):
    """The R22 line's guard fails (exit status 5): nothing is built."""


app = typer.Typer(no_args_is_help=True, add_completion=False)


def parse_form(form) -> tuple[str, str]:
    """--form: exactly one of FORMS; returns (kind, base), else Refused (controls and rep_all by name)."""
    f = str(form) if form is not None else ""
    kind = f.split("-", 1)[0]
    if kind in REFUSED_KINDS:
        raise Refused(f"--form {form!r} is refused: {kind} is {REFUSED_KINDS[kind]}")
    if f not in FORMS:
        raise Refused(f"--form {form!r} is not one of {list(FORMS)}")
    kind, base = f.split("-", 1)
    return kind, base


def frozen_paths() -> dict[str, Path]:
    return {"t3_decompose.py": bs.T3D_PATH, "rules.yaml": bs.RULES_PATH, "build_submission.py": bs.BUILDER_PATH}


def check_frozen(paths: dict[str, Path] | None = None) -> dict[str, str]:
    paths = frozen_paths() if paths is None else paths
    got = {k: bs.sha256(p) for k, p in paths.items()}
    bad = [k for k in FROZEN_SHA256 if got.get(k) != FROZEN_SHA256[k]]
    if bad:
        raise Refused(f"frozen file(s) {bad} differ from the pinned bytes: not built")
    return got


def check_clean() -> dict:
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    if commit == "unknown" or unclean is None:
        raise Refused("git cannot tell the tree's state: not built (the stamp would not say clean)")
    if unclean:
        raise Refused(f"the tree is not clean ({len(unclean)} path(s), e.g. {unclean[:3]}): commit first (the build "
                      "stamp must say clean)")
    return {"commit": commit, "clean": True}


def check_reading(arm: str, table: str, path: Path | None = None, check_tracked: bool = True) -> dict:
    """The committed K22 reading names the arm, is not VOID, and records the pinned tables; returns its record."""
    path = READING if path is None else Path(path)
    if path.is_symlink() or not path.is_file():
        raise Refused(f"no K22 reading at {path.name}: an outside form is built only when the reading names it")
    if check_tracked and not ext.tracked(path):
        raise Refused(f"{path.name} is not committed: an outside form is built only from a committed reading")
    raw = path.read_bytes()
    r = json.loads(raw)
    if r.get("incomplete") or r.get("void") is not False:
        raise Refused("the K22 reading is VOID or incomplete: it names no file")
    if arm not in (r.get("names") or []):
        raise Refused(f"the K22 reading does not name {arm} (it names {r.get('names') or 'none'}): not built")
    e = r.get("ext") or {}
    pins = {n: v for n, v in ext.K22_TABLE_SHA256.items() if v is not None}
    if e.get("listing_sha256") != ext.K22_LISTING_SHA256 or e.get("tables") != pins or table not in pins:
        raise Refused("the K22 reading records other tables than the pins: not built")
    return {"file": str(path.relative_to(HERE.parents[1])) if path.is_relative_to(HERE.parents[1]) else str(path),
            "sha256": hashlib.sha256(raw).hexdigest(), "names": r["names"], "outcome": r.get("outcome")}


def check_exponent(e) -> float:
    """--table-exponent: one of TABLE_EXPONENTS, else Refused."""
    try:
        v = float(e)
    except (TypeError, ValueError):
        raise Refused(f"--table-exponent {e!r} is not one of {list(TABLE_EXPONENTS)}") from None
    if v not in TABLE_EXPONENTS:
        raise Refused(f"--table-exponent {e!r} is not one of {list(TABLE_EXPONENTS)}")
    return v


def check_reading_k23(e: float, path: Path | None = None, check_tracked: bool = True) -> dict:
    """The committed K23 reading names exactly x_eE, is not VOID, records the pinned tables and holds E among its
    pre-gate's exponents, not dropped; returns its record."""
    path = READING_K23 if path is None else Path(path)
    arm = f"x_e{e:g}"
    if path.is_symlink() or not path.is_file():
        raise Refused(f"no K23 reading at {path.name}: a form at exponent {e:g} is built only when the reading names "
                      "it")
    if check_tracked and not ext.tracked(path):
        raise Refused(f"{path.name} is not committed: a form at exponent {e:g} is built only from a committed reading")
    raw = path.read_bytes()
    r = json.loads(raw)
    if r.get("incomplete") or r.get("void") is not False:
        raise Refused("the K23 reading is VOID or incomplete: it names no file")
    if (r.get("names") or []) != [arm]:
        raise Refused(f"the K23 reading does not name {arm} (it names {r.get('names') or 'none'}): not built")
    rec = r.get("ext") or {}
    pins = {n: v for n, v in ext.K22_TABLE_SHA256.items() if v is not None}
    if rec.get("listing_sha256") != ext.K22_LISTING_SHA256 or rec.get("tables") != pins or "all" not in pins:
        raise Refused("the K23 reading records other tables than the pins: not built")
    pre = r.get("pregate") or {}
    exps = [float(v) for v in ((pre.get("constants") or {}).get("exponents") or [])]
    if e not in exps or f"{e:g}" in (pre.get("dropped") or {}):
        raise Refused(f"the K23 reading's pre-gate does not hold exponent {e:g} undropped: not built")
    return {"file": str(path.relative_to(HERE.parents[1])) if path.is_relative_to(HERE.parents[1]) else str(path),
            "sha256": hashlib.sha256(raw).hexdigest(), "names": r["names"], "outcome": r.get("outcome")}


def captain_order_names(row: str, option: str, path: Path | None = None) -> bool:
    """Whether TRIES.md (TRIES_FILE) holds exactly one row `row` whose first cell starts with "<row> (the captain's
    order", says "ruling, no run" and names `option` (in backticks) and "`--captain-ruling <row>`" (build_test.py's
    check of rows R20 and R22). A row that only mentions the option elsewhere does not pass."""
    p = TRIES_FILE if path is None else Path(path)
    if not p.is_file():
        return False
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {row} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    return (first.startswith(f"{row} (the captain's order") and "ruling, no run" in first and option in first
            and f"`--captain-ruling {row}`" in first)


def check_captain_ruling(ruling, form: str, e: float, path: Path | None = None) -> dict:
    """--captain-ruling: R22 only, on R22_FORM at R22_EXPONENT only, and only while row R22's first cell is the
    captain's order naming R22_OPTION; returns the record, else Refused."""
    if ruling != R22_ROW:
        raise Refused(f"--captain-ruling {ruling!r} is refused: only row {R22_ROW} orders a file of this builder")
    if e != R22_EXPONENT or form != R22_FORM:
        raise Refused(f"--captain-ruling {R22_ROW} takes --form {R22_FORM} --table-exponent {R22_EXPONENT:g} only "
                      f"(not --form {form} --table-exponent {e:g})")
    p = TRIES_FILE if path is None else Path(path)
    if not captain_order_names(R22_ROW, R22_OPTION, p):
        raise Refused(f"--table-exponent {R22_EXPONENT:g} --captain-ruling {R22_ROW} is refused: "
                      f"experiments/harness/TRIES.md has no single row {R22_ROW} whose first cell is the captain's "
                      f"order (\"ruling, no run\") naming {R22_OPTION} and `--captain-ruling {R22_ROW}`")
    return {"row": R22_ROW, "option": R22_OPTION, "ledger": "experiments/harness/TRIES.md",
            "ledger_sha256": bs.sha256(p)}


def wt_den(wt: Path) -> tuple[list[str], np.ndarray]:
    """WT E8.75's genes and mean within-cell proportions (den), as row K23's pre-gate made them: load_stage with
    duplicates dropped, dense counts, each cell's counts over its total, the mean over cells."""
    from vec.prep import load_stage

    a = load_stage(Path(wt), drop_duplicates=True)
    C = np.asarray(bs.t3d.dense_counts(a), dtype=np.float64)
    return list(a.var_names), (C / np.maximum(C.sum(1, keepdims=True), 1)).mean(0)


def exponent_guard(phi_m, x: float, T, e: float, den) -> dict:
    """The line's guard, a pure function of phi_M, the table and den: the den-weighted total of phi_M ** x * T ** e
    against that of phi_M ** x; passes when the relative shift is within GUARD_MAX_DEN_SHIFT."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    base = phi_m if x == 1.0 else phi_m ** x
    try:
        prod = ext.product(phi_m, x, T, e)
        f = ext.figures(base, prod, den)
    except ext.FormRefused as err:
        raise FormRefused(f"xall: {err}") from err
    shift = f["den_total"] / f["den_total_phi_m"] - 1.0
    return {"exponent": float(e), "x": float(x), "den_total_factor": f["den_total"],
            "den_total_base": f["den_total_phi_m"], "den_shift": float(shift), "max_den_shift": GUARD_MAX_DEN_SHIFT,
            "passed": bool(abs(shift) <= GUARD_MAX_DEN_SHIFT),
            "factor_sha256": hashlib.sha256(np.ascontiguousarray(prod, dtype=np.float64).tobytes()).hexdigest()}


def run_guard(decision_rec: dict, wt: Path, run_dir: Path | None, x: float, e: float, env: dict | None) -> dict:
    """The guard from released WT E8.75 and the pinned tables only, before any build: phi_M through the frozen loader
    (checked against TABLE_SHA256) in WT E8.75's gene order, the table `all` through outside_forms; GuardRefused when
    it fails."""
    genes, den = wt_den(wt)
    phi_m, rec = bs.load_factor(run_dir, TABLE_M, genes, decision_rec)
    if rec["sha256"] != TABLE_SHA256:
        raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                      f"({TABLE_SHA256[:16]}): not built")
    try:
        listing, lsha = ext.read_listing()
        T, trec = ext.load_table(GUARD_TABLE, list(genes), listing, lsha, env=env)
    except ext.TableRefused as err:
        raise Refused(str(err)) from None
    g = exponent_guard(phi_m, x, T, e, den)
    g |= {"wt": str(wt), "wt_sha256": bs.sha256(wt), "den": "WT E8.75's mean within-cell proportions (load_stage, "
          "duplicates dropped)", "table_m_sha256": rec["sha256"], "table": GUARD_TABLE, "table_sha256": trec["sha256"],
          "n_genes": len(genes), "genes": genes}
    if not g["passed"]:
        raise GuardRefused(f"the line's guard fails: the den-weighted total of phi_M ** {x:g} * T ** {e:g} moves "
                           f"{g['den_shift']:+.2%} against phi_M ** {x:g}'s (beyond {GUARD_MAX_DEN_SHIFT:.0%}): not "
                           "built, the line ends at 2.0")
    return g


def form_factor(kind: str, x: float, phi_m: np.ndarray, T: np.ndarray | None,
                e: float = 1.0) -> tuple[np.ndarray, dict]:
    """The vector applied: phi_M ** x (identity) or outside_forms.product(phi_M, x, T) (at e = 1, the default; else
    outside_forms.product(phi_M, x, T, e)); FormRefused where the product cannot be made or changes fewer than
    MIN_CHANGED genes against phi_M ** x."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    base = phi_m if x == 1.0 else phi_m ** x
    if kind == "identity":
        if e != 1.0:
            raise Refused("identity takes no table exponent")
        return base, {"applied": ext.figures(base, base)}
    try:
        phi = ext.product(phi_m, x, T) if e == 1.0 else ext.product(phi_m, x, T, e)
        figs = ext.figures(base, phi)
    except ext.FormRefused as e:
        raise FormRefused(f"{kind}: {e}") from e
    if figs["genes_changed"] < ext.MIN_CHANGED:
        raise FormRefused(f"{kind}: the product changes {figs['genes_changed']} genes against phi_M ** x "
                          f"(< {ext.MIN_CHANGED})")
    return phi, {"applied": figs, "table": ext.table_figures(phi_m, T)}


@contextmanager
def outside_loader(kind: str, x: float, out: dict, env: dict | None = None, exponent: float = 1.0,
                   guard: dict | None = None):
    """bs.load_factor, for this process only: the frozen loader for phi_M (checked against TABLE_SHA256 too), then, for
    an outside form, the table through outside_forms (by name, against the same genes), then the form's factor. With
    guard (the R22 line), the genes must be the guard's and the factor the one it passed."""
    original = bs.load_factor

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        if guard is not None and list(genes) != guard["genes"]:
            raise Refused("the build's genes are not the guard's (WT E8.75's): not built")
        phi_m, rec = original(run_dir, TABLE_M, genes, dec)
        if rec["sha256"] != TABLE_SHA256:
            raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                          f"({TABLE_SHA256[:16]}): not built")
        T, trec = None, None
        if kind != "identity":
            table = KINDS[kind][0]
            try:
                listing, lsha = ext.read_listing()
                T, trec = ext.load_table(table, list(genes), listing, lsha, env=env)
            except ext.TableRefused as e:
                raise Refused(str(e)) from None
        phi, frec = (form_factor(kind, x, phi_m, T) if exponent == 1.0 else
                      form_factor(kind, x, phi_m, T, exponent))
        if guard is not None and hashlib.sha256(np.ascontiguousarray(phi, dtype=np.float64).tobytes()).hexdigest() \
                != guard["factor_sha256"]:
            raise RuntimeError("the factor applied is not the one the guard passed")
        out.update(table_m=rec, table=trec, factor=frec, n_loads=out.get("n_loads", 0) + 1)
        return phi, rec

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


def build_form(decision: Path, form, wt: Path, board, n_cells: int, out: Path, run_dir: Path | None = None,
               seed: int = 1_000_001, overwrite: bool = False, threads: int = 16, env: dict | None = None,
               reading: Path | None = None, table_exponent: float = 1.0, k23_reading: Path | None = None,
               captain_ruling: str | None = None, tries: Path | None = None) -> dict:
    """The named form through bs.build_submission; returns the sidecar (the frozen builder's plus `outside_build`).
    table_exponent 1.0 (the default) is K22's path unchanged; 2.0 or 4.0 build xall at that exponent from the K23
    reading (k23_reading, default READING_K23). captain_ruling 'R22' (default None, today's paths) builds the R22
    line's xall-s2 at 4.0 under the ruling row of TRIES.md (tries, default TRIES_FILE) after the line's guard."""
    kind, base = parse_form(form)
    te = check_exponent(table_exponent)
    if te != 1.0 and kind != "xall":
        raise Refused(f"--table-exponent {te:g} is for the kind xall only (K23 reads `all`), not {kind}")
    ruled = None if captain_ruling is None else check_captain_ruling(captain_ruling, f"{kind}-{base}", te, tries)
    slot, x, base_sha = BASES[base]
    frozen = check_frozen()
    tree = check_clean()
    rrec = k23rec = None
    if kind != "identity":
        try:
            ext.check_pins(ext.K22_LISTING_SHA256, ext.K22_TABLE_SHA256)
        except ext.TableRefused as e:
            raise Refused(str(e)) from None
        if te == 1.0:
            rrec = check_reading(KINDS[kind][1], KINDS[kind][0], reading)
        elif ruled is None:
            k23rec = check_reading_k23(te, k23_reading)
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)
    spec = bs.barred_table(dec, rules)[slot]
    if spec["factor"] != FROZEN_FACTOR:
        raise Refused(f"the frozen barred {slot} factor is {spec['factor']}, not {FROZEN_FACTOR}: not built")
    out_p = Path(out)
    if ruled is not None:
        if not overwrite and (out_p.exists() or out_p.with_suffix(".json").exists()):
            raise Refused(f"{out_p} or its sidecar exists; pass --overwrite")
        guard = run_guard(dec, Path(wt), run_dir, x, te, env)
    got: dict = {}
    with (outside_loader(kind, x, got, env) if te == 1.0 else
          outside_loader(kind, x, got, env, te) if ruled is None else outside_loader(kind, x, got, env, te, guard)):
        side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                   run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
    if got.get("n_loads") != 1 or side["construction"]["factor_exponent"] != 1.0 or \
            side["spec"]["factor"] != FROZEN_FACTOR:
        raise RuntimeError("the build did not apply the form's factor once at w = 1")
    side["inputs"]["factor_note"] = FACTOR_NOTE
    trec = got["table"]
    side["outside_build"] = {
        "script": SCRIPT, "script_sha256": bs.sha256(SCRIPT_PATH),
        "outside_forms": "experiments/t3_decomp/outside_forms.py",
        "outside_forms_sha256": bs.sha256(HERE / "outside_forms.py"), "form": f"{kind}-{base}",
        "base": {"name": base, "slot": slot, "x": x, "reference_sha256": base_sha},
        "table_m": {"file": got["table_m"]["file"], "sha256": got["table_m"]["sha256"], "pinned_sha256": TABLE_SHA256},
        "table": None if trec is None else {
            **{k: trec[k] for k in ("name", "file", "dir_realpath", "sha256", "listing_sha256", "pinned_sha256",
                                    "combine_sha256", "manifest_sha256", "studies", "tiers", "sign_rule",
                                    "refstats_sha256")},
            "register_accessions": {s: STUDY_SERIES.get(s, [s]) for s in trec["studies"]}},
        "exponent": ext.EXPONENT, "factor": got["factor"], "frozen_files_sha256": frozen, "tree_at_start": tree,
        "k22_reading": rrec, "list_entry": LIST_ENTRY,
        "identity": f"the {base} base file byte for byte ({base_sha[:12]})" if kind == "identity" else None,
        "depth_step": None,
    }
    if te != 1.0 and ruled is None:
        side["outside_build"] |= {"exponent": te, "k23_reading": k23rec, "priced_at": K23_PRICED_AT,
                                  "list_entry": K23_LIST_ENTRY}
    elif ruled is not None:
        side["outside_build"] |= {"exponent": te, "k23_reading": None, "captain_ruling": ruled, "line": R22_LINE,
                                  "priced_at": K23_PRICED_AT, "priced_note": R22_PRICED_NOTE,
                                  "guard": {k: v for k, v in guard.items() if k != "genes"},
                                  "list_entry": R22_LIST_ENTRY}
    bs.dump_json(side, Path(out).with_suffix(".json"))
    return side


@app.callback()
def main() -> None:
    """The barred file with its Mab21l2 factor times one outside-knockout table, in one of try K22's fixed forms."""


@app.command()
def build(
    decision: Path = typer.Argument(..., help="decision file (experiments/t3_decomp/results/gata4/"
                                    "gata4_decisions.json)"),
    form: str = typer.Option(..., help=f"one of {', '.join(FORMS)}; nothing else"),
    wt: Path = typer.Option(..., help="WT reference (WT E8.75), as for the frozen builder"),
    board: str = typer.Option(..., help="T3 board key in boards/"),
    out: Path = typer.Option(..., help="output .h5ad (sidecar: same name with .json)"),
    run_dir: Path = typer.Option(None, help="run directory holding factors/"),
    seed: int = typer.Option(1_000_001, help="1000001, the list's seed; any other is refused"),
    n_cells: str = typer.Option("max", help="'max', the list's count: the board's max_cells, read in code, never "
                                "printed; any other is refused"),
    boards_dir: Path = typer.Option(None, help="boards directory (default: the repository's boards/)"),
    overwrite: bool = typer.Option(False, help="replace an existing output"),
    threads: int = typer.Option(16, help="BLAS threads while preparing the reference (the frozen default)"),
    table_exponent: float = typer.Option(1.0, help="the table's exponent: 1.0 (the default, K22's forms), or 2.0 or "
                                         "4.0 for xall when the committed K23 reading names that exponent"),
    captain_ruling: str = typer.Option(None, help="R22 only, with --form xall-s2 --table-exponent 4.0: the T3 search "
                                       "line's file, while row R22 of TRIES.md is the captain's order naming it"),
):
    """Build one fixed form."""
    try:
        if seed != LIST_SEED or n_cells != LIST_N_CELLS:
            raise Refused(LIST_INPUTS_REFUSAL)
        parse_form(form)
        te = check_exponent(table_exponent)
        b = bs.board_spec(board, boards_dir)
        if captain_ruling is not None:
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, table_exponent=te, captain_ruling=captain_ruling)
        elif te == 1.0:
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads)
        else:
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, table_exponent=te)
    except FormRefused as e:
        bs.say(f"refused (the form's own refusal): {e}", "red")
        raise typer.Exit(EXIT_FORM_REFUSED) from None
    except GuardRefused as e:
        bs.say(f"refused (the line's guard): {e}", "red")
        raise typer.Exit(EXIT_GUARD_REFUSED) from None
    except Refused as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(EXIT_REFUSED) from None
    chk = s["output"]["check"]
    ob = s["outside_build"]
    bs.say(f"form {ob['form']}: base {ob['base']['name']} (x = {ob['base']['x']}), genes changed against phi_M ** x "
           f"{ob['factor']['applied']['genes_changed']}")
    for w in s["warnings"]:
        bs.say(f"  warning: {w}", "yellow")
    for e in chk["errors"]:
        bs.say(f"  check error: {e}", "red")
    bs.say(f"{'OK' if chk['ok'] else 'FAILED the format check'}: {s['output']['h5ad']}",
           "green" if chk["ok"] else "red")
    raise typer.Exit(0 if chk["ok"] else 1)


if __name__ == "__main__":
    app()
