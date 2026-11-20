#!/usr/bin/env python
"""Outside-knockout tables in the barred factor's form: one implementation shared by the practice read
(k22_outside_forms.py, row K22) and the thin builder (build_outside_form.py), so that a file built from a table applies
exactly the factor its practice arm read.

A table enters only by NAME. Its bytes are found under the folder the environment variable VEC_T3_EXT_TABLES names
(combine.py's --out folder, outside every checkout), at <dir>/<name>/factors/<name>.csv, and are accepted only when
their sha256 equals the line of the committed listing experiments/t3_ext/TABLES.sha256 and the pins below (filled in
the commit that places row K22). No path is ever taken from a command line.

The listing, experiments/t3_ext/TABLES.sha256: one line per pooled table, sorted by name, `<64 lowercase hex>  <name>`
(two spaces), LF line ends, a final newline, no blank or comment line; a name matches ^[A-Za-z0-9_]+$ and never ends
in `_per_group`. It must be tracked by git. `python -I -B experiments/t3_decomp/outside_forms.py listing` writes it
from COMBINE.json (each file hashed again); `... pins` prints the two constant lines for this module from it.

load_table(name, genes) reads the file's bytes ONCE, hashes them and parses the same buffer. It refuses: a name not in
TRIES, a `_per_group` name, a name not listed or listed under another sha256 than the pin, an unset, relative,
symlinked or in-repository folder (t3x_common.outside_repo), a symlinked or missing table folder, factors folder,
csv, factors.json, manifest or COMBINE.json, bytes whose sha256 differs from the listing, a COMBINE.json or
factors/factors.json or manifest that records another sha256 (all three must exist), columns other than gene,factor,
a gene list other than `genes` (order included, as the frozen build_submission.load_factor checks phi_M), duplicated
genes, and values that are not finite and positive. The values are used as loaded: t3x_tables.finish already
rescaled them, clipped them to [0.25, 4] and set absent genes to exactly 1.

product(phi_m, x, T) = phi_m ** x * T ** EXPONENT (EXPONENT 1.0; phi_m itself at x = 1, T itself at 1.0); nothing is
rescaled or clipped afterwards (as build_two_knockout.py, factor_forms.py and decision log 2026-11-20 (q)'s C1).
permuted(T, basis) is harness.t3_decile_permutation(T, basis, PERM_SEED = 20261113), 10 bins, as K18 and K21 permute.
pregate(...) gives the target-free figures and the drops fixed in row K22.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
T3_EXT = REPO / "experiments" / "t3_ext"
TABLES_FILE = T3_EXT / "TABLES.sha256"
ENV = "VEC_T3_EXT_TABLES"
TRIES = ("all", "unrelated")        # the two tables K22 may read, also the tie order
EXPONENT = 1.0                      # the table's exponent, fixed from the tables' form (row K22)
PERM_SEED = 20261113                # k7_wt_revert.PERM_SEED (a test pins the equality)
PERM_BINS = 10                      # harness.K2_PERM_BINS (a test pins the equality)
CLIP = (0.25, 4.0)                  # t3x_tables.finish's clip
CLIP_TOL = 1e-9
MIN_CHANGED = 10                    # factor_forms.MIN_CHANGED (a test pins the equality)
MIN_STUDIES = 2                     # a try's table must hold at least two studies (the sign rule applied)
CARDIAC = "cardiac"                 # a table whose studies are all of this tier is dropped
MAX_ABS_CORR = 0.5                  # |corr(ln phi_M, ln T)| over T's genes != 1 at or above this drops the try
FILE2_MAX_CORR = 0.9                # file 2 is named only when corr(ln T_all, ln T_unrel) is below this
SINGLE_TIERS = ("wnt", "cardiac")   # 'all' is dropped when its bytes equal one of these tables or a study_ table
# the frozen T3 files' bytes, as build_two_knockout.py, build_barred_exponent.py and build_factor_form.py pin them
FROZEN_SHA256 = {
    "t3_decompose.py": "3e671c3bb5ff57bba5ceec80ea36fcc739b1f111399a1e2854400b5c9f5133d2",
    "rules.yaml": "cd2d94aa6e202529a2eead2a7ad60b312d9a2fd9ff0bd989d9320662de750c13",
    "build_submission.py": "1194a63993a27dde4dedfdb9180de035c5f8f2a57aaa7d97b38e7a15ea9e0c74",
}

# Filled in the commit that places row K22, from the committed listing (`outside_forms.py pins` prints both lines).
# None: the listing's pin is not set, and every run, pre-gate and outside build refuses. A table pin of None means the
# table is not listed (its try is dropped by the pre-gate); a listed table must have its pin, and the pin its line.
K22_LISTING_SHA256: str | None = "884b294bb8ea15b4e97b8f9f3d72869c6f9a375c55a5ca11a6daf3c70b8d94ea"
K22_TABLE_SHA256: dict[str, str | None] = {"all": "9490f524e40e146480c775f85ae8829ec9c005607104ad8ace129c5c46ed6228", "unrelated": None}

_LINE = re.compile(r"([0-9a-f]{64})  ([A-Za-z0-9_]+)")
_NAME = re.compile(r"[A-Za-z0-9_]+")


class TableRefused(Exception):
    """A refusal of the whole run or build (exit status 2): pins, listing, folder, bytes or records disagree."""


class FormRefused(ValueError):
    """A form that cannot be made from these factors (the builder's distinct exit status)."""


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
# The listing


def tracked(path: Path) -> bool:
    """git ls-files --error-unmatch in the repository that holds `path` (no lock taken)."""
    path = Path(path)
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["GIT_OPTIONAL_LOCKS"] = "0"
    r = subprocess.run(["git", "-C", str(path.parent), "ls-files", "--error-unmatch", "--", path.name],
                       capture_output=True, env=env, check=False)
    return r.returncode == 0


def parse_listing(raw: bytes) -> dict[str, str]:
    """{name: sha256} from the listing's bytes; refuses any other form."""
    try:
        text_ = raw.decode("ascii")
    except UnicodeDecodeError:
        raise TableRefused("TABLES.sha256 is not plain ASCII") from None
    if not text_.endswith("\n") or "\r" in text_:
        raise TableRefused("TABLES.sha256 must end in a newline and use LF line ends only")
    out: dict[str, str] = {}
    names = []
    for i, ln in enumerate(text_[:-1].split("\n"), 1):
        m = _LINE.fullmatch(ln)
        if not m:
            raise TableRefused(f"TABLES.sha256 line {i} is not `<64 lowercase hex>  <name>`")
        sha, name = m.groups()
        if name.endswith("_per_group"):
            raise TableRefused(f"TABLES.sha256 line {i}: {name}: per-group tables are never listed")
        if name in out:
            raise TableRefused(f"TABLES.sha256 lists {name} twice")
        out[name] = sha
        names.append(name)
    if names != sorted(names):
        raise TableRefused("TABLES.sha256 is not sorted by name")
    return out


def check_pins(pin: str | None, pins: dict) -> None:
    if pin is None or any(k not in pins for k in TRIES):
        raise TableRefused("the pins K22_LISTING_SHA256 and K22_TABLE_SHA256 are not filled (outside_forms.py; they "
                           "go in with row K22's commit, from the committed TABLES.sha256): refused")


def read_listing(path: Path | None = None, pin: str | None = ..., pins: dict | None = None,
                 check_tracked: bool = True) -> tuple[dict[str, str], str]:
    """(listing, its sha256), checked: the file exists and is no symlink, is tracked by git, hashes to the pin, has the
    fixed form, and agrees with the table pins for each name of TRIES (listed iff pinned, at the pinned sha256)."""
    path = TABLES_FILE if path is None else Path(path)
    pin = K22_LISTING_SHA256 if pin is ... else pin
    pins = K22_TABLE_SHA256 if pins is None else pins
    check_pins(pin, pins)
    if path.is_symlink() or not path.is_file():
        raise TableRefused(f"{path.name} is missing or a symlink: refused")
    if check_tracked and not tracked(path):
        raise TableRefused(f"{path.name} is not tracked by git (it is committed before row K22): refused")
    raw = path.read_bytes()
    sha = sha256_bytes(raw)
    if sha != pin:
        raise TableRefused(f"{path.name} hashes to {sha[:12]}, not the pinned K22_LISTING_SHA256 {pin[:12]}: refused")
    listing = parse_listing(raw)
    for name in TRIES:
        if listing.get(name) != pins.get(name):
            raise TableRefused(f"{name}: the listing says {str(listing.get(name))[:12]}, the pin "
                               f"{str(pins.get(name))[:12]}: refused")
    return listing, sha


# ---------------------------------------------------------------------------------------------------------------------
# The folder and one table


def _cm():
    p = str(T3_EXT)
    if p not in sys.path:
        sys.path.insert(0, p)
    import t3x_common as cm

    return cm


def table_dir(env: dict | None = None) -> Path:
    """The tables folder from the environment: set, absolute, no symlink, an existing directory, outside this checkout,
    the main checkout and the DGX's reference checkout (t3x_common.outside_repo). Returns its realpath."""
    env = os.environ if env is None else env
    v = env.get(ENV)
    if not v:
        raise TableRefused(f"the environment variable {ENV} is not set: refused (no table path is taken otherwise)")
    p = Path(v)
    if not p.is_absolute():
        raise TableRefused(f"{ENV} is not an absolute path: refused")
    if p.is_symlink():
        raise TableRefused(f"{ENV} names a symlink: refused")
    if not p.is_dir():
        raise TableRefused(f"{ENV} does not name a directory: refused")
    cm = _cm()
    try:
        return cm.outside_repo(Path(os.path.realpath(p)))
    except cm.Refused as e:
        raise TableRefused(f"{ENV}: {e}") from None


def _plain_file(p: Path, what: str) -> Path:
    if p.is_symlink() or not p.is_file():
        raise TableRefused(f"{what} ({p.name}) is missing or a symlink: refused")
    return p


def _json(p: Path, what: str) -> tuple[dict, str]:
    raw = _plain_file(p, what).read_bytes()
    try:
        return json.loads(raw), sha256_bytes(raw)
    except ValueError:
        raise TableRefused(f"{what} is not JSON: refused") from None


def combine_record(d: Path) -> tuple[dict, str]:
    """COMBINE.json's pooled sha256 per table name ({name: sha}) and the file's own sha256."""
    c, sha = _json(Path(d) / "COMBINE.json", "COMBINE.json")
    out = {}
    for name, files in (c.get("tables") or {}).items():
        rec = files.get(name) if isinstance(files, dict) else None
        if not isinstance(rec, dict) or not isinstance(rec.get("sha256"), str):
            raise TableRefused(f"COMBINE.json holds no pooled record for {name}: refused")
        out[name] = rec["sha256"]
    return out, sha


def check_name(name: str) -> None:
    if not isinstance(name, str) or not _NAME.fullmatch(name) or name.endswith("_per_group"):
        raise TableRefused(f"{name!r} is not a pooled table name: refused")
    if name not in TRIES:
        raise TableRefused(f"{name!r} is not one of K22's tables {list(TRIES)}: refused")


def load_table(name: str, genes: list[str], listing: dict[str, str], listing_sha256: str,
               pins: dict | None = None, env: dict | None = None) -> tuple[np.ndarray, dict]:
    """The named pooled table in `genes`' order (which must be the file's order), with its record; see the module."""
    import pandas as pd

    check_name(name)
    pins = K22_TABLE_SHA256 if pins is None else pins
    if name not in listing:
        raise TableRefused(f"{name} is not listed in TABLES.sha256: refused")
    want = listing[name]
    if pins.get(name) != want:
        raise TableRefused(f"{name}: the listing's sha256 is not the pinned one: refused")
    d = table_dir(env)
    td, fd = d / name, d / name / "factors"
    for p, what in ((td, "the table folder"), (fd, "the factors folder")):
        if p.is_symlink() or not p.is_dir():
            raise TableRefused(f"{name}: {what} is missing or a symlink: refused")
    path = _plain_file(fd / f"{name}.csv", f"{name}'s pooled table")
    raw = path.read_bytes()  # read once: the bytes hashed are the bytes parsed
    sha = sha256_bytes(raw)
    if sha != want:
        raise TableRefused(f"{name}: the file hashes to {sha[:12]}, not the listed {want[:12]}: refused")
    comb, comb_sha = combine_record(d)
    if comb.get(name) != sha:
        raise TableRefused(f"{name}: COMBINE.json records another sha256: refused")
    fj, _ = _json(fd / "factors.json", f"{name}'s factors.json")
    if ((fj.get("files") or {}).get(name) or {}).get("sha256") != sha:
        raise TableRefused(f"{name}: factors/factors.json records another sha256: refused")
    man, man_sha = _json(td / f"{name}_manifest.json", f"{name}'s manifest")
    if man.get("name", name) != name or ((man.get("files") or {}).get(name) or {}).get("sha256") != sha:
        raise TableRefused(f"{name}: its manifest records another table or sha256: refused")
    try:
        df = pd.read_csv(io.BytesIO(raw), float_precision="round_trip")
    except Exception:  # noqa: BLE001
        raise TableRefused(f"{name}: not a readable CSV: refused") from None
    if list(df.columns) != ["gene", "factor"]:
        raise TableRefused(f"{name}: columns {list(df.columns)}, not ['gene', 'factor']: refused")
    got = list(df["gene"].astype(str))
    if len(set(got)) != len(got):
        raise TableRefused(f"{name}: duplicated genes: refused")
    if got != list(genes):
        raise TableRefused(f"{name}: genes differ from the barred factor's genes (order included): refused")
    try:
        phi = df["factor"].to_numpy(dtype=np.float64)
    except (TypeError, ValueError):
        raise TableRefused(f"{name}: factor values are not numbers: refused") from None
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise TableRefused(f"{name}: factor values must be finite and positive: refused")
    studies = [str(s) for s in man.get("studies") or []]
    tiers = sorted({str(t) for t in man.get("tiers") or []})
    rec = {"name": name, "file": f"<{ENV}>/{name}/factors/{name}.csv", "dir_realpath": str(d), "sha256": sha,
           "listing_sha256": listing_sha256, "pinned_sha256": pins.get(name), "combine_sha256": comb_sha,
           "manifest_sha256": man_sha, "studies": studies, "tiers": tiers, "sign_rule": man.get("sign_rule"),
           "refstats_sha256": man.get("refstats_sha256"), "checked": man.get("checked")}
    return phi, rec


def load_listed(genes: list[str], path: Path | None = None, pin: str | None = ..., pins: dict | None = None,
                env: dict | None = None, check_tracked: bool = True) -> dict:
    """The listing checked against its pins, then every listed table of TRIES loaded (a listed table that fails a
    check refuses the whole run), and COMBINE.json's pooled sha256 of every table (for the pre-gate's byte rule).
    Returns {'listing_sha256', 'tables': {name: (T, record)}, 'combine': {name: sha} or {}, 'combine_sha256'}."""
    pins = K22_TABLE_SHA256 if pins is None else pins
    listing, lsha = read_listing(path, pin, pins, check_tracked)
    tables = {n: load_table(n, genes, listing, lsha, pins, env) for n in TRIES if n in listing}
    comb, csha = (combine_record(table_dir(env)) if tables else ({}, None))
    refs = {r["refstats_sha256"] for _, r in tables.values()}
    if len(refs) > 1:
        raise TableRefused("the listed tables were made against different reference statistics: refused")
    return {"listing_sha256": lsha, "tables": tables, "combine": comb, "combine_sha256": csha}


def ext_record(loaded: dict) -> dict:
    """What every task file and sidecar records of the tables (the VOID check compares it across files)."""
    return {"listing_sha256": loaded["listing_sha256"], "combine_sha256": loaded["combine_sha256"],
            "tables": {n: r["sha256"] for n, (_, r) in loaded["tables"].items()}}


# ---------------------------------------------------------------------------------------------------------------------
# The forms


def _positive(v, what: str) -> np.ndarray:
    v = np.asarray(v, dtype=np.float64)
    if v.ndim != 1 or v.size == 0 or not np.all(np.isfinite(v) & (v > 0)):
        raise FormRefused(f"{what} must be a non-empty vector of finite positive values")
    return v


def product(phi_m, x: float, T, e: float = EXPONENT) -> np.ndarray:
    """phi_m ** x * T ** e, gene by gene, nothing rescaled or clipped afterwards; phi_m (T) itself at x = 1 (e = 1)."""
    phi_m = _positive(phi_m, "phi_M")
    T = _positive(T, "the outside table")
    if phi_m.shape != T.shape:
        raise FormRefused("phi_M and the outside table do not cover the same genes")
    out = (phi_m if x == 1.0 else phi_m ** x) * (T if e == 1.0 else T ** e)
    if not np.all(np.isfinite(out) & (out > 0)):
        raise FormRefused("the product is not finite and positive")
    return out


def permuted(T, basis) -> np.ndarray:
    """The control's table: harness.t3_decile_permutation within deciles of `basis`, seed PERM_SEED, PERM_BINS bins."""
    from vec import harness

    return harness.t3_decile_permutation(np.asarray(T, dtype=np.float64), basis, PERM_SEED, PERM_BINS)


def _at_clip(phi: np.ndarray) -> np.ndarray:
    return np.abs(np.log(phi)) >= math.log(CLIP[1]) - CLIP_TOL


def _corr(a: np.ndarray, b: np.ndarray) -> float | None:
    if a.size < 3 or not (np.std(a) > 0 and np.std(b) > 0):
        return None
    return float(np.corrcoef(a, b)[0, 1])


def figures(phi_m, phi, den=None) -> dict:
    """Target-free figures of a factor against phi_M (phi_m): genes changed, genes != 1, genes at the clip, SD and
    gene mean of ln, min, max, n_up, n_down and, with den (WT E8.75's mean within-cell proportions), the den-weighted
    total den . phi / sum den against phi_M's."""
    phi_m, phi = _positive(phi_m, "phi_M"), _positive(phi, "the factor")
    lp = np.log(phi)
    out = {"genes_changed": int((phi != phi_m).sum()), "genes_ne_1": int((phi != 1.0).sum()),
           "genes_at_clip": int(_at_clip(phi).sum()), "sd_ln": float(np.std(lp)), "gene_mean_ln": float(lp.mean()),
           "min": float(phi.min()), "max": float(phi.max()), "n_up": int((phi > 1).sum()),
           "n_down": int((phi < 1).sum())}
    if den is not None:
        den = np.asarray(den, dtype=np.float64)
        out["den_total"] = float(den @ phi / den.sum())
        out["den_total_phi_m"] = float(den @ phi_m / den.sum())
    return out


def table_figures(phi_m, T, den=None) -> dict:
    """A table's figures, plus corr and OLS slope of ln T on ln phi_M over T's genes != 1."""
    phi_m, T = _positive(phi_m, "phi_M"), _positive(T, "the outside table")
    out = figures(np.ones_like(T), T, den)
    out.pop("genes_changed")
    m = T != 1.0
    lm, lt = np.log(phi_m[m]), np.log(T[m])
    out["corr_ln_phi_m"] = _corr(lm, lt)
    v = float(np.var(lm)) if m.sum() >= 3 else 0.0
    out["ols_slope_ln_t_on_ln_phi_m"] = float(np.mean((lm - lm.mean()) * (lt - lt.mean())) / v) if v > 0 else None
    return out


def corr_between(T1, T2) -> float | None:
    """corr(ln T1, ln T2) over the genes != 1 in either table (file 2's rule)."""
    T1, T2 = _positive(T1, "table 1"), _positive(T2, "table 2")
    m = (T1 != 1.0) | (T2 != 1.0)
    return _corr(np.log(T1[m]), np.log(T2[m]))


def pregate(phi_m, loaded: dict, basis=None, den=None) -> dict:
    """The pre-gate of row K22, target-free, from phi_M and the loaded tables alone. The drops depend on neither
    `basis` nor `den` (both only add figures): a try's table is dropped (with its try, its control and, for 'all',
    rep_all) when it is not listed; when its manifest lists fewer than MIN_STUDIES studies or only cardiac-tier
    studies; ('all') when its bytes equal a study_ table's or the wnt or cardiac table's in COMBINE.json; when it holds
    fewer than MIN_CHANGED genes != 1 (a permutation moves values, not their count, so this covers the control too);
    when the product changes fewer than MIN_CHANGED genes against phi_M; when |corr(ln phi_M, ln T)| over its genes
    != 1 is MAX_ABS_CORR or more, or undefined; ('unrelated') when its bytes equal 'all''s and 'all' is not dropped
    for its bytes."""
    phi_m = _positive(phi_m, "phi_M")
    tables, comb = loaded["tables"], loaded.get("combine") or {}
    out: dict = {"n_genes": int(phi_m.size), "phi_m": figures(phi_m, phi_m, den), "tables": {}, "products": {},
                 "controls": {}, "dropped": {}, "corr_ln_all_ln_unrelated": None,
                 "constants": {"min_studies": MIN_STUDIES, "min_changed": MIN_CHANGED, "max_abs_corr": MAX_ABS_CORR,
                               "file2_max_corr": FILE2_MAX_CORR, "exponent": EXPONENT, "perm_seed": PERM_SEED,
                               "perm_bins": PERM_BINS}}
    for name in TRIES:
        if name not in tables:
            out["dropped"][name] = "not listed in TABLES.sha256"
            continue
        T, rec = tables[name]
        tf = table_figures(phi_m, T, den)
        tf |= {"studies": rec["studies"], "tiers": rec["tiers"], "sign_rule": rec["sign_rule"],
               "sha256": rec["sha256"]}
        out["tables"][name] = tf
        p = product(phi_m, 1.0, T)
        out["products"][name] = figures(phi_m, p, den)
        if basis is not None:
            out["controls"][name] = figures(phi_m, product(phi_m, 1.0, permuted(T, basis)), den)
        why = []
        if len(rec["studies"]) < MIN_STUDIES:
            why.append(f"its manifest lists {len(rec['studies'])} study (< {MIN_STUDIES}; no sign rule)")
        if rec["tiers"] and set(rec["tiers"]) <= {CARDIAC}:
            why.append("only cardiac-tier studies")
        if name == "all":
            same = sorted(k for k, v in comb.items() if v == rec["sha256"] and k != name
                          and (k.startswith("study_") or k in SINGLE_TIERS))
            if same:
                why.append(f"its bytes equal {same}")
        if tf["genes_ne_1"] < MIN_CHANGED:
            why.append(f"{tf['genes_ne_1']} genes != 1 (< {MIN_CHANGED})")
        if out["products"][name]["genes_changed"] < MIN_CHANGED:
            why.append(f"the product changes {out['products'][name]['genes_changed']} genes (< {MIN_CHANGED})")
        c = tf["corr_ln_phi_m"]
        if c is None or abs(c) >= MAX_ABS_CORR:
            why.append(f"|corr(ln phi_M, ln T)| {'undefined' if c is None else f'{abs(c):.3f}'} (>= {MAX_ABS_CORR})")
        if why:
            out["dropped"][name] = "; ".join(why)
    if "all" in tables and "unrelated" in tables:
        out["corr_ln_all_ln_unrelated"] = corr_between(tables["all"][0], tables["unrelated"][0])
        if tables["unrelated"][1]["sha256"] == tables["all"][1]["sha256"] and "unrelated" not in out["dropped"]:
            out["dropped"]["unrelated"] = "its bytes equal all's (one table, read once, under 'all')"
    return out


# ---------------------------------------------------------------------------------------------------------------------
# The listing from COMBINE.json and the pins (for the main session)


def listing_from_combine(env: dict | None = None) -> str:
    """TABLES.sha256's text for every pooled table COMBINE.json names, each file hashed again from disk (refused if a
    hash differs or a file is a symlink)."""
    d = table_dir(env)
    comb, _ = combine_record(d)
    lines = []
    for name in sorted(comb):
        if not _NAME.fullmatch(name) or name.endswith("_per_group"):
            raise TableRefused(f"COMBINE.json names {name!r}, not a pooled table name: refused")
        p = _plain_file(d / name / "factors" / f"{name}.csv", f"{name}'s pooled table")
        for q in (d / name, d / name / "factors"):
            if q.is_symlink():
                raise TableRefused(f"{name}: a folder on its path is a symlink: refused")
        sha = sha256_bytes(p.read_bytes())
        if sha != comb[name]:
            raise TableRefused(f"{name}: the file hashes to {sha[:12]}, COMBINE.json says {comb[name][:12]}: refused")
        lines.append(f"{sha}  {name}\n")
    return "".join(lines)


def pins_text(path: Path | None = None) -> str:
    path = TABLES_FILE if path is None else Path(path)
    raw = path.read_bytes()
    listing = parse_listing(raw)
    items = ", ".join(f'"{n}": ' + (f'"{listing[n]}"' if n in listing else "None") for n in TRIES)
    return (f'K22_LISTING_SHA256: str | None = "{sha256_bytes(raw)}"\n'
            f"K22_TABLE_SHA256: dict[str, str | None] = {{{items}}}\n")


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    try:
        if argv == ["listing"]:
            sys.stdout.write(listing_from_combine())
            return 0
        if argv == ["pins"]:
            sys.stdout.write(pins_text())
            return 0
    except TableRefused as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    print("usage: outside_forms.py listing | pins (no other argument; the folder comes from "
          f"{ENV})", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
