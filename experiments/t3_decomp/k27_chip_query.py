#!/usr/bin/env python
"""Try K27 (experiments/harness/TRIES.md, T3 candidate CHIPQ-S; the row is placed before any run, and every step
refuses without it): the knocked-out gene's own binding set, taken by code from a pinned generic ChIP database
(register row X018) in wild-type, untreated experiments, entered as a 500-gene table on top of row K23's named form
x_e2 = phi_M * two ** 2 at two committed sizes, read on the practice knockout (Gata4) against the named form and three
labelled controls per size, in K18's two practice worlds. A labelled reading for a ruling, not the adoption rule.

Cloned from k24_wider_tables.py. It imports outside_forms.py, k22_outside_forms.py, k23_outside_exponents.py,
k18_built_files.py, k21_factor_forms.py, build_outside_form.py, experiments/t3_ext/chipq.py and the harness unchanged
and, inside one process only, replaces k18's ARMS, build_arm and fold_factors, restored on exit, as rows K21 to K26 do.

THE TABLES enter only by KEY through the committed listing experiments/t3_ext/TABLES_K27.sha256 (11 lines
`<sha256>  <name>`, sorted by name: chipq_Gata4, chipq_rand_seed0 to chipq_rand_seed9) and the pins below, from the
folder the variable VEC_T3_EXT_TABLES_K27 names (chipq.tables_folder's checks). At the row's registration commit the
listing is tracked and EMPTY and every pin is None; the row's second commit (after P2, before P3) fills the listing's 11
lines, K27_LISTING_SHA256 and the 11 table pins (`pins --runs R` prints all three from P2's record). `p2` refuses when
the pins are filled; every counted step (`check`, `pregate`, `run`, `guard`) refuses (exit 2) while a pin is empty,
while the listing is absent or differs, or while a table's bytes differ from the listing, the pins, its manifest, the
fetch manifest it records or the barred factor's gene order. `two` is K22's table `all` through K22's listing and pins
(outside_forms.load_listed, K23's code path) from the folder VEC_T3_EXT_TABLES_TWO names.

Arms (phi_M = transfer_mab21l2_all as loaded; x_e2 = phi_M * two ** 2.0; for each size s in {S*, 4 S*}, keys s1 and s4,
e = S / SD(ln C), population SD over the 500 genes in file order; perm = harness.t3_decile_permutation within deciles of
each world's basis, seed 20261113, 10 bins, of the scaled table):
  barred_form     phi_M                          guard 1, reported (K18's arm value for value)
  x_e2            phi_M * two ** 2.0             guard 2 and the base of every paired difference (never named here)
  chip_q_<s>      x_e2 * C_Gata4 ** e            the try
  chip_rand_<s>   x_e2 * C_p(seed) ** e_p(seed)  control: seed s's placebo (chipq.placebo_draws), at the same size
  chip_perm_<s>   x_e2 * perm(C_Gata4 ** e)      control
  chip_size_<s>   k * phi_M ** a                 control: a = SD(ln P) / SD(ln phi_M), k = den_total(P) /
                                                 den_total(phi_M ** a), P = chip_q_<s>'s product, from P3
Nothing is rescaled or clipped after any product.

Steps (`python -I -B experiments/t3_decomp/k27_chip_query.py <step>`; the paths are the run folder and the data and run
folders the earlier rows take, nothing else):
  p2 --out R                 P2, database only, uncounted, at the registration commit: Gata4's bars, the pool, the
                             tables (chipq.tables), R/p2.json; prints counts and the listing lines only
  pins --runs R              the listing's text, K27_LISTING_SHA256 and K27_TABLE_SHA256 from R/p2.json, the tables
                             hashed again from disk (for the second commit)
  check                      the 11 tables against the listing, the pins, their manifests and the gene order
  pregate --out R            P3, COUNTED (the eight tries count from its start): R/pregate.json and its lines (no score)
  run --part guard|score --world plain --task g1|g2 [--seeds] --out R    as row K24
  run --part guard|score --world dm [--workers 8 --threads 2] --out R
  guard --runs R             the guard check: R/guard.json (exit 6 when a guard fails)
  read --runs R --out DIR    the fixed reading (never over an existing one: exit 4)

Exit codes: 0 done; 1 a crash (masked traceback) or a world whose chip_perm figures are not the pre-gate's (INCOMPLETE;
not a reading, the eight tries stand counted from P3); 2 refused (no row, a dirty tree, pins unfilled for a counted step
or filled for P2, the listing off its pin, a duplicated placebo table or a sha256 twice among the 11 listing lines, the
kept list absent, uncommitted or not P2's, a folder variable unset, bytes or records that disagree, a gene order
mismatch, a scaled SD off S by more than a relative 1e-12, the frozen files off their pins, `two` off K23's figures, a
scoring part without a passing guard.json); 3 P2: a bar not met (and, as row K24, guard and read: a file missing,
INCOMPLETE); 4 read: the reading exists; 5 pregate: the arm is dropped on either correlation; 6 guard: a guard failed.

No log line, message or result file carries a gene symbol of a bound set, a score of the database or the name of a
placebo antigen; no score in any log.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import traceback
from contextlib import contextmanager
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
T3X = REPO / "experiments" / "t3_ext"
for _p in (HERE, T3X):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
import build_outside_form as bof  # noqa: E402  (phi_M's pinned sha256)
import chipq as cq  # noqa: E402  (the database part: kept list, tables, pool, draws; unchanged)
import k7_wt_revert as k7  # noqa: E402  (the paired summary)
import k8_s2_form as k8  # noqa: E402  (the upload's central composition, for the record)
import k18_built_files as k18  # noqa: E402  (the worlds, tasks, stream, seeds, build path and scoring; unchanged)
import k21_factor_forms as k21  # noqa: E402  (the seed parser, gated count and K18's committed rows; unchanged)
import k22_outside_forms as k22  # noqa: E402  (frozen-file check, stamp check, k18's originals; unchanged)
import k23_outside_exponents as k23  # noqa: E402  (x_e2's code path, its pre-gate and constants; unchanged)
import outside_forms as ext  # noqa: E402  (the listing parser, the product, the permutation, figures; unchanged)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

cf = cq.cf

SCRIPT = "experiments/t3_decomp/k27_chip_query.py"
ROW_SCRIPT = f"`{SCRIPT}`"  # the ledger row names this script
TRIES_MD = REPO / "experiments" / "harness" / "TRIES.md"  # looked up at call time (row_present)
LISTING_FILE = T3X / "TABLES_K27.sha256"
ENV_TABLES = cq.ENV_TABLES  # "VEC_T3_EXT_TABLES_K27"
ENV_TWO = "VEC_T3_EXT_TABLES_TWO"
QUERY = cq.PRACTICE_QUERY  # "Gata4"
QUERY_TABLE = cq.TABLE_QUERY.format(q=QUERY)  # "chipq_Gata4"
PLACEBO_TABLE = {int(s): cq.TABLE_PLACEBO.format(s=s) for s in SEEDS}  # seed -> "chipq_rand_seed<s>"
TABLE_NAMES = tuple(sorted([QUERY_TABLE, *PLACEBO_TABLE.values()]))  # the listing's 11 names, sorted
BARRED_GENE = "Mab21l2"  # the barred factor's gene (transfer_mab21l2_all); P2 records only whether it is an antigen

# THE PINS. Empty (None) at the row's registration commit; the row's second commit, made after P2 and before P3, fills
# K27_LISTING_SHA256 and the 11 table pins from the committed TABLES_K27.sha256 (`pins --runs R` prints these lines).
K27_LISTING_SHA256: str | None = None
K27_TABLE_SHA256: dict[str, str | None] = {
    "chipq_Gata4": None,
    "chipq_rand_seed0": None,
    "chipq_rand_seed1": None,
    "chipq_rand_seed2": None,
    "chipq_rand_seed3": None,
    "chipq_rand_seed4": None,
    "chipq_rand_seed5": None,
    "chipq_rand_seed6": None,
    "chipq_rand_seed7": None,
    "chipq_rand_seed8": None,
    "chipq_rand_seed9": None,
}

STREAM = k18.STREAM  # "t3_s2_form"
BASE = k18.BARRED  # "barred_form": guard 1, reported
REPRO = "x_e2"  # K23's named form: guard 2 and the base of every paired difference
SIZE_KEYS = ("s1", "s4")  # S* and 4 S*
SIZE_OF_KEY = {"s1": "S*", "s4": "4S*"}  # chipq.SIZES' keys
SIZE_SHOWN = {"s1": "S*", "s4": "4 S*"}  # as the row writes them
TRY_KINDS = ("chip_q", "chip_rand", "chip_perm", "chip_size")
CONTROLS = ("chip_rand", "chip_perm", "chip_size")
ARMS = (BASE, REPRO, *(f"{k}_{s}" for s in SIZE_KEYS for k in TRY_KINDS))
GUARDS = (BASE, REPRO)
TRY_ARMS = tuple(a for a in ARMS if a not in GUARDS)
N_TRIES = 8  # the eight tries, counted from the start of P3
MAX_TRIES = 10  # plus one for each failed guard arm
TASKS = k18.TASKS
WORLDS = k18.WORLDS

# The row's constants (none new; none from a score)
S_STAR = cq.S_STAR  # 0.04949423311987974, K23's committed SD(ln two ** 2) (a test pins the equality with K24's)
SIZES = {"s1": cq.SIZES["S*"], "s4": cq.SIZES["4S*"]}
REL_TOL = 1e-12  # P2: SD(ln C ** e) equals S within this relative tolerance
P2_MIN_GENES_QUERY = 10  # Gata4's table differs from 1 at 10 genes or more
P2_MIN_GENES_PLACEBO = 1  # each placebo table differs from 1 at one gene or more
MAX_ABS_CORR = 0.5  # P3: abs corr at or above this, or undefined, drops the arm and its controls at both sizes
MIN_GAIN, MIN_AHEAD, MAX_EXTRA_GATED = k23.MIN_GAIN, k23.MIN_AHEAD, k23.MAX_EXTRA_GATED  # 0.3, 8 of 10, +1
GUARD_MEAN, GUARD_SEED = k23.GUARD_MEAN, k23.GUARD_SEED  # 0.1 and 0.5 (K17's tolerance)
DE_TERMS = ("de_score", "de_direction")
# K24's table bounds: RECORDED, NOT APPLIED (row K27)
K24_BOUNDS = {"max_abs_ln_scaled": 0.5, "max_sd_ratio": k23.MAX_SD_RATIO, "max_den_shift": k23.MAX_DEN_SHIFT,
              "max_at_clip": k23.MAX_AT_CLIP, "e_range": [0.125, 8.0]}
K18_TASKS = k23.K18_TASKS
K18_SHA256 = dict(k23.K18_SHA256)
K23_RESULTS = k23.RESULTS
RESULTS = HERE / "results" / "k27_chip_query"
EXIT_CRASH, EXIT_REFUSED, EXIT_BAR, EXIT_INCOMPLETE, EXIT_EXISTS, EXIT_DROP, EXIT_GUARD = 1, 2, 3, 3, 4, 5, 6
_K27: dict = {}  # the loaded tables (file order), `two`, the exponents and the twins' constants, set before any arm

# THE SENTENCES THE READING PRINTS (row K27), one of them
S_P1 = "Try K27: STOP at P1: the database bars are not met; no fold read, no try, no file."
S_P2 = "Try K27: STOP at P2: the practice bars or the placebo pool are not met; no fold read, no try, no file."
S_DROP = ("Try K27: dropped at the pre-gate: the table is not independent of the barred factor or of the named form's "
          "cross-knockout table; no file; tries counted: 8.")
S_VOID = "Try K27: VOID: {clause}; no file; tries counted: 8 (plus the failed guard arms)."
S_FAIL = ("Try K27: no file: the binding table does not qualify against the named form and its controls at either "
          "size; the K23 file stands alone; tries counted: 8.")
S_PASS = ("Try K27: PASS at size {size}: one candidate file is named for the captain and nothing is built: "
          "phi_M ** x * two ** 2 * C_Ctnnb1 ** e at size {size} on base B, e computed by the same formula from "
          "Ctnnb1's own table and recorded to full precision; the file exists only after the captain's ruling; "
          "tries counted: 8.")


class Refused(Exception):
    """A refusal of the whole run (exit status 2)."""


def refuse(msg: str) -> None:
    print(f"refused: {msg}", file=sys.stderr, flush=True)
    sys.exit(EXIT_REFUSED)


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def arm(kind: str, s: str) -> str:
    return f"{kind}_{s}"


def size_of_arm(a: str) -> str | None:
    return a.rsplit("_", 1)[1] if a in TRY_ARMS else None


def kind_of_arm(a: str) -> str | None:
    return a.rsplit("_", 1)[0] if a in TRY_ARMS else None


# ---------------------------------------------------------------------------------------------------------------------
# The pins and the listing


def pins_state(pin: str | None = ..., pins: dict | None = None) -> str:
    """'empty' (every pin None), 'filled' (every pin set, the 11 names exactly) or 'mixed'."""
    pin = K27_LISTING_SHA256 if pin is ... else pin
    pins = K27_TABLE_SHA256 if pins is None else pins
    if set(pins) != set(TABLE_NAMES):
        return "mixed"
    vals = [pin, *(pins[n] for n in TABLE_NAMES)]
    if all(v is None for v in vals):
        return "empty"
    if all(isinstance(v, str) and re.fullmatch(r"[0-9a-f]{64}", v) for v in vals):
        return "filled"
    return "mixed"


def check_pins_filled(pin: str | None = ..., pins: dict | None = None) -> None:
    if pins_state(pin, pins) != "filled":
        raise Refused("the pins K27_LISTING_SHA256 and K27_TABLE_SHA256 are not filled (k27_chip_query.py; the row's "
                      "second commit fills them from the committed TABLES_K27.sha256): refused (no counted step runs)")
    pins = K27_TABLE_SHA256 if pins is None else pins
    if len(set(pins.values())) != len(pins):
        raise Refused("a sha256 appears twice among the 11 table pins: refused")


def read_listing(path: Path | None = None, pin: str | None = ..., pins: dict | None = None,
                 check_tracked: bool = True) -> tuple[dict[str, str], str]:
    """(listing, its sha256), checked: pins filled, the file no symlink, tracked, at its pin, of the fixed form
    (outside_forms.parse_listing: `<sha256>  <name>` lines, sorted, LF, no name twice), names exactly the 11, each at
    its pin, no sha256 twice."""
    path = LISTING_FILE if path is None else Path(path)
    pin = K27_LISTING_SHA256 if pin is ... else pin
    pins = K27_TABLE_SHA256 if pins is None else pins
    check_pins_filled(pin, pins)
    if path.is_symlink() or not path.is_file():
        raise Refused(f"{path.name} is missing or a symlink: refused")
    if check_tracked and not ext.tracked(path):
        raise Refused(f"{path.name} is not tracked by git: refused")
    raw = path.read_bytes()
    sha = sha256_bytes(raw)
    if sha != pin:
        raise Refused(f"{path.name} hashes to {sha[:12]}, not the pinned K27_LISTING_SHA256 {pin[:12]}: refused")
    try:
        listing = ext.parse_listing(raw)
    except ext.TableRefused as e:
        raise Refused(f"{path.name}: {e}") from None
    if sorted(listing) != list(TABLE_NAMES):
        raise Refused(f"{path.name} does not list exactly the 11 tables: refused")
    if len(set(listing.values())) != len(listing):
        raise Refused(f"{path.name}: a sha256 appears twice among the 11 listing lines: refused")
    for n in TABLE_NAMES:
        if listing[n] != pins[n]:
            raise Refused(f"{n}: the listing says {listing[n][:12]}, the pin {str(pins[n])[:12]}: refused")
    return listing, sha


def listing_empty(path: Path | None = None, check_tracked: bool = True) -> None:
    """At the registration commit: the listing is tracked and empty, and every pin is None (else Refused)."""
    path = LISTING_FILE if path is None else Path(path)
    if pins_state() != "empty":
        raise Refused("the pins are already filled (or partly): P2 runs only at the registration commit, before the "
                      "second commit fills them: refused")
    if path.is_symlink() or not path.is_file():
        raise Refused(f"{path.name} is missing or a symlink (it is committed, empty, with the row): refused")
    if check_tracked and not ext.tracked(path):
        raise Refused(f"{path.name} is not tracked by git: refused")
    if path.read_bytes() != b"":
        raise Refused(f"{path.name} is not empty: its lines are filled by the second commit only: refused")


def listing_text(tables: dict[str, str]) -> str:
    """The listing's bytes as text: `<sha256>  <name>` per table, sorted by name, LF."""
    if sorted(tables) != list(TABLE_NAMES):
        raise Refused("not exactly the 11 tables: refused")
    if len(set(tables.values())) != len(tables):
        raise Refused("a sha256 appears twice among the 11 listing lines: refused")
    return "".join(f"{tables[n]}  {n}\n" for n in TABLE_NAMES)


def pins_text(tables: dict[str, str]) -> str:
    """The pin lines the second commit writes into this script, from the listing's text."""
    lt = listing_text(tables)
    body = "".join(f'    "{n}": "{tables[n]}",\n' for n in TABLE_NAMES)
    return (f'K27_LISTING_SHA256: str | None = "{sha256_bytes(lt.encode("ascii"))}"\n'
            f"K27_TABLE_SHA256: dict[str, str | None] = {{\n{body}}}\n")


# ---------------------------------------------------------------------------------------------------------------------
# The tables, `two` and their records


def folder(var: str, env: dict | None = None) -> Path:
    """The folder a variable names, through outside_forms.table_dir's checks unchanged (set, absolute, no symlink, a
    real directory, outside every checkout)."""
    env = os.environ if env is None else env
    v = env.get(var)
    if not v:
        raise Refused(f"the environment variable {var} is not set: refused (no table path is taken otherwise)")
    try:
        return ext.table_dir({ext.ENV: v})
    except ext.TableRefused as e:
        raise Refused(str(e).replace(ext.ENV, var)) from None


def tables_dir(env: dict | None = None) -> Path:
    try:
        return cq.tables_folder(env)
    except cq.Refused as e:
        raise Refused(str(e)) from None


def _plain_file(p: Path, what: str) -> Path:
    if p.is_symlink() or not p.is_file():
        raise Refused(f"{what} is missing or a symlink: refused")
    return p


def _json(p: Path, what: str) -> tuple[dict, str]:
    raw = _plain_file(p, what).read_bytes()
    try:
        return json.loads(raw), sha256_bytes(raw)
    except ValueError:
        raise Refused(f"{what} is not JSON: refused") from None


def _jsonable(v):
    return json.loads(json.dumps(v, default=float))


def fetch_manifest_sha256(root=None) -> str:
    """The sha256 of the fetch's MANIFEST.json (every table's manifest records the one it was built against)."""
    try:
        r = cf.root_of(root)
    except cf.Refused as e:
        raise Refused(str(e)) from None
    return cf.sha256_file(_plain_file(r / "MANIFEST.json", "the fetch's MANIFEST.json"))


def load_k27_table(name: str, genes: list[str], listing: dict[str, str], pins: dict | None = None,
                   d: Path | None = None, fetch_sha: str | None = None) -> tuple[np.ndarray, dict]:
    """A listed table in `genes`' order (which must be the file's), read once: its bytes at the listing's sha256 and
    the pin, its factors.json and manifest recording that sha256, the manifest of this candidate, query and seed, built
    against `fetch_sha` (the fetch manifest now) under chipq's readings; values finite and positive."""
    import pandas as pd

    pins = K27_TABLE_SHA256 if pins is None else pins
    if name not in TABLE_NAMES:
        raise Refused("not one of the 11 tables: refused")
    want = listing.get(name)
    if want is None or pins.get(name) != want:
        raise Refused(f"{name}: not listed, or the listing's sha256 is not the pinned one: refused")
    d = tables_dir() if d is None else Path(d)
    fd = d / "factors"
    if fd.is_symlink() or not fd.is_dir():
        raise Refused(f"{name}: the factors folder is missing or a symlink: refused")
    raw = _plain_file(fd / f"{name}.csv", f"{name}'s table").read_bytes()  # read once: hashed bytes are parsed
    sha = sha256_bytes(raw)
    if sha != want:
        raise Refused(f"{name}: the file hashes to {sha[:12]}, not the listed {want[:12]}: refused")
    fj, _ = _json(fd / "factors.json", f"{name}'s factors.json")
    if ((fj.get("files") or {}).get(name) or {}).get("sha256") != sha:
        raise Refused(f"{name}: factors/factors.json records another sha256: refused")
    man, man_sha = _json(d / f"{name}_manifest.json", f"{name}'s manifest")
    if man.get("name") != name or ((man.get("files") or {}).get(name) or {}).get("sha256") != sha:
        raise Refused(f"{name}: its manifest records another table or sha256: refused")
    if man.get("candidate") != "CHIPQ-S" or man.get("register_row") != "X018" or man.get("query") != QUERY:
        raise Refused(f"{name}: its manifest is not of this candidate's query (CHIPQ-S, X018): refused")
    if man.get("readings") != _jsonable(cq.readings()):
        raise Refused(f"{name}: its manifest was written under other readings than chipq's: refused")
    if fetch_sha is not None and man.get("fetch_manifest_sha256") != fetch_sha:
        raise Refused(f"{name}: its manifest records another fetch manifest than the one on disk: refused")
    seed = None
    if name != QUERY_TABLE:
        seed = next(s for s, n in PLACEBO_TABLE.items() if n == name)
        if man.get("seed") != seed or not isinstance(man.get("pool_index"), int) \
                or not isinstance(man.get("pool_size"), int):
            raise Refused(f"{name}: its manifest records another seed, or no pool index: refused")
    elif man.get("antigen") != QUERY:
        raise Refused(f"{name}: its manifest is not the query's: refused")
    try:
        df = pd.read_csv(io.BytesIO(raw), float_precision="round_trip")
    except Exception:  # noqa: BLE001
        raise Refused(f"{name}: not a readable CSV: refused") from None
    if list(df.columns) != ["gene", "factor"]:
        raise Refused(f"{name}: columns are not ['gene', 'factor']: refused")
    got = list(df["gene"].astype(str))
    if got != list(genes):
        raise Refused(f"{name}: genes differ from the barred factor's genes (order included): refused")
    try:
        T = df["factor"].to_numpy(dtype=np.float64)
    except (TypeError, ValueError):
        raise Refused(f"{name}: factor values are not numbers: refused") from None
    if not np.all(np.isfinite(T) & (T > 0)):
        raise Refused(f"{name}: factor values must be finite and positive: refused")
    rec = {"name": name, "sha256": sha, "manifest_sha256": man_sha, "seed": seed, "pool_index": man.get("pool_index"),
           "pool_size": man.get("pool_size"), "n_kept": man.get("n_kept"), "n_bound": man.get("n_bound"),
           "K": man.get("K"), "exponents": man.get("exponents"), "fetch_manifest_sha256": man.get(
               "fetch_manifest_sha256")}
    return T, rec


def load_k27_tables(genes: list[str], path: Path | None = None, pin: str | None = ..., pins: dict | None = None,
                    d: Path | None = None, fetch_sha: str | None = ..., check_tracked: bool = True) -> dict:
    """The listing checked against its pins, then all 11 tables loaded (any failure refuses the whole run); no two
    tables with the same bytes. {'listing_sha256', 'tables': {name: (T, record)}, 'fetch_manifest_sha256'}."""
    pins = K27_TABLE_SHA256 if pins is None else pins
    listing, lsha = read_listing(path, pin, pins, check_tracked)
    fetch_sha = fetch_manifest_sha256() if fetch_sha is ... else fetch_sha
    tables = {n: load_k27_table(n, genes, listing, pins, d, fetch_sha) for n in TABLE_NAMES}
    q = tables[QUERY_TABLE][0]
    if any(np.array_equal(tables[n][0], q) for n in PLACEBO_TABLE.values()):
        raise Refused("a placebo table equals the query's: refused")
    idx = [tables[n][1]["pool_index"] for n in PLACEBO_TABLE.values()]
    if len(set(idx)) != len(idx) or len({tables[n][1]["pool_size"] for n in PLACEBO_TABLE.values()}) != 1:
        raise Refused("the placebo tables' pool indices are not ten different draws from one pool: refused")
    return {"listing_sha256": lsha, "tables": tables, "fetch_manifest_sha256": fetch_sha}


def load_two(genes: list[str], env: dict | None = None) -> tuple[np.ndarray, dict]:
    """`two`: K22's table `all` through K22's listing and pins (outside_forms.load_listed, K23's code path), from the
    folder VEC_T3_EXT_TABLES_TWO names; its record must be the one K22 pinned and K23 read."""
    d = folder(ENV_TWO, env)
    if ext.K22_LISTING_SHA256 != k23.K23_LISTING_SHA256 or ext.K22_TABLE_SHA256.get("all") != k23.K23_TABLE_SHA256:
        raise Refused("outside_forms.py's pins are not the listing and `all` sha256 K22 pinned and K23 read: refused")
    try:
        loaded = ext.load_listed(list(genes), env={ext.ENV: str(d)})
    except (ext.TableRefused, ext.FormRefused) as e:
        raise Refused(str(e).replace(ext.ENV, ENV_TWO)) from None
    if "all" not in loaded["tables"]:
        raise Refused("`two` (K22's `all`) is not listed: refused")
    rec = ext.ext_record(loaded)
    if rec["listing_sha256"] != k23.K23_LISTING_SHA256 or rec["tables"].get("all") != k23.K23_TABLE_SHA256:
        raise Refused("`two` is not the table K22 pinned and K23 read: refused")
    return loaded["tables"]["all"][0], {"record": rec, "loaded": loaded}


def exponents_of(tables: dict) -> dict[str, dict[str, float]]:
    """e per table and size key (chipq.exponent, population SD over the genes in file order), each checked: SD(ln
    C ** e) equals S within REL_TOL and e equals the one the table's manifest records (else Refused)."""
    out = {}
    for n, (T, rec) in tables.items():
        out[n] = {}
        for s, S in SIZES.items():
            e = cq.exponent(T, S)
            if e is None:
                raise Refused(f"{n}: SD(ln C) is 0: no exponent: refused")
            sd = float(np.std(np.log(np.asarray(T, dtype=np.float64) ** e)))
            if not math.isclose(sd, S, rel_tol=REL_TOL, abs_tol=0.0):
                raise Refused(f"{n} at {SIZE_SHOWN[s]}: SD(ln C ** e) is not S within a relative {REL_TOL}: refused")
            got = (rec.get("exponents") or {}).get(SIZE_OF_KEY[s])
            if not isinstance(got, (int, float)) or not math.isclose(float(got), e, rel_tol=REL_TOL, abs_tol=0.0):
                raise Refused(f"{n} at {SIZE_SHOWN[s]}: e differs from the one its manifest records: refused")
            out[n][s] = float(e)
    return out


def ext_record(k27_loaded: dict, two_rec: dict, kept_sha: str | None) -> dict:
    """What every task file records of the tables (the VOID check compares it across files)."""
    return {"listing_sha256": k27_loaded["listing_sha256"],
            "tables": {n: r["sha256"] for n, (_, r) in k27_loaded["tables"].items()},
            "table_manifests": {n: r["manifest_sha256"] for n, (_, r) in k27_loaded["tables"].items()},
            "fetch_manifest_sha256": k27_loaded["fetch_manifest_sha256"], "kept_sha256": kept_sha, "two": two_rec}


def pinned_ext() -> dict:
    return {"listing_sha256": K27_LISTING_SHA256, "tables": dict(K27_TABLE_SHA256)}


# ---------------------------------------------------------------------------------------------------------------------
# P2 (database only, uncounted)


def p2_compute(genes: list[str], root=None, kept=None, man: dict | None = None) -> dict:
    """P2's bars on the database alone, in memory: Gata4's kept count and |B| at P1's bars, the pool size (at least
    POOL_MIN_SIZE), genes changed by Gata4's table (at least 10) and by each seed's placebo table (at least 1), the
    scaled SD at each size (refused off S), no placebo table equal to the query's (refused), and whether the analysis
    list holds an antigen for the barred factor's gene. Returns counts only (no name, no gene, no database score)."""
    root = cf.root_of(root)
    man = cf.load_manifest(root) if man is None else man
    rq = cq.build(QUERY, genes, root, kept, man)
    out = {"query_kept": int(rq["n_kept"]), "query_bound": int(rq["n_bound"]), "query_K": int(rq["K"]),
           "query_genes_changed": int((rq["C"] != 1.0).sum()), "pool_size": None, "placebos": {},
           "barred_gene_antigen": int(BARRED_GENE in cf._lists(root, man)[0]), "why": []}
    if rq["n_kept"] < cq.P1_MIN_EXPERIMENTS or rq["n_bound"] < cq.P1_MIN_BOUND:
        out["why"].append(f"the practice query's kept experiments or bound set are under P1's bars "
                          f"({cq.P1_MIN_EXPERIMENTS} and {cq.P1_MIN_BOUND})")
    if out["query_genes_changed"] < P2_MIN_GENES_QUERY:
        out["why"].append(f"the practice query's table differs from 1 at fewer than {P2_MIN_GENES_QUERY} genes")
    pool = cq.placebo_pool(QUERY, rq["K"], genes, root, man)
    out["pool_size"] = len(pool)
    if len(pool) < max(cq.POOL_MIN_SIZE, len(SEEDS)):
        out["why"].append(f"the placebo pool holds fewer than {cq.POOL_MIN_SIZE} antigens")
        return out | {"passes": False}
    tabs = {QUERY_TABLE: rq["C"]}
    for s, a, i in cq.placebo_draws(pool):
        rp = cq.build(a, genes, root, None, man)
        out["placebos"][str(s)] = {"pool_index": int(i), "kept": int(rp["n_kept"]), "bound": int(rp["n_bound"]),
                                   "K": int(rp["K"]), "genes_changed": int((rp["C"] != 1.0).sum())}
        if out["placebos"][str(s)]["genes_changed"] < P2_MIN_GENES_PLACEBO:
            out["why"].append(f"seed {s}'s placebo table differs from 1 at no gene")
        tabs[PLACEBO_TABLE[int(s)]] = rp["C"]
    if out["why"]:
        return out | {"passes": False}
    for n in PLACEBO_TABLE.values():
        if np.array_equal(tabs[n], tabs[QUERY_TABLE]):
            raise Refused("a placebo's finished table equals the query's: refused (nothing counted, before P3)")
    out["exponents"] = exponents_of({n: (T, {"exponents": {SIZE_OF_KEY[s]: cq.exponent(T, S)
                                                           for s, S in SIZES.items()}}) for n, T in tabs.items()})
    return out | {"passes": True}


def p2_write_tables(genes: list[str], d: Path, root=None, kept=None) -> dict[str, str]:
    """chipq.tables into the K27 folder; then the 11 written files hashed again: no sha256 twice and no placebo's bytes
    equal to the query's (else Refused). {name: sha256}."""
    written = cq.tables(genes, d, None, root, kept)
    got = {}
    for fname, sha in written:
        name = fname[:-len(".csv")]
        raw = _plain_file(Path(d) / "factors" / fname, f"{name}'s table").read_bytes()
        if sha256_bytes(raw) != sha:
            raise Refused(f"{name}: the written file does not hash to the sha256 chipq returned: refused")
        got[name] = sha
    if sorted(got) != list(TABLE_NAMES):
        raise Refused("chipq.tables did not write exactly the 11 tables: refused")
    q = (Path(d) / "factors" / f"{QUERY_TABLE}.csv").read_bytes()
    if any((Path(d) / "factors" / f"{n}.csv").read_bytes() == q for n in PLACEBO_TABLE.values()):
        raise Refused("a placebo's finished table equals the query's byte for byte: refused (nothing counted)")
    listing_text(got)  # refuses a sha256 twice among the 11 lines
    return got


def p2_lines(res: dict) -> list[str]:
    lines = [f"P2: practice query kept experiments {res['query_kept']}, bound set {res['query_bound']}, genes changed "
             f"{res['query_genes_changed']}; placebo pool {res['pool_size']}; the analysis list holds an antigen for "
             f"the barred factor's gene: {res['barred_gene_antigen']}"]
    for s, v in sorted(res["placebos"].items(), key=lambda kv: int(kv[0])):
        lines.append(f"  seed {s}: pool index {v['pool_index']}, kept {v['kept']}, bound set {v['bound']}, genes "
                     f"changed {v['genes_changed']}")
    lines.append("P2: " + ("passes" if res["passes"] else "fails: " + "; ".join(res["why"])))
    return lines


def cmd_p2(a) -> None:
    """P2 at the registration commit: refused when the pins are filled; R/p2.json; exit 0, 3 (a bar not met: the STOP
    reading follows) or 2."""
    rid = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse(DIRTY)
    out = Path(a.out)
    path = out / "p2.json"
    if path.exists():
        refuse(f"{path} exists: P2 is made once")
    try:
        listing_empty()
        d = tables_dir()
    except Refused as e:
        refuse(str(e))
    try:
        with contextlib.redirect_stdout(sys.stderr):  # a loader's lines go to stderr
            kept = cq.load_kept()
            genes = cq.barred_genes(Path(a.data), Path(a.run_dir), Path(a.decision))
            res = p2_compute(genes, None, kept)
            fsha = fetch_manifest_sha256()
            tables = p2_write_tables(genes, d, None, kept) if res["passes"] else {}
    except (Refused, cq.Refused, cf.Refused) as e:
        refuse(str(e))
    except Exception as e:  # noqa: BLE001  (another module's message could name a gene or antigen: its type only)
        refuse(f"{type(e).__name__} in P2")
    body = {**res, "row": rid, "kept_sha256": kept.sha256, "fetch_manifest_sha256": fsha, "readings": cq.readings(),
            "tables": tables, "listing": listing_text(tables) if tables else None,
            "listing_sha256": sha256_bytes(listing_text(tables).encode("ascii")) if tables else None,
            "provenance": st}
    out.mkdir(parents=True, exist_ok=True)
    k18.write_task_file(path, body)
    print("\n".join(p2_lines(res)), flush=True)
    if not res["passes"]:
        sys.exit(EXIT_BAR)
    for n in TABLE_NAMES:
        print(f"{tables[n]}  {n}", flush=True)
    print(f"listing sha256 {body['listing_sha256']}", flush=True)


def load_p2(runs: Path) -> tuple[dict, str]:
    p = Path(runs) / "p2.json"
    if not p.is_file():
        raise Refused(f"{p} is missing: P2 runs first, at the registration commit")
    return _json(p, "p2.json")


def cmd_pins(a) -> None:
    """The second commit's lines from R/p2.json, with every table hashed again from the K27 folder: the listing's text,
    then the pin lines. sha256 and table names only."""
    try:
        p2, _ = load_p2(Path(a.runs))
        if p2.get("passes") is not True or sorted(p2.get("tables") or {}) != list(TABLE_NAMES):
            raise Refused("P2 did not pass: there is no listing to commit")
        d = tables_dir()
        for n in TABLE_NAMES:
            raw = _plain_file(d / "factors" / f"{n}.csv", f"{n}'s table").read_bytes()
            if sha256_bytes(raw) != p2["tables"][n]:
                raise Refused(f"{n}: the table on disk is not the one P2 wrote: refused")
        lt = listing_text(p2["tables"])
        if lt != p2.get("listing") or sha256_bytes(lt.encode("ascii")) != p2.get("listing_sha256"):
            raise Refused("p2.json's listing is not its tables': refused")
    except Refused as e:
        refuse(str(e))
    sys.stdout.write("# experiments/t3_ext/TABLES_K27.sha256 (exactly these 11 lines)\n" + lt
                     + "# experiments/t3_decomp/k27_chip_query.py (the pin lines)\n" + pins_text(p2["tables"]))


# ---------------------------------------------------------------------------------------------------------------------
# P3, the pre-gate (counted)


def _world_vectors(v, genes, wd) -> np.ndarray:
    return k23._world_vectors(np.asarray(v, dtype=np.float64), genes, wd)


def _with_max_abs_ln(f: dict, v) -> dict:
    return f | {"max_abs_ln": float(np.max(np.abs(np.log(np.asarray(v, dtype=np.float64)))))}


def pregate_k27(phi_m, two, tabs: dict[str, np.ndarray], e: dict[str, dict[str, float]], bases_by_world: dict, den,
                genes: list[str] | None = None, pool_index: dict[str, int] | None = None) -> dict:
    """P3's figures, a pure function of phi_M, `two`, the 11 raw tables (file order), their exponents, each world's
    basis and WT E8.75's mean within-cell proportions: both correlations (over genes != 1 in either), the drop, and per
    size the scaled table's and the product's figures, the twin's a and k, each world's chip_perm figures and K24's
    bounds (recorded, not applied); the placebos' e by pool index."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    two = np.asarray(two, dtype=np.float64)
    C = np.asarray(tabs[QUERY_TABLE], dtype=np.float64)
    x_e2 = ext.product(phi_m, 1.0, two, 2.0)
    sd_phi = float(np.std(np.log(phi_m)))
    sd_x = float(np.std(np.log(x_e2)))
    den_phi = float(den @ phi_m / den.sum())
    den_x = float(den @ x_e2 / den.sum())
    c1 = ext.corr_between(C, phi_m)
    c2 = ext.corr_between(C, two ** 2.0)
    why = []
    for label, c in (("corr(ln C_Gata4, ln phi_M)", c1), ("corr(ln C_Gata4, ln two ** 2)", c2)):
        if c is None or abs(c) >= MAX_ABS_CORR:
            why.append(f"abs {label} {'undefined' if c is None else f'{abs(c):.3f}'} (>= {MAX_ABS_CORR} or undefined)")
    out: dict = {"n_genes": int(phi_m.size), "s_star": S_STAR, "sizes": dict(SIZES),
                 "corr": {"ln_phi_m": c1, "ln_two_sq": c2}, "dropped": "; ".join(why) or None,
                 "phi_m": ext.figures(phi_m, phi_m, den), "x_e2": _with_max_abs_ln(ext.figures(phi_m, x_e2, den), x_e2),
                 "e": {n: dict(v) for n, v in e.items()}, "arms": {}, "k24_bounds": {}, "placebos": {},
                 "constants": {"max_abs_corr": MAX_ABS_CORR, "rel_tol": REL_TOL, "perm_seed": ext.PERM_SEED,
                               "perm_bins": ext.PERM_BINS, "worlds": sorted(bases_by_world),
                               "k24_bounds_recorded_not_applied": K24_BOUNDS}}
    for s in SIZE_KEYS:
        eq = e[QUERY_TABLE][s]
        Ce = C ** eq
        P = ext.product(x_e2, 1.0, Ce, 1.0)
        sd_prod = float(np.std(np.log(P)))
        a = sd_prod / sd_phi
        k = float((den @ P) / (den @ (phi_m ** a)))
        twin = k * phi_m ** a
        pf = _with_max_abs_ln(ext.figures(phi_m, P, den), P)
        pf["genes_changed_vs_x_e2"] = int((P != x_e2).sum())
        controls = {}
        for w in sorted(bases_by_world):
            wd = bases_by_world[w]
            phi_w = _world_vectors(phi_m, genes, wd)
            x_w = ext.product(phi_w, 1.0, _world_vectors(two, genes, wd), 2.0)
            Ce_w = _world_vectors(Ce, genes, wd)
            den_w = _world_vectors(den, genes, wd)
            controls[w] = ext.figures(phi_w, ext.product(x_w, 1.0, ext.permuted(Ce_w, wd["basis"]), 1.0), den_w)
        shift = pf["den_total"] / den_phi - 1.0
        out["arms"][s] = {"size": SIZES[s], "e": eq, "table": _with_max_abs_ln(ext.table_figures(phi_m, Ce, den), Ce),
                          "product": pf, "a": float(a), "k": float(k), "twin": ext.figures(phi_m, twin, den),
                          "controls": controls, "sd_ln_phi_m": sd_phi, "sd_ln_x_e2": sd_x, "sd_ln_product": sd_prod,
                          "den_shift_vs_phi_m": float(shift),
                          "den_shift_vs_x_e2": float(pf["den_total"] / den_x - 1.0)}
        mx = float(np.max(np.abs(np.log(Ce))))
        out["k24_bounds"][s] = {
            "max_abs_ln_scaled": {"value": mx, "bound": K24_BOUNDS["max_abs_ln_scaled"],
                                  "beyond": mx > K24_BOUNDS["max_abs_ln_scaled"]},
            "sd_ratio": {"value": sd_prod / sd_phi, "bound": K24_BOUNDS["max_sd_ratio"],
                         "beyond": sd_prod > K24_BOUNDS["max_sd_ratio"] * sd_phi},
            "den_shift": {"value": float(shift), "bound": K24_BOUNDS["max_den_shift"],
                          "beyond": abs(shift) > K24_BOUNDS["max_den_shift"]},
            "genes_at_clip": {"value": pf["genes_at_clip"], "bound": K24_BOUNDS["max_at_clip"],
                              "beyond": pf["genes_at_clip"] > K24_BOUNDS["max_at_clip"]},
            "e": {"value": eq, "bound": K24_BOUNDS["e_range"],
                  "beyond": not (K24_BOUNDS["e_range"][0] <= eq <= K24_BOUNDS["e_range"][1])},
            "applied": False}
    for sd, n in PLACEBO_TABLE.items():
        idx = (pool_index or {}).get(n)
        out["placebos"][str(idx if idx is not None else f"seed{sd}")] = {"e": dict(e[n])}
    return out


def constants_of(pre: dict) -> dict:
    """What every arm needs from the pre-gate, at full precision, applied unchanged in both worlds: e per table and
    size, and the twins' a and k per size."""
    return {"e": {n: dict(v) for n, v in pre["e"].items()},
            "twins": {s: {"a": pre["arms"][s]["a"], "k": pre["arms"][s]["k"]} for s in SIZE_KEYS}}


def dropped(pre: dict) -> bool:
    return pre.get("dropped") is not None


def _close(a, b, rel: float = REL_TOL) -> bool:
    if isinstance(a, (float, int)) and isinstance(b, (float, int)) and not isinstance(a, bool) \
            and not isinstance(b, bool) and (isinstance(a, float) or isinstance(b, float)):
        return math.isclose(float(a), float(b), rel_tol=rel, abs_tol=0.0) or a == b
    return a == b


def _same(got, want) -> bool:
    if isinstance(want, dict):
        return isinstance(got, dict) and set(got) == set(want) and all(_same(got[k], want[k]) for k in want)
    return _close(got, want)


def check_two(phi_m, two, two_loaded: dict, bases_by_world: dict, den, genes, k23_pre: dict) -> dict:
    """`two` at e = 2 through K23's own pre-gate (k23.pregate_exponents): S* and the table, product, controls, twin, a
    and k equal K23's committed pregate.json figures to a relative 1e-12, else Refused (row K24's check)."""
    p = k23.pregate_exponents(phi_m, two_loaded, bases_by_world, den, genes)
    got, want = (p.get("exponents") or {}).get("2"), ((k23_pre or {}).get("exponents") or {}).get("2")
    if got is None or want is None:
        raise Refused("`two` at e = 2: K23's pre-gate figures are missing: refused")
    sd = float(np.std(np.log(np.asarray(two, dtype=np.float64) ** 2.0)))
    if not _close(sd, S_STAR) or not _close(float(want["table"]["sd_ln"]), S_STAR):
        raise Refused(f"S* recomputed from `two` is not {S_STAR!r} (relative {REL_TOL}): refused")
    bad = [k for k in ("table", "product", "controls", "size", "a", "k", "den_shift", "sd_ln_phi_m", "sd_ln_product")
           if not _same(got.get(k), want.get(k))]
    if bad:
        raise Refused(f"`two` at e = 2 differs from K23's committed pregate.json in {bad}: refused")
    return got


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def fold_factors(bs, ref, per: dict, phi: np.ndarray, group: np.ndarray, basis: np.ndarray, phi_barred: np.ndarray,
                 spec: dict) -> dict:
    """K18's fold_factors (the original, called unchanged) plus `two` and every scaled table in the reference's gene
    order, and the query's scaled tables permuted within deciles of the `basis` this function receives (one per world:
    perm of the scaled vector)."""
    if not _K27.get("T") or _K27.get("two") is None:
        raise RuntimeError("the K27 tables are not loaded (k27_chip_query.py loads them before any arm)")
    F = k22._K18_ORIGINAL["fold_factors"](bs, ref, per, phi, group, basis, phi_barred, spec)
    order = lambda v: k18.in_gene_order(v, _K27["genes"], list(ref.genes))  # noqa: E731
    F["two"] = order(_K27["two"])
    F["scaled"], F["q_perm"] = {}, {}
    e = (_K27.get("consts") or {}).get("e")
    if e is not None:
        for n, T in _K27["T"].items():
            F["scaled"][n] = {s: order(T) ** e[n][s] for s in SIZE_KEYS}
        F["q_perm"] = {s: ext.permuted(F["scaled"][QUERY_TABLE][s], basis) for s in SIZE_KEYS}
    return F


def x_e2_factor(F: dict) -> np.ndarray:
    return ext.product(F["phi_barred"], 1.0, F["two"], 2.0)


def arm_factor(a: str, F: dict, seed: int | None = None) -> np.ndarray:
    """The factor an arm passes in place of phi, in the reference's gene order; chip_rand takes seed's placebo."""
    if a == BASE:
        return F["phi_barred"]
    if a == REPRO:
        return x_e2_factor(F)
    if a not in TRY_ARMS:
        raise RuntimeError(f"unknown arm {a!r}")
    kind, s = kind_of_arm(a), size_of_arm(a)
    if kind == "chip_size":
        tw = ((_K27.get("consts") or {}).get("twins") or {}).get(s)
        if tw is None:
            raise RuntimeError(f"the twin's constants at {s} are not loaded from pregate.json")
        return tw["k"] * F["phi_barred"] ** tw["a"]
    if not F.get("scaled"):
        raise RuntimeError("the scaled tables are not made (the exponents are not loaded from pregate.json)")
    if kind == "chip_q":
        t = F["scaled"][QUERY_TABLE][s]
    elif kind == "chip_perm":
        t = F["q_perm"][s]
    else:
        if seed is None or int(seed) not in PLACEBO_TABLE:
            raise RuntimeError("chip_rand needs the seed whose placebo it takes")
        t = F["scaled"][PLACEBO_TABLE[int(seed)]][s]
    return ext.product(x_e2_factor(F), 1.0, t, 1.0)


def build_arm(bs, ref, a: str, n: int, seed: int, F: dict):
    """One arm's drawn cells, as K22's to K24's build_arm: every arm draws from stable_rng(seed, STREAM) through
    bs.construct with the barred spec's keywords (central composition, w = 1, no depth or stage step); only the factor
    passed differs. barred_form passes them unchanged."""
    if a not in ARMS:
        raise RuntimeError(f"unknown arm {a!r}")
    rng = bs.t3d.stable_rng(seed, STREAM)
    kw, _ = bs.prediction_inputs(ref, F["barred_spec"], {}, F["phi_barred"])
    if kw["ratio"] is not None or kw["alphas"] is not None or kw["w"] != 1.0:
        raise RuntimeError("the barred s1 is the Mab21l2 factor at exponent 1 with no depth ratio or stage mixing")
    if a != BASE:
        kw["phi"] = arm_factor(a, F, seed)
    return bs.construct(ref, n, rng, **kw)


_K18_ORIGINAL = dict(k22._K18_ORIGINAL)  # k18's own ARMS, build_arm and fold_factors


@contextmanager
def k18_arms(arms: tuple[str, ...]):
    """For this process only: k18's ARMS, build_arm and fold_factors replaced by this try's, restored on exit and on an
    exception."""
    if not set(arms) <= set(ARMS) or BASE not in arms:
        raise RuntimeError(f"arms {arms} are not a set of this try's arms with the base")
    k18.ARMS, k18.build_arm, k18.fold_factors = tuple(arms), build_arm, fold_factors
    k18._FOLD.clear()
    try:
        yield
    finally:
        k18.ARMS = _K18_ORIGINAL["ARMS"]
        k18.build_arm = _K18_ORIGINAL["build_arm"]
        k18.fold_factors = _K18_ORIGINAL["fold_factors"]
        k18._FOLD.clear()


def candidates(arms: tuple[str, ...]) -> list:
    def one(a):
        return harness.Candidate(a, frozenset({"T3"}), (), lambda ctx, a=a: k18._build_plain(ctx, a),
                                 params={"form": a, "stream": STREAM, "composition": k8.spec()["composition"],
                                         "factor": ("the barred s1's" if a == BASE else f"k27 ({a})"),
                                         "perm_seed": ext.PERM_SEED if kind_of_arm(a) == "chip_perm" else None},
                                 requires=lambda t: t.key in k18.KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in arms)]


def try_count(failed_guards=()) -> int:
    n = N_TRIES + len([x for x in failed_guards if x in GUARDS])
    if n > MAX_TRIES:
        raise RuntimeError(f"{n} tries: never more than {MAX_TRIES}")
    return n


# ---------------------------------------------------------------------------------------------------------------------
# The ledger row, the preconditions and the guard references


def row_present(path: Path | None = None) -> str:
    """The id of the TRIES.md row that names this script; refuses (exit 2) when there is none or more than one."""
    path = TRIES_MD if path is None else path
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.startswith("| ") and ROW_SCRIPT in ln]
    ids = sorted({ln[2:].split(" ", 1)[0] for ln in rows})
    if len(ids) != 1:
        refuse(f"{'no' if not ids else 'more than one'} row of {Path(path).name} names {ROW_SCRIPT} ({ids}); the row "
               "is committed before any run")
    return ids[0]


DIRTY = ("the tree is not clean (git status --porcelain --untracked-files=no prints a tracked file, or git failed): "
         "not started, nothing written. Row K27: commit or restore first; INCOMPLETE, not VOID")


def is_ancestor(older: str, newer: str) -> bool:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["GIT_OPTIONAL_LOCKS"] = "0"
    r = subprocess.run(["git", "-C", str(REPO), "merge-base", "--is-ancestor", str(older), str(newer)],
                       capture_output=True, env=env, check=False)
    return r.returncode == 0


def p2_matches(p2: dict, st: dict, kept_sha: str) -> None:
    """P2's record agrees with the pins and the kept list, and was made at an earlier commit (the registration commit,
    an ancestor of this one): else Refused."""
    if p2.get("passes") is not True:
        raise Refused("P2 did not pass: no counted step runs")
    if p2.get("tables") != dict(K27_TABLE_SHA256) or p2.get("listing_sha256") != K27_LISTING_SHA256:
        raise Refused("p2.json's tables or listing are not the pinned ones: refused")
    if p2.get("kept_sha256") != kept_sha:
        raise Refused("the kept list is not the one P2 read (chipq_kept.yaml changed after P2): refused")
    pc = (p2.get("provenance") or {}).get("commit")
    if not pc or pc == st["commit"] or not is_ancestor(pc, st["commit"]):
        raise Refused("p2.json was not made at an earlier commit of this line (the registration commit): the pins "
                      "are filled at the wrong commit: refused")
    if (p2.get("provenance") or {}).get("machine") != st["machine"] or (p2.get("provenance") or {}).get(
            "tree_clean") is not True:
        raise Refused("p2.json is from another machine or a dirty tree: refused")


def preconditions(runs: Path | None = None) -> tuple[str, dict]:
    """Before anything loads in a counted step: the row, a clean tree, the pins filled, the listing tracked and at its
    pins, both folder variables valid, the frozen files, K23's committed pregate.json and reading.json tracked, the kept
    list committed and P2's (when the run folder is given)."""
    rid = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse(DIRTY)
    try:
        read_listing()
        tables_dir()
        folder(ENV_TWO)
        kept = cq.load_kept()
        if runs is not None:
            p2, _ = load_p2(runs)
            p2_matches(p2, st, kept.sha256)
    except (Refused, cq.Refused) as e:
        refuse(str(e))
    k22.frozen_ok()
    for f in ("pregate.json", "reading.json"):
        p = K23_RESULTS / f
        if not p.is_file() or not ext.tracked(p):
            refuse(f"K23's committed {f} is missing or not tracked by git: refused")
    return rid, st


def k23_committed(k23_dir: Path | None = None) -> tuple[dict, dict]:
    k23_dir = K23_RESULTS if k23_dir is None else Path(k23_dir)
    out = []
    for f in ("pregate.json", "reading.json"):
        p = Path(k23_dir) / f
        if not p.is_file():
            raise Refused(f"K23's committed {f} is missing: refused")
        out.append(json.loads(p.read_text()))
    return out[0], out[1]


def k23_runs_default(k23_read: dict) -> Path:
    """K23's run folder in this checkout, as run_k23.sh named it: runs/k23_<first 7 of the reading's commit>."""
    return REPO / "runs" / f"k23_{str(k23_read.get('commit'))[:7]}"


def k23_seed_rows(k23_read: dict, runs: Path) -> dict:
    """K23's per-seed x_e2 rows per cell, from K23's task files at the sha256 K23's reading records (else Refused)."""
    want = k23_read.get("files_sha256") or {}
    out: dict = {}
    for w in WORLDS:
        for n in TASKS:
            rows = []
            for s in SEEDS:
                f = f"plain_{n}/seed_{s}.json" if w == "plain" else f"dm/{n}_seed_{s}.json"
                p = Path(runs) / f
                if not p.is_file() or want.get(f) is None:
                    raise Refused(f"K23's task file {f} is missing (K23's run folder {Path(runs).name}): refused")
                raw = p.read_bytes()
                if sha256_bytes(raw) != want[f]:
                    raise Refused(f"K23's task file {f} is not at the sha256 K23's reading records: refused")
                r = json.loads(raw)
                row = r["rows"][0] if w == "plain" else r
                rows.append({"seed": s, "scores": {REPRO: row["scores"][REPRO]}, "t3": {REPRO: row["t3"][REPRO]}})
            out[(w, n)] = rows
    return out


def guard_refs(k23_dir: Path | None = None, refs: Path = K18_TASKS) -> dict:
    """K18's committed barred_form rows (at their pins), K23's committed x_e2 cells and its per-seed x_e2 rows from
    K23's run folder in this checkout (no option overrides it)."""
    _, k23_read = k23_committed(k23_dir)
    saved = dict(k21.K18_SHA256)
    k21.K18_SHA256.clear()
    k21.K18_SHA256.update(K18_SHA256)
    try:
        k18r, why = k21.k18_rows(refs)
    finally:
        k21.K18_SHA256.clear()
        k21.K18_SHA256.update(saved)
    if why:
        raise Refused("; ".join(why) + ": refused")
    cells = {}
    for w in WORLDS:
        for n in TASKS:
            c = ((((k23_read.get("worlds") or {}).get(w) or {}).get(n) or {}).get("arms") or {}).get(REPRO)
            if not c:
                raise Refused(f"K23's committed reading.json holds no x_e2 in {w} {n}: refused")
            cells[(w, n)] = {"mean": float(c["mean"]), "severity_gated": int(c["severity_gated"])}
    return {"k18": k18r, "k23_cells": cells, "k23_rows": k23_seed_rows(k23_read, k23_runs_default(k23_read))}


# ---------------------------------------------------------------------------------------------------------------------
# Loading for a counted step


def load_k27(genes: list[str], kept_sha: str | None = None) -> dict:
    """The 11 tables and `two` against the barred factor's file genes (order included), held in _K27 for fold_factors.
    A refusal exits 2."""
    try:
        loaded = load_k27_tables(list(genes))
        two, two_info = load_two(list(genes))
        e = exponents_of(loaded["tables"])
    except Refused as e_:
        refuse(str(e_))
    _K27.clear()
    _K27.update(genes=list(genes), T={n: T for n, (T, _) in loaded["tables"].items()}, two=two, loaded=loaded,
                two_loaded=two_info["loaded"], e=e, record=ext_record(loaded, two_info["record"], kept_sha),
                table_records={n: r for n, (_, r) in loaded["tables"].items()})
    return loaded


def load_pregate(out: Path, record: dict) -> dict:
    """The run folder's pregate.json: its tables must be the loaded ones and its e the loaded tables'; fills _K27's
    constants. Refuses (exit 2) when it is missing, disagrees, or drops the arm."""
    p = Path(out) / "pregate.json"
    if not p.is_file():
        refuse(f"{p} is missing: the pre-gate runs first (`pregate --out {out}`)")
    body = json.loads(p.read_text())
    if body.get("ext") != _jsonable(record):
        refuse("pregate.json records other tables than the ones loaded: refused")
    pre = body["pregate"]
    if dropped(pre):
        refuse("the pre-gate dropped the arm and its controls: nothing runs")
    consts = constants_of(pre)
    if body.get("constants") != consts or not _same(consts["e"], _K27.get("e") or {}):
        refuse("pregate.json's constants are not its pre-gate's, or its e not the loaded tables': refused")
    _K27["consts"] = consts
    return pre


def _pregate_record(pre: dict) -> dict:
    return {"arms": list(ARMS), "ext": _jsonable(_K27["record"]), "constants": constants_of(pre),
            "sizes": dict(SIZES)}


def _world_figures(F: dict, arms: tuple[str, ...], seeds) -> dict:
    phi = F["phi_barred"]
    out = {}
    for a in arms:
        if a == BASE:
            continue
        if kind_of_arm(a) == "chip_rand":
            out[a] = {str(s): ext.figures(phi, arm_factor(a, F, s)) for s in seeds}
        else:
            out[a] = ext.figures(phi, arm_factor(a, F))
    return out


class _Mismatch(RuntimeError):
    pass


def check_world(pre: dict, world: str, F: dict) -> None:
    """The run's own chip_perm figures equal pregate.json's in this world (else no task file: INCOMPLETE)."""
    for s in SIZE_KEYS:
        want = ((pre["arms"].get(s) or {}).get("controls") or {}).get(world)
        got = ext.figures(F["phi_barred"], arm_factor(arm("chip_perm", s), F))
        if want is None or any(not _close(got[k], want.get(k)) for k in got):
            raise _Mismatch(f"the {world} world's chip_perm at {SIZE_SHOWN[s]} is not the pre-gate's: no task file "
                            "written (INCOMPLETE)")


GUARD_FILE = "guard.json"
GUARD_DIR = "guard"


def guard_passed(out: Path, st: dict) -> None:
    p = Path(out) / GUARD_FILE
    if not p.is_file():
        refuse(f"{p} is missing: the guard part runs and passes first (`guard --runs {out}`)")
    g = json.loads(p.read_text())
    prov = g.get("provenance") or {}
    if g.get("passes") is not True:
        refuse("the guard part did not pass: the run is VOID and no scoring part runs")
    if prov.get("commit") != st["commit"] or prov.get("machine") != st["machine"]:
        refuse("guard.json is from another commit or machine: refused")


def progress_log(part: str, world: str, task: str):
    lead = "K27 guard" if part == "guard" else "K27"

    def log(line: str) -> None:
        m = re.match(r".* seed (\d+):", str(line))
        print(f"{lead} {world} {task} seed {m.group(1) if m else '?'}: scored", flush=True)

    return log


def run_plain(a, st: dict, out: Path, part: str, kept_sha: str) -> None:
    root = out / GUARD_DIR if part == "guard" else out
    folder_ = root / f"plain_{a.task}"
    folder_.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder_, st)
    rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
    load_k27(k18._BARRED["genes"], kept_sha)
    pre = load_pregate(out, _K27["record"])
    run_arms = GUARDS if part == "guard" else ARMS
    lead = "K27 guard" if part == "guard" else "K27"
    with k18_arms(run_arms):
        for s in k21.parse_seeds(a.seeds):
            path = folder_ / f"seed_{s}.json"
            if path.exists():
                print(f"{lead} plain {a.task} seed {s}: present, skipped", flush=True)
                continue
            res = harness.run_task(TASKS[a.task], candidates(run_arms), seeds=[s], data_dir=Path(a.data),
                                   log=progress_log(part, "plain", a.task))
            got = res["provenance"]
            if got.get("tree_clean") is not True or got.get("commit") != st["commit"]:
                raise k22._Unclean(f"the harness's stamp says commit {got.get('commit')}, clean "
                                   f"{got.get('tree_clean')}: no task file written (INCOMPLETE)")
            k22._stamp_ok(st)
            fold = k18._FOLD.get("F")
            if fold is None:
                raise _Mismatch("the fold's factors were not made: no task file written")
            if part != "guard":
                check_world(pre, "plain", fold)
            res["provenance"]["machine"] = st["machine"]
            res |= {"world": "plain", "part": part, "task_source": SCRIPT, "seed": s, "pregate": _pregate_record(pre),
                    "factor_record": k7._PREP["fold"]["record"], "barred_factor": rec,
                    "world_figures": _world_figures(fold, run_arms, [s])}
            prov = k18.write_task_file(path, res)
            print(f"wrote {path}: commit {prov['commit']}, clean {prov['tree_clean']}, machine {prov['machine']}",
                  flush=True)


def run_dm(a, st: dict, out: Path, part: str, kept_sha: str) -> None:
    """The depth-matched world, set up exactly as k24.run_dm sets it up, then each (replicate, seed) row through
    k18.dm_row with this part's arms; rows already written at this commit on this machine are skipped."""
    import multiprocessing as mp

    import orth_diag as od

    root = out / GUARD_DIR if part == "guard" else out
    folder_ = root / "dm"
    folder_.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder_, st)
    seeds = k21.parse_seeds(a.seeds)
    bs = harness._t3_head()
    t3d = bs.t3d
    spec, dec, rules = k18.barred_spec(Path(a.decision))
    if list(rules["protocol"]["seeds"]) != list(SEEDS):
        raise RuntimeError("the frozen script's seeds are not vec.adoption.SEEDS")
    paths = od.context_paths(Path(a.data), with_targets=True)
    n_par = a.workers * a.threads
    by_name = {v: k for k, v in k18.DM_TARGETS.items()}
    with bs.thread_limits(n_par), tempfile.TemporaryDirectory() as tmp:
        ctx, _ = t3d.build_context("gata4", paths, rules, n_jobs=n_par, threads=a.threads, out=Path(tmp))
        if sorted(T.name for T in ctx.targets) != sorted(k18.DM_TARGETS.values()):
            raise RuntimeError(f"targets {[T.name for T in ctx.targets]}, expected {sorted(k18.DM_TARGETS.values())}")
        k18.quiet_selfcheck(t3d, ctx)
        genes = list(ctx.genes)
        phi_b, frec = k18.barred_factor(Path(a.run_dir), dec, spec, genes)
        load_k27(genes, kept_sha)
        pre = load_pregate(out, _K27["record"])
        run_arms = GUARDS if part == "guard" else ARMS
        ref_C, ref_lab, ref_P, gap = od.reference_counts(ctx, paths["wt"], ctx.preset["label_mode"])
        wt_path = Path(paths["wt"])
        series = bs.wt_revert_series(wt_path, {n: Path(a.data) / "T2" / f"{n}.h5ad" for n in k7.FULL_PANEL
                                               if n != wt_path.stem})
        phi, per, group, rrec = bs.wt_revert_factors(wt_path, series, genes, ctx.wt.pool.C)
        refw = bs.Reference(preset=ctx.preset_name, rules=rules, genes=genes, V=ctx.V, pool=ctx.wt.pool,
                            coords=np.zeros((ctx.wt.pool.C.shape[0], 3)))
        files = [*(Path(a.data) / f for f in od.WT_FILES.values()), *(Path(a.data) / f for f in od.KO_FILES.values()),
                 *series.values(), Path(a.decision), Path(a.run_dir) / "factors" / f"{od.BARRED}.csv"]
        setup = {"world": "dm", "part": part, "barred_factor": frec,
                 "wt_revert": {k: v for k, v in rrec.items() if k != "series"},
                 "reference_counts_largest_gap": gap, "pregate": _pregate_record(pre),
                 "files_sha256": {str(f): t3d.sha256(f) for f in dict.fromkeys(files)}}
        sp = folder_ / "setup.json"
        if sp.exists():
            old = json.loads(sp.read_text())
            if {k: old.get(k) for k in setup} != json.loads(json.dumps(setup, default=float)):
                sys.exit(f"{sp} records another set-up (inputs, factor, tables or pre-gate): not run")
        else:
            k18.write_task_file(sp, setup | {"provenance": st})
        with k18_arms(run_arms):
            F = k18.fold_factors(bs, refw, per, phi, group, np.asarray(ctx.wt.X, dtype=np.float64).mean(0), phi_b,
                                 spec)
            if part != "guard":
                check_world(pre, "dm", F)
            wf = _world_figures(F, run_arms, SEEDS)
            k18._STATE.clear()
            k18._STATE.update(ctx=ctx, phi_b=phi_b, ref_counts=ref_C, ref_lab=ref_lab, ref_props=ref_P, F=F,
                              threads=a.threads)
            todo = [(ti, s) for ti, T in enumerate(ctx.targets) for s in seeds
                    if not (folder_ / f"{by_name[T.name]}_seed_{s}.json").exists()]

            def write(row: dict) -> None:
                k22._stamp_ok(st)
                path = folder_ / f"{by_name[row['target']]}_seed_{row['seed']}.json"
                k18.write_task_file(path, row | {"world": "dm", "part": part, "task": by_name[row["target"]],
                                                 "arms": list(run_arms), "pregate": _pregate_record(pre),
                                                 "provenance": st, "world_figures": wf, "task_source": SCRIPT})
                print(f"wrote {path}", flush=True)

            if a.workers > 1 and len(todo) > 1:
                if not bs.openmp_limited():
                    raise RuntimeError("the OpenMP runtime is not limited to one thread: forked workers would deadlock")
                with mp.get_context("fork").Pool(min(a.workers, len(todo))) as p:
                    for row in p.imap_unordered(k18._dm_task, todo):
                        write(row)
            else:
                for t in todo:
                    write(k18._dm_task(t))


def cmd_run(a) -> None:
    if a.world == "plain" and a.task not in TASKS:
        refuse("--world plain needs --task g1 or g2")
    out = Path(a.out)
    _, st = preconditions(out)
    kept_sha = cq.load_kept().sha256
    if a.part == "guard":
        try:
            guard_refs()
        except Refused as e:
            refuse(str(e))
    else:
        guard_passed(out, st)
    out.mkdir(parents=True, exist_ok=True)
    try:
        (run_plain if a.world == "plain" else run_dm)(a, st, out, a.part, kept_sha)
    except Exception:  # noqa: BLE001  (any crash: masked traceback, exit 1; seeds written so far stay)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(EXIT_CRASH)


# ---------------------------------------------------------------------------------------------------------------------
# The pre-gate command


def _fmt(v) -> str:
    return f"{v:.6g}" if isinstance(v, float) else str(v)


def _figs(r: dict) -> str:
    return ", ".join(f"{k} {_fmt(v)}" for k, v in r.items() if not isinstance(v, (dict, list)))


def pregate_lines(pre: dict) -> list[str]:
    c = pre["corr"]
    lines = [f"pre-gate P3 (counted; the tables, phi_M, `two` and the worlds' bases alone; no score): {pre['n_genes']} "
             f"genes; S* {pre['s_star']!r}; corr(ln C_Gata4, ln phi_M) {_fmt(c['ln_phi_m'])}, corr(ln C_Gata4, ln "
             f"two ** 2) {_fmt(c['ln_two_sq'])} (each over genes != 1 in either)",
             f"  x_e2: {_figs(pre['x_e2'])}"]
    for s, v in pre["arms"].items():
        lines.append(f"  {SIZE_SHOWN[s]}: e {v['e']!r}, a {v['a']!r}, k {v['k']!r}, den shift against phi_M "
                     f"{v['den_shift_vs_phi_m']:+.4%}, against x_e2 {v['den_shift_vs_x_e2']:+.4%}")
        lines.append(f"    scaled table: {_figs(v['table'])}")
        lines.append(f"    product: {_figs(v['product'])}")
        for w, cc in v["controls"].items():
            lines.append(f"    chip_perm ({w} basis): {_figs(cc)}")
        lines.append(f"    size twin: {_figs(v['twin'])}")
        b = pre["k24_bounds"][s]
        lines.append("    K24's bounds (recorded, not applied): " + ", ".join(
            f"{k} {_fmt(x['value'])} beyond {x['beyond']}" for k, x in b.items() if isinstance(x, dict)))
    lines.append("  placebos' e by pool index: " + "; ".join(
        f"{i}: S* {v['e']['s1']!r}, 4 S* {v['e']['s4']!r}" for i, v in pre["placebos"].items()))
    lines.append("  dropped: " + (pre["dropped"] or "no"))
    return lines


def cmd_pregate(a) -> None:
    """P3, COUNTED from its start: phi_M through the builder's loader (WT E8.75's gene order), the 11 tables and `two`
    through this module's loaders (every refusal before any figure), `two` checked against K23's committed pre-gate,
    each world's basis (K23's) and the den-weighted totals on WT E8.75's mean within-cell proportions; R/pregate.json.
    Exit 5 when either correlation drops the arm."""
    out = Path(a.out)
    rid, st = preconditions(out)
    path = out / "pregate.json"
    if path.exists():
        refuse(f"{path} exists")
    try:
        guard_refs()
        k23_pre, _ = k23_committed()
        p2, p2_sha = load_p2(out)
    except Refused as e:
        refuse(str(e))
    try:
        from vec.prep import load_stage

        bs = harness._t3_head()
        rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
        if rec.get("sha256") != bof.TABLE_SHA256:
            refuse("the barred factor is not transfer_mab21l2_all at its pinned sha256: refused")
        phi_m, genes = k18._BARRED["phi"], k18._BARRED["genes"]
        wt = load_stage(Path(a.data) / "T2" / "E8.75.h5ad", drop_duplicates=True)
        if list(wt.var_names) != list(genes):
            raise RuntimeError("WT E8.75's genes are not the barred factor's")
        C = np.asarray(bs.t3d.dense_counts(wt), dtype=np.float64)
        den = (C / np.maximum(C.sum(1, keepdims=True), 1)).mean(0)
        load_k27(genes, cq.load_kept().sha256)
        bases = k23.world_bases(Path(a.data), list(genes), bs)
        try:
            two_k23 = check_two(phi_m, _K27["two"], _K27["two_loaded"], bases, den, list(genes), k23_pre["pregate"])
        except Refused as e:
            refuse(str(e))
        pool_index = {n: r["pool_index"] for n, r in _K27["table_records"].items() if r["seed"] is not None}
        pre = pregate_k27(phi_m, _K27["two"], _K27["T"], _K27["e"], bases, den, list(genes), pool_index)
        res = {"pregate": pre, "arms": list(ARMS), "constants": constants_of(pre), "ext": _jsonable(_K27["record"]),
               "table_records": _K27["table_records"], "two_at_e2_as_k23": two_k23, "row": rid,
               "barred_factor": rec, "tries_counted_from_this_step": N_TRIES,
               "p2": {"sha256": p2_sha, "commit": (p2.get("provenance") or {}).get("commit")},
               "provenance": st}
        k18.write_task_file(path, res)
    except SystemExit:
        raise
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(EXIT_CRASH)
    print("\n".join(pregate_lines(pre)), flush=True)
    if dropped(pre):
        print("the arm and its controls are dropped at both sizes: nothing runs; the eight tries count", flush=True)
        sys.exit(EXIT_DROP)


def cmd_check(a) -> None:
    """The 11 tables against the listing, the pins, their manifests, the fetch manifest and the gene order; `two`
    against K22's pins; the exponents against S. Exit 0 or 2; prints one line."""
    out = Path(a.runs) if a.runs else None
    preconditions(out)
    try:
        with contextlib.redirect_stdout(sys.stderr):
            genes = cq.barred_genes(Path(a.data), Path(a.run_dir), Path(a.decision))
    except (Refused, cq.Refused) as e:
        refuse(str(e))
    except Exception as e:  # noqa: BLE001
        refuse(f"{type(e).__name__} in the loader")
    load_k27(genes)
    print(f"tables: {len(TABLE_NAMES)} checked against the listing, the pins, their manifests and the gene order",
          flush=True)


# ---------------------------------------------------------------------------------------------------------------------
# The guard part


def guard_files() -> list[str]:
    return ([f"{GUARD_DIR}/plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"{GUARD_DIR}/dm/{n}_seed_{s}.json" for n in TASKS for s in SEEDS])


def expected_files() -> list[str]:
    return ([f"plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"dm/{n}_seed_{s}.json" for n in TASKS for s in SEEDS])


def _rows_of(runs: Path, files: list[str], prefix: str, want_arms: tuple[str, ...], why: list[str], sha: dict,
             prov: dict, want_rec: dict | None = None) -> dict:
    """Rows per cell from task files (plain: rows[0]; dm: the row), each checked for its task, seed and arms, and
    (with want_rec) for the pre-gate's arms, sizes, constants and tables."""
    rows: dict = {(w, n): [] for w in WORLDS for n in TASKS}
    for f in files:
        r, sha[f] = k21._read_json(Path(runs) / f)
        prov[f] = r.get("provenance") or {}
        rel = f[len(prefix):]
        world, rest = ("plain", rel.split("/")[0][len("plain_"):]) if rel.startswith("plain_") else ("dm", rel[3:5])
        seed = int(rel.rsplit("_", 1)[1][:-len(".json")])
        if want_rec is not None:
            rec = r.get("pregate") or {}
            for k in ("arms", "sizes", "constants", "ext"):
                if not _same(rec.get(k), want_rec[k]):
                    why.append(f"{f}: its {k} differ from pregate.json's")
        if world == "plain":
            if r["task"]["key"] != TASKS[rest].key or [q["seed"] for q in r["rows"]] != [seed]:
                why.append(f"{f}: another task or seed")
                continue
            row = r["rows"][0]
            want = (harness.FLOOR, *want_arms)
        else:
            if r.get("task") != rest or r.get("seed") != seed:
                why.append(f"{f}: another task or seed")
                continue
            row = r
            want = want_arms
        miss = [x for x in want if x not in row.get("scores", {}) or x not in row.get("skills", {})
                or x not in row.get("t3", {})]
        if miss:
            why.append(f"{f}: missing arms {miss}")
            continue
        rows[(world, rest)].append(row)
    return {k: sorted(v, key=lambda q: q["seed"]) for k, v in rows.items()}


def _gated(rows: list[dict], a: str) -> int:
    return k21._gated([q["t3"][a]["severity"] for q in rows])


def guard_eval(rows: dict, refs: dict) -> dict:
    """GUARD 1: barred_form within K17's tolerance of K18's committed rows in all four cells. GUARD 2: x_e2 within the
    same tolerance of K23's per-seed x_e2 and of its committed ten-seed mean, its gated count equal to K23's, in all
    four cells."""
    g1, g2 = {}, {}
    for w in WORLDS:
        for n in TASKS:
            rr = rows.get((w, n)) or []
            seeds_ok = [q["seed"] for q in rr] == list(SEEDS)
            d1 = k18.distance(rr, BASE, refs["k18"].get((w, n)) or [], BASE) if seeds_ok else None
            g1[f"{w} {n}"] = {"distance": d1, "passes": bool(seeds_ok and k18.guard_passes(d1))}
            d2 = k18.distance(rr, REPRO, refs["k23_rows"].get((w, n)) or [], REPRO) if seeds_ok else None
            cell = refs["k23_cells"][(w, n)]
            m = float(np.mean([q["scores"][REPRO] for q in rr])) if rr else None
            gated = _gated(rr, REPRO) if rr else None
            ok = bool(seeds_ok and k18.guard_passes(d2) and m is not None
                      and abs(m - cell["mean"]) <= GUARD_MEAN and gated == cell["severity_gated"])
            g2[f"{w} {n}"] = {"distance": d2, "mean": m, "k23_mean": cell["mean"], "gated": gated,
                              "k23_gated": cell["severity_gated"], "passes": ok}
    failed = ([BASE] if not all(v["passes"] for v in g1.values()) else []) + \
        ([REPRO] if not all(v["passes"] for v in g2.values()) else [])
    return {"guard1": g1, "guard2": g2, "failed": failed, "passes": not failed}


def cmd_guard(a) -> None:
    """The guard check over R/guard/: R/guard.json; exit 0 when both guards hold, 6 when either fails (VOID), 3 when a
    file is missing."""
    runs = Path(a.runs)
    rid, st = preconditions(runs)
    p = runs / GUARD_FILE
    if p.exists():
        refuse(f"{p} exists: the guard check is made once")
    missing = [f for f in guard_files() if not (runs / f).exists()]
    if missing:
        print(f"INCOMPLETE: {len(missing)} guard file(s) missing, e.g. {missing[:4]}", flush=True)
        sys.exit(EXIT_INCOMPLETE)
    try:
        refs = guard_refs()
    except Refused as e:
        refuse(str(e))
    try:
        why: list[str] = []
        sha, prov = {}, {}
        rows = _rows_of(runs, guard_files(), f"{GUARD_DIR}/", GUARDS, why, sha, prov)
        bad_prov = [f for f, q in prov.items() if q.get("tree_clean") is not True or q.get("commit") != st["commit"]
                    or q.get("machine") != st["machine"]]
        if bad_prov or why:
            refuse(f"the guard files are not this commit's and machine's, or not well formed ({(bad_prov + why)[:4]})")
        res = guard_eval(rows, refs)
        body = {**res, "row": rid, "files_sha256": sha, "provenance": st}
        k18.write_task_file(p, body)
    except SystemExit:
        raise
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(EXIT_CRASH)
    print(f"guard part: {'passes' if res['passes'] else 'FAILS'} (failed guard arms: {res['failed'] or 'none'})",
          flush=True)
    if not res["passes"]:
        sys.exit(EXIT_GUARD)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _cells():
    return [(w, n) for w in WORLDS for n in TASKS]


def summarise(rows: list[dict], arms: tuple[str, ...]) -> dict:
    """Per arm: mean and sd of the score, mean skill per metric, severity states, gated seeds. Paired (k7._paired):
    every try arm against x_e2; chip_q of each size against its three controls; chip_q(4 S*) - chip_q(S*); every arm
    against barred_form (reported); barred_form against the floor (plain). Per-term skill differences of every try arm
    against x_e2 and against barred_form; the de_score plus de_direction skill change of every try arm against x_e2."""
    out: dict = {"arms": {}, "paired": {}, "per_term": {}, "per_term_vs_barred": {}, "de_change": {}}
    for x in arms:
        s = np.array([q["scores"][x] for q in rows], dtype=float)
        sev = [q["t3"][x]["severity"] for q in rows]
        out["arms"][x] = {
            "mean": float(s.mean()), "sd": float(s.std(ddof=1)),
            "skills": {m: float(np.mean([q["skills"][x][m] for q in rows])) for m in rows[0]["skills"][x]},
            "severity_states": {v: sum(z["state"] == v for z in sev) for v in sorted({z["state"] for z in sev})},
            "severity_gated": k21._gated(sev)}
    pairs = [(x, REPRO) for x in TRY_ARMS if x in arms]
    for s in SIZE_KEYS:
        pairs += [(arm("chip_q", s), arm(c, s)) for c in CONTROLS]
    pairs.append((arm("chip_q", "s4"), arm("chip_q", "s1")))
    pairs += [(x, BASE) for x in arms if x not in (BASE, harness.FLOOR)]
    if harness.FLOOR in arms:
        pairs.append((BASE, harness.FLOOR))
    out["paired"] = {f"{x} - {y}": k7._paired(rows, x, y) for x, y in pairs if x in arms and y in arms}
    for x in TRY_ARMS:
        if x not in arms:
            continue
        ms = rows[0]["skills"][REPRO]
        out["per_term"][x] = {m: float(np.mean([q["skills"][x][m] - q["skills"][REPRO][m] for q in rows])) for m in ms}
        out["per_term_vs_barred"][x] = {m: float(np.mean([q["skills"][x][m] - q["skills"][BASE][m] for q in rows]))
                                        for m in rows[0]["skills"][BASE]}
        out["de_change"][x] = float(np.mean([sum(q["skills"][x][m] - q["skills"][REPRO][m] for m in DE_TERMS)
                                             for q in rows]))
    return out


def _beats(v: dict | None) -> bool:
    return bool(v is not None and v["mean"] > MIN_GAIN and v["ahead"] >= MIN_AHEAD)


def size_status(per: dict, s: str) -> dict:
    """Clauses (a) to (g) of one size, against its own three controls (row K27, QUALIFIED)."""
    q = arm("chip_q", s)
    key = lambda c: f"{c[0]} {c[1]}"  # noqa: E731
    pair = lambda c, x, y: per[c[0]][c[1]]["paired"].get(f"{x} - {y}")  # noqa: E731
    gated_of = lambda c, x: per[c[0]][c[1]]["arms"][x]["severity_gated"]  # noqa: E731
    d = {c: pair(c, q, REPRO) for c in _cells()}
    ahead = {c: _beats(v) for c, v in d.items()}
    behind = {c: bool(v["mean"] < -MIN_GAIN) for c, v in d.items()}
    gated = {c: (gated_of(c, q), gated_of(c, REPRO)) for c in _cells()}
    up = [c for c in _cells() if ahead[c]]
    a = (any(ahead[("plain", n)] for n in TASKS) and any(ahead[("dm", n)] for n in TASKS)
         and any(ahead[(w, "g1")] for w in WORLDS) and any(ahead[(w, "g2")] for w in WORLDS)
         and any(ahead[("dm", n)] and gated[("dm", n)][0] == gated[("dm", n)][1] for n in TASKS))
    b = not any(behind.values())
    cc = all(g <= g0 + MAX_EXTRA_GATED for g, g0 in gated.values())
    dd = (all(_beats(pair(c, q, arm("chip_rand", s))) for c in up)
          and any(c[0] == "dm" and _beats(pair(c, q, arm("chip_rand", s))) for c in up))
    d2 = all(_beats(pair(c, q, arm("chip_perm", s))) for c in up)
    ee = all(_beats(pair(c, q, arm("chip_size", s))) for c in up)
    ff = all(per[c[0]][c[1]]["de_change"][q] > 0 for c in up)
    ctrl_ahead = {key(c): [x for x in CONTROLS if _beats(pair(c, arm(x, s), REPRO))] for c in up}
    gg = not any(ctrl_ahead.values())
    qualified = bool(a and b and cc and dd and d2 and ee and ff and gg)
    return {"size": SIZE_SHOWN[s], "qualified": qualified, "a": bool(a), "b": bool(b), "c": bool(cc), "d": bool(dd),
            "d2": bool(d2), "e": bool(ee), "f": bool(ff), "g": bool(gg),
            "ahead_cells": [key(c) for c in up], "behind_cells": [key(c) for c in _cells() if behind[c]],
            "means_vs_x_e2": {key(c): float(v["mean"]) for c, v in d.items()},
            "gated_arm_x_e2": {key(c): list(v) for c, v in gated.items()},
            "de_change": {key(c): float(per[c[0]][c[1]]["de_change"][q]) for c in _cells()},
            "controls_ahead_of_x_e2_in_ahead_cells": ctrl_ahead}


def pivot(per: dict) -> dict:
    """K23's fixed pivot for two sizes: chip_q(4 S*) - chip_q(S*), paired by seed, meets the AHEAD bar in at least one
    plain cell and at least one dm cell and has a mean below -0.3 in no cell."""
    k = f"{arm('chip_q', 's4')} - {arm('chip_q', 's1')}"
    d = {c: per[c[0]][c[1]]["paired"][k] for c in _cells()}
    ok = (any(_beats(d[("plain", n)]) for n in TASKS) and any(_beats(d[("dm", n)]) for n in TASKS)
          and not any(v["mean"] < -MIN_GAIN for v in d.values()))
    return {"holds": bool(ok), "means": {f"{w} {n}": float(v["mean"]) for (w, n), v in d.items()},
            "ahead": {f"{w} {n}": _beats(v) for (w, n), v in d.items()}}


def decide(per: dict) -> dict:
    """EITHER SIZE QUALIFYING: each size on its own clauses (a) to (g); PASS when at least one qualifies. WHICH SIZE A
    PASS NAMES: only 4 S*, 4 S*; only S*, S*; both, 4 S* when the pivot holds, else S*."""
    status = {s: size_status(per, s) for s in SIZE_KEYS}
    q = [s for s in SIZE_KEYS if status[s]["qualified"]]
    pv = pivot(per)
    if not q:
        named, rule, outcome = None, "nothing qualifies at either size", S_FAIL
    elif q == ["s4"]:
        named, rule = "s4", "only 4 S* qualifies"
    elif q == ["s1"]:
        named, rule = "s1", "only S* qualifies"
    else:
        named = "s4" if pv["holds"] else "s1"
        rule = "both qualify; " + ("the pivot holds: 4 S*" if pv["holds"] else "the pivot does not hold: S*")
    if named is not None:
        outcome = S_PASS.format(size=SIZE_SHOWN[named])
    return {"status": status, "qualified": [SIZE_SHOWN[s] for s in q], "pivot": pv, "rule": rule,
            "passes": named is not None, "named_size": SIZE_SHOWN[named] if named else None,
            "named_size_value": SIZES[named] if named else None, "outcome": outcome}


INCOMPLETE = ("Try K27: not a reading (INCOMPLETE): a file of the run is missing (a crash, a refusal such as a dirty "
              "tree, out of memory or a kill); nothing was read. The eight tries stand counted from P3 when P3 ran; a "
              "rerun past it happens only under rows K21 to K24's read-crash and code-fix rules and makes no new "
              "reading and no new count")


def reading_code() -> dict:
    st = k18.stamp()
    h = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()  # noqa: E731
    return {"script": SCRIPT, "script_sha256": h(__file__), "chipq_sha256": h(T3X / "chipq.py"),
            "chipq_fetch_sha256": h(T3X / "chipq_fetch.py"), "outside_forms_sha256": h(HERE / "outside_forms.py"),
            "k22_sha256": h(HERE / "k22_outside_forms.py"), "k23_sha256": h(HERE / "k23_outside_exponents.py"),
            "k21_sha256": h(HERE / "k21_factor_forms.py"), "k18_sha256": h(HERE / "k18_built_files.py"),
            "commit": st["commit"], "tree_clean": st["tree_clean"], "machine": st["machine"]}


def _p1_record(runs: Path) -> dict | None:
    """P1's two integers and its exit code, as the launcher saved them (R/p1.out, R/p1.exit); None when absent."""
    out, ex = Path(runs) / "p1.out", Path(runs) / "p1.exit"
    if not ex.is_file():
        return None
    rec: dict = {"exit": ex.read_text().strip()}
    m = re.fullmatch(r"(\d+) (\d+)\n?", out.read_text()) if out.is_file() else None
    rec |= {"kept": int(m.group(1)), "bound": int(m.group(2))} if m else {"kept": None, "bound": None}
    return rec


def _p2_record(runs: Path) -> tuple[dict | None, str | None]:
    p = Path(runs) / "p2.json"
    if not p.is_file():
        return None, None
    raw = p.read_bytes()
    return json.loads(raw), sha256_bytes(raw)


def _tables_ok(erec: dict, k23_read: dict | None) -> list[str]:
    why = []
    if {"listing_sha256": erec.get("listing_sha256"), "tables": erec.get("tables")} != pinned_ext():
        why.append("pregate.json: its tables are not the pinned ones (listing or table sha256)")
    two = (erec.get("two") or {})
    if two.get("listing_sha256") != ext.K22_LISTING_SHA256 or (two.get("tables") or {}).get("all") != \
            ext.K22_TABLE_SHA256.get("all"):
        why.append("pregate.json: `two` is not at K22's pins")
    if k23_read is not None:
        ke = k23_read.get("ext") or {}
        if (ke.get("tables") or {}).get("all") != (two.get("tables") or {}).get("all") or \
                ke.get("listing_sha256") != two.get("listing_sha256") or \
                ke.get("combine_sha256") != two.get("combine_sha256"):
            why.append("pregate.json: `two` is not the table K23's committed reading.json records")
    else:
        why.append("K23's committed reading.json is missing")
    return why


def _void(common: dict, clauses: list[str], tries: int) -> dict:
    clause = clauses[0].replace(";", ",") if clauses else "unrecorded"
    return {"incomplete": False, **common, "void": True, "why": clauses, "tries_counted": tries,
            "outcome": S_VOID.format(clause=clause), "passes": False}


def reading(runs: Path, refs: Path = K18_TASKS, u3d: Path = k18.U3D_RECORDS, k23_dir: Path = K23_RESULTS) -> dict:
    """Apply the reading of row K27 to a run folder: a STOP at P1 or P2 (nothing counted), a drop at P3 (the eight
    tries), a VOID (the eight tries plus each failed guard arm), INCOMPLETE (not a reading), or the clauses read."""
    runs = Path(runs)
    p1 = _p1_record(runs)
    p2, p2_sha = _p2_record(runs)
    base_rec = {"p1": p1, "p2": None if p2 is None else {k: p2.get(k) for k in (
        "query_kept", "query_bound", "query_K", "query_genes_changed", "pool_size", "placebos",
        "barred_gene_antigen", "passes", "why", "exponents", "kept_sha256", "listing_sha256", "tables")},
                "p2_sha256": p2_sha}
    if not (runs / "pregate.json").exists():
        if p1 is not None and p1["exit"] == "3":
            return {"incomplete": False, "void": False, "stop": "P1", **base_rec, "tries_counted": 0,
                    "outcome": S_P1, "passes": False}
        if p2 is not None and p2.get("passes") is False:
            return {"incomplete": False, "void": False, "stop": "P2", **base_rec, "tries_counted": 0,
                    "outcome": S_P2, "passes": False}
        return {"incomplete": True, "void": False, "missing": ["pregate.json"], "outcome": INCOMPLETE, **base_rec}
    if pins_state() != "filled":
        raise Refused("a pre-gate is in the run folder but the pins are not filled: refused")
    why: list[str] = []
    sha, prov = {}, {}
    pre_file, sha["pregate.json"] = k21._read_json(runs / "pregate.json")
    prov["pregate.json"] = pre_file.get("provenance") or {}
    pre = pre_file["pregate"]
    erec = pre_file.get("ext") or {}
    try:
        _, k23_read = k23_committed(k23_dir)
    except Refused:
        k23_read = None
    why += _tables_ok(erec, k23_read)
    if pre_file.get("constants") != constants_of(pre):
        why.append("pregate.json: its constants are not its pre-gate's")
    if pre_file.get("arms") != list(ARMS):
        why.append("pregate.json: its arms are not this row's")
    if (pre_file.get("p2") or {}).get("sha256") != p2_sha:
        why.append("p2.json is not the one the pre-gate read")
    if p2 is not None and (p2.get("provenance") or {}).get("commit") == prov["pregate.json"].get("commit"):
        why.append("the pre-gate ran at the registration commit, not the row's second commit")
    want_rec = {"arms": list(ARMS), "ext": erec, "constants": constants_of(pre), "sizes": dict(SIZES)}
    common = {"files_sha256": sha, "provenance": prov, "arms": list(ARMS), "pregate": pre, "ext": erec, **base_rec}
    if dropped(pre):
        if prov["pregate.json"].get("tree_clean") is not True:
            why.append("pregate.json: a dirty tree")
        if why:
            return _void(common, why, N_TRIES)
        return {"incomplete": False, **common, "void": False, "stop": "P3", "tries_counted": N_TRIES,
                "outcome": S_DROP, "passes": False}
    if not (runs / GUARD_FILE).exists():
        return {"incomplete": True, "void": False, "missing": [GUARD_FILE], "outcome": INCOMPLETE, **base_rec}
    g, sha[GUARD_FILE] = k21._read_json(runs / GUARD_FILE)
    prov[GUARD_FILE] = g.get("provenance") or {}
    gmiss = [f for f in guard_files() if not (runs / f).exists()]
    if gmiss:
        return {"incomplete": True, "void": False, "missing": gmiss, "outcome": INCOMPLETE, **base_rec}
    gwhy: list[str] = []
    gsha, gprov = {}, {}
    grows = _rows_of(runs, guard_files(), f"{GUARD_DIR}/", GUARDS, gwhy, gsha, gprov, want_rec)
    if gsha != g.get("files_sha256"):
        why.append("guard.json: its files' sha256 are not the guard files'")
    why += gwhy
    sha |= gsha
    prov |= gprov
    common |= {"guard_part": {k: g.get(k) for k in ("guard1", "guard2", "failed", "passes")}}
    if g.get("passes") is not True:
        failed = [x for x in (g.get("failed") or []) if x in GUARDS]
        return _void(common, [f"the guard part failed ({', '.join(failed) or 'unrecorded'})", *why],
                     try_count(failed))
    missing = [f for f in expected_files() if not (runs / f).exists()]
    if missing:
        return {"incomplete": True, "void": False, "missing": missing, "outcome": INCOMPLETE, **base_rec}
    rows = _rows_of(runs, expected_files(), "", ARMS, why, sha, prov, want_rec)
    for f, p in prov.items():
        if p.get("tree_clean") is not True:
            why.append(f"{f}: a dirty tree")
    commits = {p.get("commit") for p in prov.values()}
    machines = {p.get("machine") for p in prov.values()}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(map(str, commits))}")
    if len(machines) > 1:
        why.append(f"results from different machines {sorted(map(str, machines))}")
    saved = dict(k21.K18_SHA256)
    k21.K18_SHA256.clear()
    k21.K18_SHA256.update(K18_SHA256)
    try:
        _, kwhy = k21.k18_rows(refs)
    finally:
        k21.K18_SHA256.clear()
        k21.K18_SHA256.update(saved)
    why += kwhy
    per: dict = {"plain": {}, "dm": {}}
    match: dict = {}
    for (w, n), rr in rows.items():
        if [q["seed"] for q in rr] != list(SEEDS):
            why.append(f"{w} {n}: not the ten seeds")
            continue
        per[w][n] = summarise(rr, (harness.FLOOR, *ARMS) if w == "plain" else ARMS)
        for x in GUARDS:
            d = k18.distance(rr, x, grows.get((w, n)) or [], x)
            match[f"{w} {n} {x}"] = d
            if not k18.guard_passes(d):
                why.append(f"guard match: {w} {n}: {x} in the scoring part is off the guard part's (limits "
                           f"{GUARD_MEAN} and {GUARD_SEED})")
    check = k18.dm_check({n: rows[("dm", n)] for n in TASKS}, u3d)
    common |= {"guard_match": match, "dm_check_u3d_report_only": check}
    if why:
        return _void(common, why, N_TRIES) | {"worlds": per}
    return ({"incomplete": False, **common, "void": False, "why": [], "tries_counted": N_TRIES,
             "commit": next(iter(commits)), "machine": next(iter(machines)), "worlds": per} | decide(per))


def text(res: dict) -> str:
    """reading.txt: the outcome sentence on a line of its own beginning OUTCOME, then the reported figures."""
    lines = [f"OUTCOME {res['outcome']}"]
    if res.get("incomplete"):
        miss = res.get("missing") or []
        return "\n".join([*lines, f"missing: {len(miss)} file(s), e.g. " + ", ".join(miss[:6])]) + "\n"
    lines.append(f"tries counted: {res['tries_counted']}")
    p1 = res.get("p1")
    lines.append("P1 (Ctnnb1, database only): " + ("not recorded" if p1 is None else
                                                    f"kept {p1['kept']}, bound set {p1['bound']}, exit {p1['exit']}"))
    p2 = res.get("p2")
    if p2 is not None:
        lines.append(f"P2: practice query kept {p2['query_kept']}, bound set {p2['query_bound']}, K {p2['query_K']}, "
                     f"genes changed {p2['query_genes_changed']}; pool {p2['pool_size']}; the analysis list holds an "
                     f"antigen for the barred factor's gene: {p2['barred_gene_antigen']}; passes {p2['passes']}"
                     + (f" ({'; '.join(p2['why'])})" if p2.get("why") else ""))
        for s, v in sorted((p2.get("placebos") or {}).items(), key=lambda kv: int(kv[0])):
            lines.append(f"  P2 seed {s}: pool index {v['pool_index']}, kept {v['kept']}, bound set {v['bound']}, K "
                         f"{v['K']}, genes changed {v['genes_changed']}")
    if res.get("stop") in ("P1", "P2"):
        return "\n".join(lines) + "\n"
    lines.append(f"arms: {', '.join(res['arms'])}")
    e = res.get("ext") or {}
    lines.append(f"tables: listing {e.get('listing_sha256')}; fetch manifest {e.get('fetch_manifest_sha256')}; kept "
                 f"list {e.get('kept_sha256')}; " + ", ".join(f"{n} {s}" for n, s in (e.get("tables") or {}).items()))
    lines.append(f"`two`: {(e.get('two') or {})}")
    lines += pregate_lines(res["pregate"])
    gp = res.get("guard_part") or {}
    if gp:
        lines.append(f"guard part: passes {gp.get('passes')}, failed {gp.get('failed') or 'none'}")
        for name in ("guard1", "guard2"):
            for k, v in (gp.get(name) or {}).items():
                extra = (f", mean {v['mean']:.2f} against {v['k23_mean']:.2f}, gated {v['gated']} against "
                         f"{v['k23_gated']}" if name == "guard2" and v.get("mean") is not None else "")
                lines.append(f"  {name} {k}: passes {v['passes']}; {k21._dist(v.get('distance'))}{extra}")
    if res.get("void"):
        lines.append("VOID clauses: " + "; ".join(res["why"]))
    elif res.get("status"):
        lines.append(f"commit {res['commit']}, machine {res['machine']}.")
        for st in res["status"].values():
            lines.append(f"  {st['size']}: qualified {'yes' if st['qualified'] else 'no'} ((a) {st['a']}, (b) "
                         f"{st['b']}, (c) {st['c']}, (d) {st['d']}, (d2) {st['d2']}, (e) {st['e']}, (f) {st['f']}, "
                         f"(g) {st['g']}); ahead in {st['ahead_cells'] or 'none'}; behind in "
                         f"{st['behind_cells'] or 'none'}; controls ahead of x_e2 in its ahead cells "
                         f"{st['controls_ahead_of_x_e2_in_ahead_cells'] or 'none'}; gated (arm, x_e2) "
                         f"{st['gated_arm_x_e2']}")
        pv = res["pivot"]
        lines.append(f"  pivot chip_q(4 S*) - chip_q(S*): holds {pv['holds']}; ahead {pv['ahead']}")
        lines.append(f"  qualified: {res['qualified'] or 'none'}; rule: {res['rule']}; named size: "
                     f"{res['named_size'] or 'none'}")
    for k, d in (res.get("guard_match") or {}).items():
        lines.append(f"  scoring part against the guard part, {k}: {k21._dist(d)}")
    if res.get("dm_check_u3d_report_only"):
        lines.append("U3d's check arm in the dm world, reported only: " + "; ".join(
            f"{n} {k21._dist(c['distance'])}" for n, c in res["dm_check_u3d_report_only"].items()))
    lines.append("task files: " + ", ".join(f"{f} {h[:12]}" for f, h in (res.get("files_sha256") or {}).items()))
    lines.append("g1 and g2 are two replicates of one knockout run; the controls, x_e2 and barred_form are never "
                 "named; K24's bounds and every reported figure never change the named size, the pass or the file.")
    for w in WORLDS:
        for n, t in TASKS.items():
            v = (res.get("worlds") or {}).get(w, {}).get(n)
            if v is None:
                continue
            lines.append(f"\n{w} {n} {t.key}:")
            for a, s in v["arms"].items():
                lines.append(f"    {a:<14} {s['mean']:.2f} (sd {s['sd']:.2f}); gated {s['severity_gated']} "
                             f"{s['severity_states']}; skills " + ", ".join(f"{m} {x:.3f}" for m, x in
                                                                          s["skills"].items()))
            for k, d in v["paired"].items():
                lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})")
            for a, terms in v["per_term"].items():
                lines.append(f"    per-term skill, {a} - {REPRO}: " + ", ".join(f"{m} {x:+.3f}" for m, x in
                                                                            terms.items())
                             + f"; de_score + de_direction {v['de_change'][a]:+.3f}")
    if res.get("reading_code"):
        rc = res["reading_code"]
        lines.append(f"\nread by {rc['script']} {rc['script_sha256'][:12]} (chipq {rc['chipq_sha256'][:12]}) at commit "
                     f"{rc['commit']}, clean {rc['tree_clean']}, machine {rc['machine']}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    out = Path(a.out)
    targets = [out / "reading.json", out / "reading.txt"]
    if any(p.exists() for p in targets):
        print(f"refused: {out} holds a reading already; it is never overwritten", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    try:
        res = reading(Path(a.runs))
    except Refused as e:
        refuse(str(e))
    try:
        if not res.get("incomplete"):
            res["reading_code"] = reading_code()
        t = text(res)
        js = json.dumps(res, indent=1, default=float)
    except Exception:  # noqa: BLE001  (masked traceback, nothing written)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(EXIT_CRASH)
    if res.get("incomplete"):
        print(t, end="")
        sys.exit(EXIT_INCOMPLETE)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(js)
    (out / "reading.txt").write_text(t)
    print(t, end="")


# ---------------------------------------------------------------------------------------------------------------------


def main(argv=None) -> None:
    argv = sys.argv[1:] if argv is None else list(argv)
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p2 = sub.add_parser("p2")
    p2.add_argument("--out", required=True)
    pn = sub.add_parser("pins")
    pn.add_argument("--runs", required=True)
    ck = sub.add_parser("check")
    ck.add_argument("--runs", default=None)
    p = sub.add_parser("pregate")
    p.add_argument("--out", required=True)
    r = sub.add_parser("run")
    r.add_argument("--part", choices=("guard", "score"), default="score")
    r.add_argument("--world", choices=WORLDS, required=True)
    r.add_argument("--task", choices=sorted(TASKS))
    r.add_argument("--seeds", default=None)
    r.add_argument("--out", required=True)
    r.add_argument("--workers", type=int, default=1)
    r.add_argument("--threads", type=int, default=8)
    for q in (p2, ck, p, r):
        q.add_argument("--data", default=str(harness.DATA_DIR))
        q.add_argument("--run-dir", default=str(k18.RUN_DIR))
        q.add_argument("--decision", default=str(k18.DECISION))
    gd = sub.add_parser("guard")
    gd.add_argument("--runs", required=True)
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"p2": cmd_p2, "pins": cmd_pins, "check": cmd_check, "pregate": cmd_pregate, "run": cmd_run, "guard": cmd_guard,
     "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
