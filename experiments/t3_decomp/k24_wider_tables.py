#!/usr/bin/env python
"""Try K24 (experiments/harness/TRIES.md; the row is placed before any run, and `pregate`, `run` and `guard` refuse
without it): the barred factor times each of two wider outside-knockout tables and times their equal-weight joint, each
brought to the table size of K23's named form by a target-free exponent, read against the barred s1 on the two Gata4
replicates in K18's two practice worlds. A labelled reading that may name at most one T3 candidate file, and only one
that beats K23's named form head to head, by the rule fixed in the row; not the adoption rule.

Built on outside_forms.py, k22_outside_forms.py, k23_outside_exponents.py, k18_built_files.py and k21_factor_forms.py,
all imported unchanged. The tables come only by KEY through the committed listing experiments/t3_ext/TABLES_K24.sha256
(keys bulk, sc3, two) at the pins below, which the main session fills in the row's commit (None refuses every run):

  two   K22's and K23's two-study single-cell table `all`, from the folder VEC_T3_EXT_TABLES_TWO names, against its
        COMBINE.json; used for the guard arm x_e2, for S* and for the corr-with-two drop only
  sc3   the three-study single-cell table `all`, from the folder VEC_T3_EXT_TABLES names, against COMBINE.json
  bulk  the bulk whole-embryo table `bulk_all` (study MGP-DMDD), same folder, against COMBINE_bulk.json

Each folder goes through outside_forms.table_dir's checks unchanged; the bytes are read once and must match the listing,
the pins, the COMBINE record, factors/factors.json and the manifest; the genes are the barred factor's (order included);
every value finite and positive and used as loaded.

The size rule: S* = 0.04949423311987974, K23's committed SD(ln two ** 2); for T in {sc3, bulk} e_T = S* / SD(ln T)
(population SD over all genes in file order, genes exactly 1 included); the joint J = sc3 ** (c e_sc3) * bulk ** (c
e_bulk), c = S* / SD(e_sc3 ln sc3 + e_bulk ln bulk). Nothing is rescaled or clipped after any product.

Arms (phi_M = transfer_mab21l2_all as loaded; perm = outside_forms.permuted within deciles of each world's basis, seed
20261113, 10 bins, an index permutation):
  barred_form    phi_M                       base and guard 1 (K18's arm value for value)
  x_e2           phi_M * two ** 2.0          guard 2, K23's named form through K23's code path (never named here)
  x_sc3          phi_M * sc3 ** e_sc3        try           x_sc3_perm   phi_M * perm(sc3 ** e_sc3)   control
  x_bulk         phi_M * bulk ** e_bulk      try           x_bulk_perm  phi_M * perm(bulk ** e_bulk) control
  x_joint        phi_M * J                   try           x_joint_perm phi_M * perm(J)              control
  x_<X>_size     k * phi_M ** a              size twins (labelled controls, counted, never named), a = SD(ln product)
                                             / SD(ln phi_M), k = den_total(product) / den_total(phi_M ** a), from the
                                             pre-gate (K23's formula)

  pregate --out DIR                                     the pre-gate: DIR/pregate.json and its lines (no score)
  run --part guard --world plain --task g1|g2 [--seeds] --out DIR   the guard part (barred_form and x_e2 only) into
                                                        DIR/guard/ (resumable)
  run --part guard --world dm [--workers 8 --threads 2] --out DIR
  guard --runs DIR                                      the guard check: DIR/guard.json (exit 6 when a guard fails)
  run [--part score] --world plain --task g1|g2 [--seeds] --out DIR   a scoring part (every arm that passed the
                                                        pre-gate); refused without a passing guard.json of this
                                                        commit and machine
  run [--part score] --world dm [--workers 8 --threads 2] --out DIR
  read --runs DIR --out DIR                             the fixed reading (never over an existing one: exit 4)
  listing | pins                                        TABLES_K24.sha256's text, and the two pin lines (main session)

Exit codes: 0 done; 1 a crash (masked traceback; seeds written so far stay) or a world whose controls are not the
pre-gate's (INCOMPLETE); 2 refused (no row, a dirty tree, pins not filled or `two`'s not K22's, the listing untracked
or off its pin, K23's committed files untracked, a folder variable unset or refused, bytes or records that disagree, a
missing or unpinned table, the frozen T3 files off their pins, `two`'s figures off K23's, a scoring part without a
passing guard.json); 3 INCOMPLETE (a file missing; nothing written); 4 read: the reading files exist already; 5
pregate: no try left (both tables dropped, or bulk dropped and sc3 a repeat of the named form), nothing runs; 6
guard: the guard part failed (VOID; the failed guard arms counted).

The try count (row K24 as restated on 2026-11-20): both tables pass 9; sc3 a repeat of the named form only and the
joint passes its own drops 6 (bulk 3, joint 3), the joint failing them 3; sc3 dropped for another reason 3; bulk
dropped and sc3 passing 3; bulk dropped and sc3 a repeat 0 (exit 5); both dropped 0 (exit 5); the joint dropped on its
own with both tables valid 6. Never more than nine.

K23's run folder is runs/k23_<first 7 of K23's committed reading's commit> in this checkout; no option overrides it.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import re
import sys
import tempfile
import traceback
from contextlib import contextmanager
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k7_wt_revert as k7  # noqa: E402  (the paired summary)
import k8_s2_form as k8  # noqa: E402  (the upload's central composition, for the record)
import k18_built_files as k18  # noqa: E402  (the worlds, tasks, stream, seeds, build path and scoring; unchanged)
import k21_factor_forms as k21  # noqa: E402  (the seed parser, the summary and K18's committed rows; unchanged)
import k22_outside_forms as k22  # noqa: E402  (frozen-file check, stamp check, constants; unchanged)
import k23_outside_exponents as k23  # noqa: E402  (x_e2's code path, its pre-gate and constants; unchanged)
import outside_forms as ext  # noqa: E402  (the table checks, the product and the permutation; unchanged)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

SCRIPT = "experiments/t3_decomp/k24_wider_tables.py"
ROW_SCRIPT = f"`{SCRIPT}`"  # the ledger row names this script
TRIES_MD = REPO / "experiments" / "harness" / "TRIES.md"  # looked up at call time (row_present)
LISTING_FILE = REPO / "experiments" / "t3_ext" / "TABLES_K24.sha256"
ENV_TABLES = ext.ENV  # "VEC_T3_EXT_TABLES"
ENV_TWO = "VEC_T3_EXT_TABLES_TWO"
KEYS = ("bulk", "sc3", "two")  # the listing's keys, sorted
SOURCES = {"two": (ENV_TWO, "all", "COMBINE.json"), "sc3": (ENV_TABLES, "all", "COMBINE.json"),
           "bulk": (ENV_TABLES, "bulk_all", "COMBINE_bulk.json")}
# Filled in the row's commit from the committed listing (`k24_wider_tables.py pins` prints both lines). None refuses
# every run, pre-gate and build.
K24_LISTING_SHA256: str | None = None
K24_TABLE_SHA256: dict[str, str | None] = {
    "two": "9490f524e40e146480c775f85ae8829ec9c005607104ad8ace129c5c46ed6228", "sc3": None, "bulk": None}

STREAM = k18.STREAM  # "t3_s2_form"
BASE = k18.BARRED  # "barred_form"
REPRO = "x_e2"  # K23's named form, guard 2
ARMS = (BASE, REPRO, "x_sc3", "x_sc3_perm", "x_sc3_size", "x_bulk", "x_bulk_perm", "x_bulk_size", "x_joint",
        "x_joint_perm", "x_joint_size")
GUARDS = (BASE, REPRO)
NAMEABLE = ("x_sc3", "x_bulk", "x_joint")  # also the tie order
KIND_OF = {"x_sc3": "sc3", "x_sc3_perm": "sc3", "x_sc3_size": "sc3", "x_bulk": "bulk", "x_bulk_perm": "bulk",
           "x_bulk_size": "bulk", "x_joint": "joint", "x_joint_perm": "joint", "x_joint_size": "joint"}
ARMS_OF_KIND = {"sc3": ("x_sc3", "x_sc3_perm", "x_sc3_size"), "bulk": ("x_bulk", "x_bulk_perm", "x_bulk_size"),
                "joint": ("x_joint", "x_joint_perm", "x_joint_size")}
CONTROL_OF = {"x_sc3": "x_sc3_perm", "x_bulk": "x_bulk_perm", "x_joint": "x_joint_perm"}
SIZE_OF = {"x_sc3": "x_sc3_size", "x_bulk": "x_bulk_size", "x_joint": "x_joint_size"}
MAX_TRIES = 9
TASKS = k18.TASKS
WORLDS = k18.WORLDS

# The size rule and the drops (row K24; fixed before any figure of sc3 at another exponent or of bulk existed)
S_STAR = 0.04949423311987974  # K23's committed pregate.json exponents['2']['table']['sd_ln'] (a test pins it)
E_MAX, E_MIN = 8.0, 0.125  # s1
MIN_GENES = 100  # s2: genes != 1
MAX_ABS_LN_SCALED = 0.5  # s3: max |ln T ** e_T|
JOINT_CORR = (0.0, 0.9)  # the joint runs when 0 <= corr(ln sc3, ln bulk) < 0.9
SC3_TWO_MAX_CORR = ext.FILE2_MAX_CORR  # 0.9: sc3 is dropped as a repeat of the named form at or above it
BULK_MIN_LINES = 12  # (b2)
BULK_MIN_LIBS_UNIT, BULK_MIN_LIBS_TOTAL = 2, 36  # (b3)
# (b4) as restated by the cold challenger on 2026-11-20, before any table value, agreement count or score was read:
# at least max(AGREE_FLOOR, E0 + AGREE_EXCESS G) of the G genes != 1 with an agreement value and held in at least
# AGREE_MIN_LINES lines agree at AGREE or more; E0 the exact sum of chance_agreement(n_g) over those genes
AGREE, AGREE_FLOOR, AGREE_EXCESS = 2.0 / 3.0, 50, 0.11
AGREE_MIN_LINES = BULK_MIN_LINES  # 12, the (b2) floor
AGREE_TOL = 1e-12  # a share of exactly 2/3 counts
BULK_SIGN_RULE = "single study"  # combine.py's sign rule of a one-study bulk table
BULK_STUDY = "MGP-DMDD"
BULK_MODALITY, BULK_TIER = "bulk", "unrelated"
REL_TOL = 1e-12
# K22's and K23's constants, from their modules
MIN_CHANGED, MAX_ABS_CORR, MIN_STUDIES = ext.MIN_CHANGED, ext.MAX_ABS_CORR, ext.MIN_STUDIES  # 10, 0.5, 2
MAX_AT_CLIP, MAX_DEN_SHIFT, MAX_SD_RATIO = k23.MAX_AT_CLIP, k23.MAX_DEN_SHIFT, k23.MAX_SD_RATIO  # 45, 0.05, 1.05
MIN_GAIN, MIN_AHEAD, MAX_EXTRA_GATED = k23.MIN_GAIN, k23.MIN_AHEAD, k23.MAX_EXTRA_GATED  # 0.3, 8 of 10, +1
GUARD_MEAN, GUARD_SEED = k23.GUARD_MEAN, k23.GUARD_SEED  # 0.1 and 0.5 (K17's tolerance)
K18_TASKS = k23.K18_TASKS
K18_SHA256 = dict(k23.K18_SHA256)
K23_RESULTS = k23.RESULTS
RESULTS = HERE / "results" / "k24_wider_tables"
EXIT_REFUSED, EXIT_INCOMPLETE, EXIT_EXISTS, EXIT_NOTHING_RUNS, EXIT_GUARD = 2, 3, 4, 5, 6
_K24: dict = {}  # the loaded tables (file order) and the pre-gate's constants, set before any arm builds
NO_FILE = "no file: no wider table beats the named form at a like size; the K23 file stands alone"
NOT_AHEAD = "qualified against the barred form, not ahead of the named form"


class Refused(Exception):
    """A refusal of the whole run (exit status 2)."""


def refuse(msg: str) -> None:
    print(f"refused: {msg}", file=sys.stderr, flush=True)
    sys.exit(EXIT_REFUSED)


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
# The listing, the folders and the tables


def check_pins(pin: str | None, pins: dict) -> None:
    if pin is None or set(pins) != set(KEYS) or any(pins.get(k) is None for k in KEYS):
        raise Refused("the pins K24_LISTING_SHA256 and K24_TABLE_SHA256 are not filled (k24_wider_tables.py; they go "
                      "in with row K24's commit, from the committed TABLES_K24.sha256): refused (a missing or "
                      "unpinned table stops the row)")
    if pins["two"] != ext.K22_TABLE_SHA256.get("all") or pins["two"] != k23.K23_TABLE_SHA256:
        raise Refused("the pin of `two` is not the `all` sha256 K22 and K23 pinned and read: refused")


def read_listing(path: Path | None = None, pin: str | None = ..., pins: dict | None = None,
                 check_tracked: bool = True) -> tuple[dict[str, str], str]:
    """(listing, its sha256), checked: pins filled, the file no symlink, tracked, at its pin, of the fixed form
    (outside_forms.parse_listing), keys exactly bulk, sc3 and two, each at its pin."""
    path = LISTING_FILE if path is None else Path(path)
    pin = K24_LISTING_SHA256 if pin is ... else pin
    pins = K24_TABLE_SHA256 if pins is None else pins
    check_pins(pin, pins)
    if path.is_symlink() or not path.is_file():
        raise Refused(f"{path.name} is missing or a symlink: refused")
    if check_tracked and not ext.tracked(path):
        raise Refused(f"{path.name} is not tracked by git (it is committed with row K24): refused")
    raw = path.read_bytes()
    sha = sha256_bytes(raw)
    if sha != pin:
        raise Refused(f"{path.name} hashes to {sha[:12]}, not the pinned K24_LISTING_SHA256 {pin[:12]}: refused")
    try:
        listing = ext.parse_listing(raw)
    except ext.TableRefused as e:
        raise Refused(f"{path.name}: {e}") from None
    if sorted(listing) != list(KEYS):
        raise Refused(f"{path.name} lists {sorted(listing)}, not exactly {list(KEYS)}: refused")
    for k in KEYS:
        if listing[k] != pins[k]:
            raise Refused(f"{k}: the listing says {listing[k][:12]}, the pin {str(pins[k])[:12]}: refused")
    return listing, sha


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


def combine_record(d: Path, fname: str) -> tuple[dict, str]:
    """A COMBINE file's pooled sha256 per table name and the file's own sha256 (outside_forms.combine_record's reading,
    for COMBINE.json or COMBINE_bulk.json: {'tables': {name: {name: {'sha256': ...}}}, 'record': ...})."""
    p = Path(d) / fname
    if p.is_symlink() or not p.is_file():
        raise Refused(f"{fname} is missing or a symlink: refused")
    raw = p.read_bytes()
    try:
        c = json.loads(raw)
    except ValueError:
        raise Refused(f"{fname} is not JSON: refused") from None
    if not isinstance(c, dict) or not isinstance(c.get("tables"), dict) or "record" not in c:
        raise Refused(f"{fname} is not combine.py's form {{'tables', 'record'}}: refused")
    out = {}
    for name, files in c["tables"].items():
        rec = files.get(name) if isinstance(files, dict) else None
        if not isinstance(rec, dict) or not isinstance(rec.get("sha256"), str):
            raise Refused(f"{fname} holds no pooled record for {name}: refused")
        out[name] = rec["sha256"]
    return out, sha256_bytes(raw)


RECORD_FIELDS = ("modality", "studies", "tiers", "sign_rule", "genes_by_study", "units_by_study", "libraries_by_unit",
                 "k_by_unit", "guards_by_unit", "line_agreement_by_study", "refstats_sha256", "checked")
OPTIONAL_FIELDS = ("lines_held_by_study",)  # combine.py writes it since the restated (b4); older manifests lack it


def _plain_file(p: Path, what: str) -> Path:
    if p.is_symlink() or not p.is_file():
        raise Refused(f"{what} ({p.name}) is missing or a symlink: refused")
    return p


def _json(p: Path, what: str) -> tuple[dict, str]:
    raw = _plain_file(p, what).read_bytes()
    try:
        return json.loads(raw), sha256_bytes(raw)
    except ValueError:
        raise Refused(f"{what} is not JSON: refused") from None


def load_k24_table(key: str, genes: list[str], listing: dict[str, str], listing_sha256: str,
                   pins: dict | None = None, env: dict | None = None) -> tuple[np.ndarray, dict, dict]:
    """The keyed table in `genes`' order (which must be the file's), its record and its COMBINE file's pooled sha256
    per name; mirrors outside_forms.load_table (see the module)."""
    import pandas as pd

    if key not in SOURCES:
        raise Refused(f"{key!r} is not one of K24's tables {list(KEYS)}: refused")
    pins = K24_TABLE_SHA256 if pins is None else pins
    var, name, cfile = SOURCES[key]
    want = listing.get(key)
    if want is None or pins.get(key) != want:
        raise Refused(f"{key}: not listed, or the listing's sha256 is not the pinned one: refused")
    d = folder(var, env)
    td, fd = d / name, d / name / "factors"
    for p, what in ((td, "the table folder"), (fd, "the factors folder")):
        if p.is_symlink() or not p.is_dir():
            raise Refused(f"{key}: {what} is missing or a symlink: refused")
    raw = _plain_file(fd / f"{name}.csv", f"{key}'s pooled table").read_bytes()  # read once: hashed bytes are parsed
    sha = sha256_bytes(raw)
    if sha != want:
        raise Refused(f"{key}: the file hashes to {sha[:12]}, not the listed {want[:12]}: refused")
    comb, comb_sha = combine_record(d, cfile)
    if comb.get(name) != sha:
        raise Refused(f"{key}: {cfile} records another sha256 for {name}: refused")
    fj, _ = _json(fd / "factors.json", f"{key}'s factors.json")
    if ((fj.get("files") or {}).get(name) or {}).get("sha256") != sha:
        raise Refused(f"{key}: factors/factors.json records another sha256: refused")
    man, man_sha = _json(td / f"{name}_manifest.json", f"{key}'s manifest")
    if man.get("name", name) != name or ((man.get("files") or {}).get(name) or {}).get("sha256") != sha:
        raise Refused(f"{key}: its manifest records another table or sha256: refused")
    try:
        df = pd.read_csv(io.BytesIO(raw), float_precision="round_trip")
    except Exception:  # noqa: BLE001
        raise Refused(f"{key}: not a readable CSV: refused") from None
    if list(df.columns) != ["gene", "factor"]:
        raise Refused(f"{key}: columns {list(df.columns)}, not ['gene', 'factor']: refused")
    got = list(df["gene"].astype(str))
    if len(set(got)) != len(got):
        raise Refused(f"{key}: duplicated genes: refused")
    if got != list(genes):
        raise Refused(f"{key}: genes differ from the barred factor's genes (order included): refused")
    try:
        T = df["factor"].to_numpy(dtype=np.float64)
    except (TypeError, ValueError):
        raise Refused(f"{key}: factor values are not numbers: refused") from None
    if not np.all(np.isfinite(T) & (T > 0)):
        raise Refused(f"{key}: factor values must be finite and positive: refused")
    rec = {"key": key, "name": name, "file": f"<{var}>/{name}/factors/{name}.csv", "dir_realpath": str(d),
           "sha256": sha, "listing_sha256": listing_sha256, "pinned_sha256": pins.get(key), "combine_file": cfile,
           "combine_sha256": comb_sha, "manifest_sha256": man_sha}
    rec |= {f: man.get(f) for f in RECORD_FIELDS}
    rec |= {f: man.get(f) for f in OPTIONAL_FIELDS}  # (b4)'s n_g; absent means the fallback n_g = U, never (b5)
    rec["studies"] = [str(s) for s in rec["studies"] or []]
    rec["tiers"] = sorted({str(t) for t in rec["tiers"] or []})
    return T, rec, comb


def load_tables(genes: list[str], path: Path | None = None, pin: str | None = ..., pins: dict | None = None,
                env: dict | None = None, check_tracked: bool = True) -> dict:
    """The listing checked against its pins, then every table loaded (any failure refuses the whole run).
    Returns {'listing_sha256', 'tables': {key: (T, record)}, 'combine': {key: {name: sha}}}."""
    pins = K24_TABLE_SHA256 if pins is None else pins
    listing, lsha = read_listing(path, pin, pins, check_tracked)
    tables, comb = {}, {}
    for k in KEYS:
        T, rec, c = load_k24_table(k, genes, listing, lsha, pins, env)
        tables[k], comb[k] = (T, rec), c
    return {"listing_sha256": lsha, "tables": tables, "combine": comb}


def ext_record(loaded: dict) -> dict:
    """What every task file records of the tables (the VOID check compares it across files)."""
    return {"listing_sha256": loaded["listing_sha256"],
            "tables": {k: r["sha256"] for k, (_, r) in loaded["tables"].items()},
            "combine_sha256": {k: r["combine_sha256"] for k, (_, r) in loaded["tables"].items()}}


def pinned_ext() -> dict:
    return {"listing_sha256": K24_LISTING_SHA256, "tables": dict(K24_TABLE_SHA256)}


# ---------------------------------------------------------------------------------------------------------------------
# The size rule and the scaled tables (one function, shared with build_k24_form.py)


def exponent(T) -> float | None:
    """e_T = S* / SD(ln T), population SD over all genes in file order (genes exactly 1 included); None when SD is 0.
    A table of one constant value has SD exactly 0 by the rule, whatever rounding np.std returns for it (the second
    review's blocker: all 4.0 gave 2.2e-16 and an exponent near 1e15), so a zero range of ln T also gives None."""
    lt = np.log(np.asarray(T, dtype=np.float64))
    sd = float(np.std(lt))
    return S_STAR / sd if sd > 0 and float(np.ptp(lt)) > 0 else None


def joint_c(sc3, bulk, e_sc3: float, e_bulk: float) -> float | None:
    """c = S* / SD(e_sc3 ln sc3 + e_bulk ln bulk); None when that SD is 0."""
    v = e_sc3 * np.log(np.asarray(sc3, dtype=np.float64)) + e_bulk * np.log(np.asarray(bulk, dtype=np.float64))
    sd = float(np.std(v))
    return S_STAR / sd if sd > 0 and float(np.ptp(v)) > 0 else None   # a constant sum has SD exactly 0, as exponent()


def k24_table(kind: str, tables: dict, consts: dict) -> np.ndarray:
    """The scaled table of a kind (sc3, bulk or joint) from raw tables in one gene order: sc3 ** e_sc3, bulk ** e_bulk
    or sc3 ** (c e_sc3) * bulk ** (c e_bulk); nothing rescaled or clipped."""
    if kind == "sc3":
        return np.asarray(tables["sc3"], dtype=np.float64) ** float(consts["e_sc3"])
    if kind == "bulk":
        return np.asarray(tables["bulk"], dtype=np.float64) ** float(consts["e_bulk"])
    if kind == "joint":
        c = float(consts["c"])
        return (np.asarray(tables["sc3"], dtype=np.float64) ** (c * float(consts["e_sc3"]))
                * np.asarray(tables["bulk"], dtype=np.float64) ** (c * float(consts["e_bulk"])))
    raise RuntimeError(f"unknown kind {kind!r}")


# ---------------------------------------------------------------------------------------------------------------------
# The pre-gate


def _close(a, b, rel: float = REL_TOL) -> bool:
    if isinstance(a, (float, int)) and isinstance(b, (float, int)) and not isinstance(a, bool) \
            and not isinstance(b, bool) and (isinstance(a, float) or isinstance(b, float)):
        return math.isclose(float(a), float(b), rel_tol=rel, abs_tol=0.0) or a == b
    return a == b


def _same(got, want) -> bool:
    """Equal figures, floats to a relative REL_TOL, nested dicts key for key."""
    if isinstance(want, dict):
        return isinstance(got, dict) and set(got) == set(want) and all(_same(got[k], want[k]) for k in want)
    return _close(got, want)


def check_two(phi_m, two: tuple, bases_by_world: dict, den, genes, k23_pre: dict) -> dict:
    """(i) `two` at e = 2 through K23's own pre-gate (k23.pregate_exponents): S* and the table, product, controls,
    twin, a and k must equal K23's committed pregate.json figures to a relative 1e-12, else Refused."""
    T, rec = two[0], two[1]
    loaded = {"tables": {"all": (np.asarray(T, dtype=np.float64), rec)}, "combine": two[2] if len(two) > 2 else {}}
    p = k23.pregate_exponents(phi_m, loaded, bases_by_world, den, genes)
    got, want = (p.get("exponents") or {}).get("2"), ((k23_pre or {}).get("exponents") or {}).get("2")
    if got is None or want is None:
        raise Refused("`two` at e = 2: K23's pre-gate figures are missing (a drop, or K23's committed pregate.json "
                      "lacks them): refused")
    sd = float(np.std(np.log(np.asarray(T, dtype=np.float64) ** 2.0)))
    if not _close(sd, S_STAR) or not _close(float(want["table"]["sd_ln"]), S_STAR):
        raise Refused(f"S* recomputed from `two` is {sd!r}, not {S_STAR!r} (relative {REL_TOL}): refused")
    bad = [k for k in ("table", "product", "controls", "size", "a", "k", "den_shift", "sd_ln_phi_m", "sd_ln_product")
           if not _same(got.get(k), want.get(k))]
    if bad:
        raise Refused(f"`two` at e = 2 differs from K23's committed pregate.json in {bad} (relative {REL_TOL}): "
                      "refused")
    return got


def arm_figures(phi_m, scaled, bases_by_world: dict, den, genes=None) -> tuple[dict, list[str]]:
    """One try's target-free figures and K23's product bounds: the product phi_M * scaled, each world's control
    (phi_M * perm(scaled), perm in the world's gene order and basis), the twin k * phi_M ** a (K23's formula)."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    sd_phi = float(np.std(np.log(phi_m)))
    den_phi = float(den @ phi_m / den.sum())
    prod = ext.product(phi_m, 1.0, scaled)
    sd_prod = float(np.std(np.log(prod)))
    a = sd_prod / sd_phi
    k = float((den @ prod) / (den @ (phi_m ** a)))
    twin = k * phi_m ** a
    pf = ext.figures(phi_m, prod, den)
    controls = {}
    for w in sorted(bases_by_world):
        wd = bases_by_world[w]
        phi_w = k23._world_vectors(phi_m, genes, wd)
        s_w = k23._world_vectors(scaled, genes, wd)
        den_w = k23._world_vectors(den, genes, wd)
        controls[w] = ext.figures(phi_w, ext.product(phi_w, 1.0, ext.permuted(s_w, wd["basis"])), den_w)
    shift = pf["den_total"] / den_phi - 1.0
    why = []
    if pf["genes_at_clip"] > MAX_AT_CLIP:
        why.append(f"the product has {pf['genes_at_clip']} genes at or beyond the clip (> {MAX_AT_CLIP})")
    for w, c in controls.items():
        if c["genes_at_clip"] > MAX_AT_CLIP:
            why.append(f"the control in world {w} has {c['genes_at_clip']} genes at or beyond the clip "
                       f"(> {MAX_AT_CLIP})")
    if abs(shift) > MAX_DEN_SHIFT:
        why.append(f"the product's den-weighted total moves {shift:+.2%} (beyond {MAX_DEN_SHIFT:.0%})")
    if sd_prod > MAX_SD_RATIO * sd_phi:
        why.append(f"SD(ln product) {sd_prod:.4f} exceeds {MAX_SD_RATIO} x SD(ln phi_M) {sd_phi:.4f}")
    if pf["genes_changed"] < MIN_CHANGED:
        why.append(f"the product changes {pf['genes_changed']} genes (< {MIN_CHANGED})")
    figs = {"product": pf, "controls": controls, "size": ext.figures(phi_m, twin, den), "a": float(a), "k": float(k),
            "den_shift": float(shift), "sd_ln_phi_m": sd_phi, "sd_ln_product": sd_prod}
    return figs, why


def chance_agreement(n: int) -> Fraction:
    """2 P(Bin(n, 1/2) >= ceil(2n/3)): a gene's chance of a line agreement of 2/3 or more over n lines, an exact
    fraction (math.comb, integers only)."""
    n = int(n)
    if n < 1:
        raise ValueError(f"chance_agreement needs n >= 1, not {n}")
    j0 = -(-2 * n // 3)  # ceil(2n/3) in integers
    return Fraction(2 * sum(math.comb(n, j) for j in range(j0, n + 1)), 2 ** n)


def lines_held(rec: dict, size: int, U: int, share: list) -> tuple[list[int] | None, str, str | None]:
    """(n_g per gene in the table's gene order, its source, the reason it is malformed or None). n_g comes from the
    manifest's lines_held_by_study[MGP-DMDD] when it is there: one integer per gene, 0 <= n_g <= U, and n_g >= 1
    wherever the gene has an agreement value; anything else, a field that records nothing for MGP-DMDD included, is
    malformed (the caller drops the bulk table). Only when the manifest records no lines_held_by_study at all (a table
    built before combine.py recorded it) is n_g = U for every gene ('fallback U')."""
    field = rec.get("lines_held_by_study")
    if field is None:
        return [U] * size, "fallback U (the manifest records no lines_held_by_study)", None
    if not isinstance(field, dict):
        return None, "malformed", "lines_held_by_study is not a mapping of study to counts"
    if BULK_STUDY not in field:
        # combine.py writes this field and line_agreement_by_study under one condition, and (b4) reaches n_g only when
        # an agreement value exists for the study: a field that then records nothing for it contradicts the manifest
        return None, "malformed", (f"lines_held_by_study records nothing for {BULK_STUDY} although "
                                   "line_agreement_by_study does")
    v = field[BULK_STUDY]
    if not isinstance(v, list) or len(v) != size:
        return None, "malformed", f"lines_held_by_study[{BULK_STUDY}] is not one count per gene of the table"
    if any(isinstance(x, bool) or not isinstance(x, int) for x in v):
        return None, "malformed", f"lines_held_by_study[{BULK_STUDY}] holds a value that is not an integer"
    if any(x < 0 or x > U for x in v):
        return None, "malformed", f"lines_held_by_study[{BULK_STUDY}] holds a count outside [0, {U}]"
    if any(share[i] is not None and v[i] < 1 for i in range(size)):
        return None, "malformed", (f"lines_held_by_study[{BULK_STUDY}] holds 0 lines for a gene with an agreement "
                                   "value")
    return list(v), f"manifest lines_held_by_study[{BULK_STUDY}]", None


def somite_report(guards: dict | None) -> dict:
    """Reported only: lines whose knockout median somite count is below the wild type's, and the median differences,
    where a unit's guard record carries G6's somite record; G12's shares where recorded."""
    below, diffs, shares, recorded = 0, {}, {}, 0
    for unit, g in (guards or {}).items():
        stages = g if isinstance(g, dict) else {}
        for st, v in stages.items():
            if not isinstance(v, dict):
                continue
            som = v.get("somites") or {}
            if som.get("recorded"):
                recorded += 1
                d = som.get("median_difference")
                diffs[f"{unit} {st}"] = d
                below += int(d is not None and d < 0)
            if v.get("composition") is not None:
                shares[f"{unit} {st}"] = v.get("composition")
    return {"units_with_somites": recorded, "lines_ko_below_wt": below, "median_differences": diffs,
            "composition": shares,
            "note": "reported only" + ("" if recorded or shares else "; the manifest's guards_by_unit carries no "
                                       "somite or composition record")}


def bulk_drops(T, rec: dict) -> tuple[list[str], dict]:
    """(b1) to (b5) of a bulk table, target-free; returns the reasons and the figures."""
    why, figs = [], {}
    missing = [f for f in RECORD_FIELDS if rec.get(f) is None]
    if missing:
        why.append(f"(b5) the manifest lacks {missing}")
    studies = list(rec.get("studies") or [])
    if rec.get("modality") != BULK_MODALITY or list(rec.get("tiers") or []) != [BULK_TIER]:
        why.append(f"(b1) modality {rec.get('modality')!r}, tiers {rec.get('tiers')}, not {BULK_MODALITY!r} and "
                   f"[{BULK_TIER!r}]")
    # a bulk table of two or more studies is dropped: row K24 reads the one-study MGP-DMDD table only
    if len(studies) != 1 or studies != [BULK_STUDY] or rec.get("sign_rule") != BULK_SIGN_RULE:
        why.append(f"(b1) studies {studies}, sign rule {rec.get('sign_rule')!r}, not exactly [{BULK_STUDY!r}] with "
                   f"{BULK_SIGN_RULE!r} (a bulk table of two or more studies is dropped)")
    gbs = rec.get("genes_by_study") or {}
    lines = gbs.get(BULK_STUDY) if isinstance(gbs, dict) else None
    U = len(lines) if isinstance(lines, list) else 0
    figs["lines"] = U
    if U < BULK_MIN_LINES:
        why.append(f"(b2) {U} lines in {BULK_STUDY} (< {BULK_MIN_LINES})")
    libs = rec.get("libraries_by_unit") or {}
    n_ko = n_ctl = 0
    short = []
    if isinstance(libs, dict) and libs:
        for unit, v in libs.items():
            ko = (v or {}).get("ko") if isinstance(v, dict) else None
            ctl = (v or {}).get("control") if isinstance(v, dict) else None
            nk = len(ko) if isinstance(ko, list) else 0
            nc = len(ctl) if isinstance(ctl, list) else 0
            n_ko, n_ctl = n_ko + nk, n_ctl + nc
            if nk < BULK_MIN_LIBS_UNIT or nc < BULK_MIN_LIBS_UNIT:
                short.append(unit)
    else:
        short.append("(none recorded)")
    figs |= {"libraries_ko": n_ko, "libraries_wt": n_ctl, "units_short": short}
    if short:
        why.append(f"(b3) units with under {BULK_MIN_LIBS_UNIT} libraries in an arm: {short[:5]}")
    if n_ko < BULK_MIN_LIBS_TOTAL or n_ctl < BULK_MIN_LIBS_TOTAL:
        why.append(f"(b3) {n_ko} knockout and {n_ctl} wild-type libraries (each must be >= {BULK_MIN_LIBS_TOTAL})")
    T = np.asarray(T, dtype=np.float64)
    las = rec.get("line_agreement_by_study") or {}
    share = las.get(BULK_STUDY) if isinstance(las, dict) else None
    empty = {"agree_G": None, "agree_count": None, "E0": None, "agree_needed": None, "n_g_source": None,
             "n_g_histogram": None, "agree_left_out_few_lines": None}
    if not isinstance(share, list) or len(share) != T.size:
        why.append(f"(b4) no line agreement for {BULK_STUDY} over the table's {T.size} genes")
        figs |= empty
    else:
        n_g, source, bad = lines_held(rec, T.size, U, share)
        if bad is not None:
            why.append(f"(b4) {bad}: malformed")
            figs |= empty | {"n_g_source": source}
        else:
            ne1 = T != 1.0
            valued = [i for i in range(T.size) if ne1[i] and share[i] is not None]
            held = [i for i in valued if n_g[i] >= AGREE_MIN_LINES]
            G = len(held)
            n_agree = sum(1 for i in held if float(share[i]) >= AGREE - AGREE_TOL)
            E0 = sum((chance_agreement(n_g[i]) for i in held), Fraction(0))
            need = max(Fraction(AGREE_FLOOR), E0 + Fraction(str(AGREE_EXCESS)) * G)
            hist: dict[str, int] = {}
            for i in valued:
                hist[str(n_g[i])] = hist.get(str(n_g[i]), 0) + 1
            figs |= {"agree_G": G, "agree_count": n_agree, "E0": float(E0), "agree_needed": float(need),
                     "n_g_source": source, "n_g_histogram": dict(sorted(hist.items(), key=lambda kv: int(kv[0]))),
                     "agree_left_out_few_lines": len(valued) - G}
            if n_agree < need:
                why.append(f"(b4) {n_agree} of {G} genes != 1 held in >= {AGREE_MIN_LINES} lines with line agreement "
                           f">= 2/3 (< max({AGREE_FLOOR}, E0 {float(E0):.3f} + {AGREE_EXCESS} x {G}) = "
                           f"{float(need):.3f})")
    figs["somites_and_composition"] = somite_report(rec.get("guards_by_unit"))
    return why, figs


def _record_figs(rec: dict) -> dict:
    """The manifest's figures the pre-gate records (no gene-level value)."""
    gbs = rec.get("genes_by_study") or {}
    ubs = rec.get("units_by_study") or {}
    return {"studies": rec.get("studies"), "tiers": rec.get("tiers"), "modality": rec.get("modality"),
            "sign_rule": rec.get("sign_rule"), "refstats_sha256": rec.get("refstats_sha256"),
            "sha256": rec.get("sha256"), "checked": rec.get("checked"),
            "lines_by_study": {s: len(v) for s, v in gbs.items()} if isinstance(gbs, dict) else None,
            "units_by_study": {s: len(v) for s, v in ubs.items()} if isinstance(ubs, dict) else None,
            "libraries_by_unit": rec.get("libraries_by_unit"), "k_by_unit": rec.get("k_by_unit")}


def pregate_k24(phi_m, loaded: dict, bases_by_world: dict, den, genes: list[str] | None = None,
                k23_pre: dict | None = None) -> dict:
    """The pre-gate of row K24, target-free, from phi_M, the loaded tables, each world's basis and WT E8.75's mean
    within-cell proportions (den) alone; a pure function. Order: (i) `two` at e = 2 against K23's committed figures
    (when k23_pre is given; a mismatch raises Refused); (ii) e_T; (iii) the per-table drops; (iv) the joint; (v) per
    arm the product, the world controls and the twin."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    tabs = loaded["tables"]
    comb = loaded.get("combine") or {}
    two_T, two_rec = tabs["two"]
    out: dict = {"n_genes": int(phi_m.size), "s_star": S_STAR, "two": None, "two_checked": k23_pre is not None,
                 "phi_m": ext.figures(phi_m, phi_m, den), "tables": {}, "scaled": {}, "records": {}, "e": {},
                 "c": None, "corr": {}, "arms": {}, "bulk": None, "dropped": {}, "sc3_repeat_only": False,
                 "constants": {"s_star": S_STAR, "e_max": E_MAX, "e_min": E_MIN, "min_genes": MIN_GENES,
                               "max_abs_ln_scaled": MAX_ABS_LN_SCALED, "joint_corr": list(JOINT_CORR),
                               "sc3_two_max_corr": SC3_TWO_MAX_CORR, "bulk_min_lines": BULK_MIN_LINES,
                               "bulk_min_libs_unit": BULK_MIN_LIBS_UNIT, "bulk_min_libs_total": BULK_MIN_LIBS_TOTAL,
                               "agree": AGREE, "agree_floor": AGREE_FLOOR, "agree_excess": AGREE_EXCESS,
                               "agree_min_lines": AGREE_MIN_LINES,
                               "bulk_sign_rule": BULK_SIGN_RULE, "bulk_study": BULK_STUDY,
                               "min_changed": MIN_CHANGED, "max_abs_corr": MAX_ABS_CORR, "min_studies": MIN_STUDIES,
                               "max_at_clip": MAX_AT_CLIP, "max_den_shift": MAX_DEN_SHIFT,
                               "max_sd_ratio": MAX_SD_RATIO, "perm_seed": ext.PERM_SEED, "perm_bins": ext.PERM_BINS,
                               "worlds": sorted(bases_by_world)}}
    # (i)
    if k23_pre is not None:
        out["two"] = check_two(phi_m, (two_T, two_rec, comb.get("two") or {}), bases_by_world, den, genes, k23_pre)
    out["tables"]["two"] = ext.table_figures(phi_m, two_T, den)
    out["records"]["two"] = _record_figs(two_rec)
    # (ii) and (iii)
    scaled: dict = {}
    for key in ("sc3", "bulk"):
        T, rec = tabs[key]
        T = np.asarray(T, dtype=np.float64)
        tf = ext.table_figures(phi_m, T, den)
        out["tables"][key] = tf
        out["records"][key] = _record_figs(rec)
        e = exponent(T)
        out["e"][key] = e
        why = []
        repeat = None  # sc3 only: the repeat clause's reason, kept apart from every other drop
        if e is None:
            why.append("SD(ln T) is 0: no exponent")
        elif e > E_MAX or e < E_MIN:
            # (s1) first (the third review's blocker): an exponent outside [E_MIN, E_MAX] drops the table BEFORE
            # T ** e_T is formed; for a nearly constant table the exponent is huge, T ** e_T overflows to inf or
            # underflows to 0, and table_figures would raise where the rule says the table drops. The scaled table's
            # figures and its arm figures are then not made; (s3) is computed from e_T and max |ln T| directly.
            why.append(f"(s1) e_T {e:.4f} outside [{E_MIN}, {E_MAX}]")
            mx = float(e * np.max(np.abs(np.log(T))))
            if mx > MAX_ABS_LN_SCALED:
                why.append(f"(s3) max |ln T ** e_T| {mx:.4f} (> {MAX_ABS_LN_SCALED})")
        else:
            Te = T ** e
            sf = ext.table_figures(phi_m, Te, den)
            sf["max_abs_ln"] = float(np.max(np.abs(np.log(Te))))
            out["scaled"][key] = sf
            if not _close(sf["sd_ln"], S_STAR):
                raise Refused(f"{key}: SD(ln T ** e_T) {sf['sd_ln']!r} is not S* (relative {REL_TOL}): refused")
            if sf["max_abs_ln"] > MAX_ABS_LN_SCALED:
                why.append(f"(s3) max |ln T ** e_T| {sf['max_abs_ln']:.4f} (> {MAX_ABS_LN_SCALED})")
            scaled[key] = Te
        if tf["genes_ne_1"] < MIN_GENES:
            why.append(f"(s2) {tf['genes_ne_1']} genes != 1 (< {MIN_GENES})")
        c0 = tf["corr_ln_phi_m"]
        if c0 is None or abs(c0) >= MAX_ABS_CORR:
            why.append(f"|corr(ln phi_M, ln T)| {'undefined' if c0 is None else f'{abs(c0):.3f}'} "
                       f"(>= {MAX_ABS_CORR})")
        tiers = list(rec.get("tiers") or [])
        if tiers and set(tiers) <= {ext.CARDIAC}:
            why.append("only cardiac-tier studies")
        if "wnt" in tiers:
            why.append("a wnt-tier study")
        if rec.get("refstats_sha256") is None or rec.get("refstats_sha256") != two_rec.get("refstats_sha256"):
            why.append("its reference statistics are not `two`'s")
        if key == "sc3":
            if len(rec.get("studies") or []) < MIN_STUDIES:
                why.append(f"its manifest lists {len(rec.get('studies') or [])} study (< {MIN_STUDIES})")
            if rec.get("sign_rule") != "applied":
                why.append(f"sign rule {rec.get('sign_rule')!r}, not 'applied'")
            same = sorted(n for n, v in (comb.get("sc3") or {}).items() if v == rec["sha256"] and n != "all"
                          and (n.startswith("study_") or "__loso_" in n or n in ext.SINGLE_TIERS))
            if same:
                why.append(f"its bytes equal {same}")
            if rec["sha256"] == two_rec["sha256"]:
                why.append("its bytes equal `two`'s")
            ct = ext.corr_between(T, two_T)
            out["corr"]["sc3_two"] = ct
            # restated clause (2026-11-20): an undefined corr is a drop of sc3 like any other (it takes the joint); a
            # corr of 0.9 or more is the repeat clause, which drops x_sc3, its control and its twin only
            if ct is None:
                why.append("corr(ln sc3, ln two) undefined")
            elif ct >= SC3_TWO_MAX_CORR:
                repeat = f"corr(ln sc3, ln two) {ct:.3f} (>= {SC3_TWO_MAX_CORR}: a repeat of the named form)"
        else:
            bwhy, bfigs = bulk_drops(T, rec)
            out["bulk"] = bfigs
            why += bwhy
        if key in scaled:
            figs, pwhy = arm_figures(phi_m, scaled[key], bases_by_world, den, genes)
            out["arms"][f"x_{key}"] = figs
            why += pwhy
        if repeat is not None:
            out["sc3_repeat_only"] = not why
            why.append(repeat)  # kept in the record of sc3's drop, for reporting
        if why:
            out["dropped"][key] = "; ".join(why)
    # (iv) the joint: dropped when bulk is dropped, or sc3 for any reason other than the repeat clause alone
    jwhy = []
    if "bulk" in out["dropped"]:
        jwhy.append("bulk is dropped")
    if "sc3" in out["dropped"] and not out["sc3_repeat_only"]:
        jwhy.insert(0, "sc3 is dropped")
    out["corr"]["sc3_bulk"] = cj = ext.corr_between(tabs["sc3"][0], tabs["bulk"][0])
    if cj is None or cj < JOINT_CORR[0] or cj >= JOINT_CORR[1]:
        jwhy.append(f"corr(ln sc3, ln bulk) {'undefined' if cj is None else f'{cj:.3f}'} outside "
                    f"[{JOINT_CORR[0]}, {JOINT_CORR[1]})")
    if out["e"].get("sc3") is not None and out["e"].get("bulk") is not None:
        c = joint_c(tabs["sc3"][0], tabs["bulk"][0], out["e"]["sc3"], out["e"]["bulk"])
        out["c"] = c
        if c is None:
            jwhy.append("SD(e_sc3 ln sc3 + e_bulk ln bulk) is 0: no c")
        elif not jwhy:
            J = k24_table("joint", {"sc3": tabs["sc3"][0], "bulk": tabs["bulk"][0]},
                          {"e_sc3": out["e"]["sc3"], "e_bulk": out["e"]["bulk"], "c": c})
            sf = ext.table_figures(phi_m, J, den)
            sf["max_abs_ln"] = float(np.max(np.abs(np.log(J))))
            out["scaled"]["joint"] = sf
            if not _close(sf["sd_ln"], S_STAR):
                raise Refused(f"joint: SD(ln J) {sf['sd_ln']!r} is not S* (relative {REL_TOL}): refused")
            figs, pwhy = arm_figures(phi_m, J, bases_by_world, den, genes)
            out["arms"]["x_joint"] = figs
            jwhy += pwhy
    else:
        jwhy.append("no exponent for a table")
    if jwhy:
        out["dropped"]["joint"] = "; ".join(jwhy)
    return out


def dropped_arms(pre: dict) -> dict[str, str]:
    """The pre-gate's drops as arm drops: a dropped table drops its try, its control and its twin, and the joint with
    them, except that sc3 dropped by the repeat clause alone (pre['sc3_repeat_only']) leaves the joint in."""
    out = {}
    sc3_takes_joint = "sc3" in pre["dropped"] and not pre.get("sc3_repeat_only")
    for kind in ("sc3", "bulk", "joint"):
        why = pre["dropped"].get(kind)
        if why is None and kind == "joint" and (sc3_takes_joint or "bulk" in pre["dropped"]):
            why = "a table is dropped"
        if why is not None:
            for a in ARMS_OF_KIND[kind]:
                out[a] = f"{kind}: {why}"
    return out


def active_arms(pre: dict) -> tuple[str, ...]:
    d = dropped_arms(pre)
    return tuple(a for a in ARMS if a not in d)


def nothing_runs(arms: tuple[str, ...]) -> bool:
    return not any(t in arms for t in NAMEABLE)


def constants_of(pre: dict) -> dict:
    """e_sc3, e_bulk, c and each undropped try's twin constants (a, k), at full precision. e_sc3 is carried whenever
    the pre-gate made it, so also when x_sc3 is dropped as a repeat of the named form and the joint runs on it."""
    arms = active_arms(pre)
    return {"e_sc3": pre["e"].get("sc3"), "e_bulk": pre["e"].get("bulk"),
            "c": pre.get("c") if "x_joint" in arms else None,
            "twins": {t: {"a": pre["arms"][t]["a"], "k": pre["arms"][t]["k"]} for t in NAMEABLE if t in arms}}


def try_count(arms) -> int:
    n = sum(a not in GUARDS for a in arms)
    if n > MAX_TRIES:
        raise RuntimeError(f"{n} tries: never more than {MAX_TRIES}")
    return n


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def fold_factors(bs, ref, per: dict, phi: np.ndarray, group: np.ndarray, basis: np.ndarray, phi_barred: np.ndarray,
                 spec: dict) -> dict:
    """K22's body over this try's loaded tables: K18's fold_factors (the original, called unchanged) plus each table in
    the reference's gene order and its permutation within deciles of `basis`; then the scaled tables of the kinds the
    pre-gate kept and their permutations (one per world: perm of the scaled vector). Each kind needs only its own
    constants (sc3: e_sc3; bulk: e_bulk; joint: e_sc3, e_bulk and c), so a table with no exponent leaves the other
    table's arms whole."""
    if not _K24.get("T"):
        raise RuntimeError("the K24 tables are not loaded (k24_wider_tables.py loads them before any arm)")
    F = k22._K18_ORIGINAL["fold_factors"](bs, ref, per, phi, group, basis, phi_barred, spec)
    F["ext"] = {n: k18.in_gene_order(T, _K24["genes"], list(ref.genes)) for n, T in _K24["T"].items()}
    F["ext_perm"] = {n: ext.permuted(T, basis) for n, T in F["ext"].items()}
    consts = _K24.get("consts")
    F["scaled"], F["scaled_perm"] = {}, {}
    if consts is not None:
        for kind in ("sc3", "bulk", "joint"):
            need = ("e_sc3", "e_bulk", "c") if kind == "joint" else (f"e_{kind}",)
            if any(consts.get(k) is None for k in need):
                continue
            F["scaled"][kind] = k24_table(kind, F["ext"], consts)
            F["scaled_perm"][kind] = ext.permuted(F["scaled"][kind], basis)
    return F


def arm_factor(arm: str, F: dict) -> np.ndarray:
    """The factor an arm passes in place of phi, in the reference's gene order."""
    if arm == BASE:
        return F["phi_barred"]
    if arm == REPRO:
        return ext.product(F["phi_barred"], 1.0, F["ext"]["two"], 2.0)
    if arm not in KIND_OF:
        raise RuntimeError(f"unknown arm {arm!r}")
    kind = KIND_OF[arm]
    if arm.endswith("_size"):
        tw = ((_K24.get("consts") or {}).get("twins") or {}).get(f"x_{kind}")
        if tw is None:
            raise RuntimeError(f"the twin's constants for x_{kind} are not loaded from pregate.json")
        return tw["k"] * F["phi_barred"] ** tw["a"]
    src = F["scaled_perm"] if arm.endswith("_perm") else F["scaled"]
    if kind not in src:
        raise RuntimeError(f"the scaled table of {kind} is not made (its constants are not loaded)")
    return ext.product(F["phi_barred"], 1.0, src[kind], 1.0)


def build_arm(bs, ref, arm: str, n: int, seed: int, F: dict):
    """One arm's drawn cells, as K22's and K23's build_arm: every arm draws from stable_rng(seed, STREAM) through
    bs.construct with the barred spec's keywords; only the factor passed differs. barred_form passes them unchanged."""
    if arm not in ARMS:
        raise RuntimeError(f"unknown arm {arm!r}")
    rng = bs.t3d.stable_rng(seed, STREAM)
    kw, _ = bs.prediction_inputs(ref, F["barred_spec"], {}, F["phi_barred"])
    if kw["ratio"] is not None or kw["alphas"] is not None or kw["w"] != 1.0:
        raise RuntimeError("the barred s1 is the Mab21l2 factor at exponent 1 with no depth ratio or stage mixing")
    if arm != BASE:
        kw["phi"] = arm_factor(arm, F)
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
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: k18._build_plain(ctx, arm),
                                 params={"form": arm, "stream": STREAM, "composition": k8.spec()["composition"],
                                         "factor": ("the barred s1's" if arm == BASE else f"k24 ({arm})"),
                                         "perm_seed": ext.PERM_SEED if arm.endswith("_perm") else None},
                                 requires=lambda t: t.key in k18.KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in arms)]


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
         "not started, nothing written. Row K24: commit or restore first; INCOMPLETE, not VOID")


def preconditions() -> tuple[str, dict]:
    """Before anything loads: the row, a clean tree, the listing tracked and at its pins (all filled; `two`'s K22's),
    the frozen files, K23's committed pregate.json and reading.json tracked, both folder variables valid."""
    rid = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse(DIRTY)
    try:
        read_listing()
        for var in (ENV_TABLES, ENV_TWO):
            folder(var)
    except Refused as e:
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
    """K23's per-seed x_e2 rows per cell, from K23's task files, each at the sha256 K23's committed reading.json
    records (else Refused). Reads the scores of x_e2 and its severity only."""
    want = k23_read.get("files_sha256") or {}
    out: dict = {}
    for w in WORLDS:
        for n in TASKS:
            rows = []
            for s in SEEDS:
                f = f"plain_{n}/seed_{s}.json" if w == "plain" else f"dm/{n}_seed_{s}.json"
                p = Path(runs) / f
                if not p.is_file() or want.get(f) is None:
                    raise Refused(f"K23's task file {f} is missing (K23's run folder {Path(runs).name}): guard 2's "
                                  "per-seed check cannot be made: refused")
                raw = p.read_bytes()
                if sha256_bytes(raw) != want[f]:
                    raise Refused(f"K23's task file {f} is not at the sha256 K23's reading records: refused")
                r = json.loads(raw)
                row = r["rows"][0] if w == "plain" else r
                rows.append({"seed": s, "scores": {REPRO: row["scores"][REPRO]},
                             "t3": {REPRO: row["t3"][REPRO]}})
            out[(w, n)] = rows
    return out


def guard_refs(k23_dir: Path | None = None, refs: Path = K18_TASKS) -> dict:
    """What the guards compare with: K18's committed barred_form rows (at their pins), K23's committed x_e2 cells and
    its per-seed x_e2 rows from K23's run folder in this checkout (k23_runs_default; no option overrides it). Refused
    when any is missing or off its sha256."""
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
# Running


def progress_log(part: str, world: str, task: str):
    """The harness's per-seed line lists every arm's score; this prints the part (guard only), world, task and seed."""
    lead = "K24 guard" if part == "guard" else "K24"

    def log(line: str) -> None:
        m = re.match(r".* seed (\d+):", str(line))
        print(f"{lead} {world} {task} seed {m.group(1) if m else '?'}: scored", flush=True)

    return log


class _Mismatch(RuntimeError):
    pass


def load_k24(genes: list[str], phi_m) -> dict:
    """The tables against the barred factor's file genes (order included), held in _K24 for fold_factors. A refusal
    exits 2."""
    try:
        loaded = load_tables(list(genes))
    except Refused as e:
        refuse(str(e))
    _K24.clear()
    _K24.update(genes=list(genes), T={k: T for k, (T, _) in loaded["tables"].items()}, loaded=loaded,
                record=ext_record(loaded), table_records={k: r for k, (_, r) in loaded["tables"].items()})
    return loaded


def load_pregate(out: Path, record: dict) -> tuple[dict, tuple[str, ...]]:
    """The run folder's pregate.json: its tables must be the loaded ones; fills _K24's constants. Refuses (exit 2)
    when it is missing or disagrees, or when it drops both tables."""
    p = Path(out) / "pregate.json"
    if not p.is_file():
        refuse(f"{p} is missing: the pre-gate runs first (`pregate --out {out}`)")
    body = json.loads(p.read_text())
    if body.get("ext") != record:
        refuse("pregate.json records other tables than the ones loaded: refused")
    pre = body["pregate"]
    arms = active_arms(pre)
    if list(arms) != body.get("arms"):
        refuse("pregate.json's arms are not its pre-gate's: refused")
    if nothing_runs(arms):
        refuse("the pre-gate leaves no try (both tables dropped, or bulk dropped and sc3 a repeat of the named form): "
               "nothing runs")
    _K24["consts"] = constants_of(pre)
    return pre, arms


def _pregate_record(pre: dict, arms: tuple[str, ...]) -> dict:
    return {"dropped": dropped_arms(pre), "arms": list(arms), "ext": dict(_K24["record"]),
            "constants": constants_of(pre)}


def _world_figures(F: dict, arms: tuple[str, ...]) -> dict:
    phi = F["phi_barred"]
    return {a: ext.figures(phi, arm_factor(a, F)) for a in arms if a != BASE}


def check_world(pre: dict, world: str, F: dict, arms: tuple[str, ...]) -> None:
    """The run's own controls equal pregate.json's controls in this world."""
    for t, c in CONTROL_OF.items():
        if c not in arms:
            continue
        want = ((pre["arms"].get(t) or {}).get("controls") or {}).get(world)
        got = ext.figures(F["phi_barred"], arm_factor(c, F))
        if want is None or any(not _close(got[k], want.get(k)) for k in got):
            raise _Mismatch(f"the {world} world's control {c} is not the pre-gate's (the pre-gate's basis is not this "
                            "world's): no task file written")


GUARD_FILE = "guard.json"
GUARD_DIR = "guard"


def guard_passed(out: Path, st: dict) -> None:
    """A scoring part starts only after a passing guard.json of this commit and machine."""
    p = Path(out) / GUARD_FILE
    if not p.is_file():
        refuse(f"{p} is missing: the guard part runs and passes first (`guard --runs {out}`)")
    g = json.loads(p.read_text())
    prov = g.get("provenance") or {}
    if g.get("passes") is not True:
        refuse("the guard part did not pass: the run is VOID and no scoring part runs")
    if prov.get("commit") != st["commit"] or prov.get("machine") != st["machine"]:
        refuse("guard.json is from another commit or machine: refused")


def run_plain(a, st: dict, out: Path, part: str) -> None:
    root = out / GUARD_DIR if part == "guard" else out
    folder_ = root / f"plain_{a.task}"
    folder_.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder_, st)
    rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
    load_k24(k18._BARRED["genes"], k18._BARRED["phi"])
    pre, arms = load_pregate(out, _K24["record"])
    run_arms = GUARDS if part == "guard" else arms
    lead = "K24 guard" if part == "guard" else "K24"
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
                check_world(pre, "plain", fold, run_arms)
            res["provenance"]["machine"] = st["machine"]
            res |= {"world": "plain", "part": part, "task_source": SCRIPT, "seed": s,
                    "pregate": _pregate_record(pre, arms), "factor_record": k7._PREP["fold"]["record"],
                    "barred_factor": rec, "world_figures": _world_figures(fold, run_arms)}
            prov = k18.write_task_file(path, res)
            print(f"wrote {path}: commit {prov['commit']}, clean {prov['tree_clean']}, machine {prov['machine']}",
                  flush=True)


def run_dm(a, st: dict, out: Path, part: str) -> None:
    """The depth-matched world, set up exactly as k23.run_dm sets it up, then each (replicate, seed) row through
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
        load_k24(genes, phi_b)
        pre, arms = load_pregate(out, _K24["record"])
        run_arms = GUARDS if part == "guard" else arms
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
                 "reference_counts_largest_gap": gap, "pregate": _pregate_record(pre, arms),
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
                check_world(pre, "dm", F, run_arms)
            wf = _world_figures(F, run_arms)
            k18._STATE.clear()
            k18._STATE.update(ctx=ctx, phi_b=phi_b, ref_counts=ref_C, ref_lab=ref_lab, ref_props=ref_P, F=F,
                              threads=a.threads)
            todo = [(ti, s) for ti, T in enumerate(ctx.targets) for s in seeds
                    if not (folder_ / f"{by_name[T.name]}_seed_{s}.json").exists()]

            def write(row: dict) -> None:
                k22._stamp_ok(st)
                path = folder_ / f"{by_name[row['target']]}_seed_{row['seed']}.json"
                k18.write_task_file(path, row | {"world": "dm", "part": part, "task": by_name[row["target"]],
                                                 "arms": list(run_arms), "pregate": _pregate_record(pre, arms),
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
    _, st = preconditions()
    out = Path(a.out)
    if a.part == "guard":
        try:
            guard_refs()
        except Refused as e:
            refuse(str(e))
    else:
        guard_passed(out, st)
    out.mkdir(parents=True, exist_ok=True)
    try:
        (run_plain if a.world == "plain" else run_dm)(a, st, out, a.part)
    except Exception:  # noqa: BLE001  (any crash: masked traceback, exit 1; seeds written so far stay)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)


def _fmt(v) -> str:
    return f"{v:.6g}" if isinstance(v, float) else str(v)


def _figs(r: dict) -> str:
    return ", ".join(f"{k} {_fmt(v)}" for k, v in r.items() if not isinstance(v, dict))


def pregate_lines(pre: dict) -> list[str]:
    lines = [f"pre-gate (the tables, phi_M and the worlds' bases alone; no score): {pre['n_genes']} genes; S* "
             f"{pre['s_star']!r}; `two` at e = 2 checked against K23's committed figures: {pre['two_checked']}"]
    for k, v in pre["tables"].items():
        lines.append(f"  table {k}: {_figs(v)}")
        lines.append(f"    record: {_figs(pre['records'].get(k) or {})}")
    for k, v in pre["scaled"].items():
        lines.append(f"  scaled {k}: {_figs(v)}")
    lines.append(f"  e_sc3 {pre['e'].get('sc3')!r}, e_bulk {pre['e'].get('bulk')!r}, c {pre['c']!r}; corr(ln sc3, ln "
                 f"two) {_fmt(pre['corr'].get('sc3_two'))}, corr(ln sc3, ln bulk) {_fmt(pre['corr'].get('sc3_bulk'))}")
    if pre.get("bulk"):
        b = pre["bulk"]
        lines.append("  bulk: " + _figs(b) + "; somites and composition (reported only): "
                     + _figs(b.get("somites_and_composition") or {}))
        lines.append(f"  bulk (b4): genes != 1 with an agreement value per n_g (counts only): {b.get('n_g_histogram')}")
    for t, v in pre["arms"].items():
        lines.append(f"  {t}: a {v['a']!r}, k {v['k']!r}, den shift {v['den_shift']:+.4%}")
        lines.append(f"    product: {_figs(v['product'])}")
        for w, c in v["controls"].items():
            lines.append(f"    control ({w} basis): {_figs(c)}")
        lines.append(f"    size twin: {_figs(v['size'])}")
    lines.append("  dropped: " + ("; ".join(f"{k} ({v})" for k, v in pre["dropped"].items()) or "none"))
    if pre.get("sc3_repeat_only"):
        lines.append("  sc3 is a repeat of the named form only: x_sc3, x_sc3_perm and x_sc3_size are dropped; sc3 "
                     "enters the joint unchanged")
    return lines


def world_bases(data: Path, genes: list[str], bs) -> dict:
    return k23.world_bases(data, genes, bs)


def cmd_pregate(a) -> None:
    """The pre-gate before any scoring: phi_M through the builder's loader (WT E8.75's gene order), the tables through
    this module's loader, each world's basis from WT E8.75 (K23's) and the den-weighted totals on its mean within-cell
    proportions; `two` checked against K23's committed pregate.json; DIR/pregate.json."""
    rid, st = preconditions()
    out = Path(a.out)
    path = out / "pregate.json"
    if path.exists():
        refuse(f"{path} exists")
    try:
        guard_refs()
        k23_pre, _ = k23_committed()
    except Refused as e:
        refuse(str(e))
    out.mkdir(parents=True, exist_ok=True)
    try:
        from vec.prep import load_stage

        bs = harness._t3_head()
        rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
        phi_m, genes = k18._BARRED["phi"], k18._BARRED["genes"]
        wt = load_stage(Path(a.data) / "T2" / "E8.75.h5ad", drop_duplicates=True)
        if list(wt.var_names) != list(genes):
            raise RuntimeError("WT E8.75's genes are not the barred factor's")
        C = np.asarray(bs.t3d.dense_counts(wt), dtype=np.float64)
        den = (C / np.maximum(C.sum(1, keepdims=True), 1)).mean(0)
        loaded = load_k24(genes, phi_m)
        bases = world_bases(Path(a.data), list(genes), bs)
        try:
            pre = pregate_k24(phi_m, loaded, bases, den, list(genes), k23_pre["pregate"])
        except Refused as e:
            refuse(str(e))
        arms = active_arms(pre)
        res = {"pregate": pre, "arms": list(arms), "dropped_arms": dropped_arms(pre), "constants": constants_of(pre),
               "tries_if_run": try_count(arms), "row": rid, "barred_factor": rec, "ext": _K24["record"],
               "table_records": _K24["table_records"], "provenance": st}
        k18.write_task_file(path, res)
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print("\n".join(pregate_lines(pre)), flush=True)
    if nothing_runs(arms):
        print("no try left (both tables dropped, or bulk dropped and sc3 a repeat of the named form): nothing runs and "
              "nothing is counted", flush=True)
        sys.exit(EXIT_NOTHING_RUNS)


# ---------------------------------------------------------------------------------------------------------------------
# The guard part


def guard_files() -> list[str]:
    return ([f"{GUARD_DIR}/plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"{GUARD_DIR}/dm/{n}_seed_{s}.json" for n in TASKS for s in SEEDS])


def _rows_of(runs: Path, files: list[str], prefix: str, want_arms: tuple[str, ...], why: list[str], sha: dict,
             prov: dict, want_rec: dict | None = None) -> dict:
    """Rows per cell from task files (plain: rows[0]; dm: the row), each checked for its task, seed and arms."""
    rows: dict = {(w, n): [] for w in WORLDS for n in TASKS}
    for f in files:
        r, sha[f] = k21._read_json(Path(runs) / f)
        prov[f] = r["provenance"]
        rel = f[len(prefix):]
        world, rest = ("plain", rel.split("/")[0][len("plain_"):]) if rel.startswith("plain_") else ("dm", rel[3:5])
        seed = int(rel.rsplit("_", 1)[1][:-len(".json")])
        if want_rec is not None:
            rec = r.get("pregate") or {}
            for k in ("arms", "dropped", "constants", "ext"):
                if rec.get(k) != want_rec[k]:
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


def _gated(rows: list[dict], arm: str) -> int:
    return k21._gated([q["t3"][arm]["severity"] for q in rows])


def guard_eval(rows: dict, refs: dict) -> dict:
    """GUARD 1: barred_form within K17's tolerance of K18's committed barred_form in all four cells. GUARD 2: x_e2
    within the same tolerance of K23's per-seed x_e2 and of its committed ten-seed mean, its gated count equal to
    K23's, in all four cells. Returns the per-cell records, the failed guard arms and passes."""
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
    """The guard check over DIR/guard/: writes DIR/guard.json; exit 0 when both guards hold, 6 when either fails (the
    run is VOID: the launcher reads it, and the reading counts the failed guard arms only), 3 when a file is missing."""
    rid, st = preconditions()
    runs = Path(a.runs)
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
        sys.exit(1)
    print(f"guard part: {'passes' if res['passes'] else 'FAILS'} (failed guard arms: {res['failed'] or 'none'})",
          flush=True)
    if not res["passes"]:
        sys.exit(EXIT_GUARD)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _cells():
    return [(w, n) for w in WORLDS for n in TASKS]


def extra_pairs(arms: tuple[str, ...]) -> list[tuple[str, str]]:
    """Beyond K21's pairs (each arm against barred_form, each try against its control): each try against its size
    twin and against x_e2."""
    out = [(t, SIZE_OF[t]) for t in NAMEABLE if t in arms and SIZE_OF[t] in arms]
    out += [(t, REPRO) for t in NAMEABLE if t in arms and REPRO in arms]
    return out


def _summarise(rows: list[dict], arms: tuple[str, ...]) -> dict:
    saved = dict(k21.CONTROL_OF)
    k21.CONTROL_OF.clear()
    k21.CONTROL_OF.update(CONTROL_OF)
    try:
        out = k21.summarise(rows, arms)
    finally:
        k21.CONTROL_OF.clear()
        k21.CONTROL_OF.update(saved)
    for x, y in extra_pairs(arms):
        out["paired"][f"{x} - {y}"] = k7._paired(rows, x, y)
    return out


def _beats(v: dict | None) -> bool:
    return bool(v is not None and v["mean"] > MIN_GAIN and v["ahead"] >= MIN_AHEAD)


def _status(per: dict, t: str, arms: tuple[str, ...]) -> dict:
    """K23's (a) to (e) unchanged and (f) the head to head against x_e2."""
    if t not in arms:
        return {"ran": False, "qualified": False, "qualified_a_to_e": False}
    d = {(w, n): per[w][n]["paired"][f"{t} - {BASE}"] for w, n in _cells()}
    ahead = {k: _beats(v) for k, v in d.items()}
    behind = {k: bool(v["mean"] < -MIN_GAIN) for k, v in d.items()}
    arm_gated = lambda k, x: per[k[0]][k[1]]["arms"][x]["severity_gated"]  # noqa: E731
    gated = {k: (arm_gated(k, t), arm_gated(k, BASE)) for k in d}
    pair = lambda k, other: per[k[0]][k[1]]["paired"].get(f"{t} - {other}")  # noqa: E731
    a = (any(ahead[("plain", n)] for n in TASKS) and any(ahead[("dm", n)] for n in TASKS)
         and any(ahead[(w, "g1")] for w in WORLDS) and any(ahead[(w, "g2")] for w in WORLDS)
         and any(ahead[("dm", n)] and gated[("dm", n)][0] == gated[("dm", n)][1] for n in TASKS))
    b = not any(behind.values())
    cc = all(g <= g0 + MAX_EXTRA_GATED for g, g0 in gated.values())
    up = [k for k in d if ahead[k]]
    dd = (all(_beats(pair(k, CONTROL_OF[t])) for k in up)
          and any(k[0] == "dm" and _beats(pair(k, CONTROL_OF[t])) for k in up))
    ee = all(_beats(pair(k, SIZE_OF[t])) for k in up)
    h2h = {k: pair(k, REPRO) for k in d}
    ff = (REPRO in arms and all(v is not None for v in h2h.values())
          and any(_beats(h2h[("plain", n)]) for n in TASKS)
          and any(_beats(h2h[("dm", n)]) and arm_gated(("dm", n), t) == arm_gated(("dm", n), REPRO) for n in TASKS)
          and not any(v["mean"] < -MIN_GAIN for v in h2h.values()))
    a_to_e = bool(a and b and cc and dd and ee)
    return {"ran": True, "qualified": bool(a_to_e and ff), "qualified_a_to_e": a_to_e, "a": bool(a), "b": bool(b),
            "c": bool(cc), "d": bool(dd), "e": bool(ee), "f": bool(ff),
            "ahead_cells": [f"{w} {n}" for (w, n), v in ahead.items() if v],
            "behind_cells": [f"{w} {n}" for (w, n), v in behind.items() if v],
            "means": {f"{w} {n}": float(v["mean"]) for (w, n), v in d.items()},
            "head_to_head": {f"{w} {n}": (None if v is None else float(v["mean"])) for (w, n), v in h2h.items()},
            "mean_of_four": float(np.mean([v["mean"] for v in d.values()])),
            "gated_vs_base": {f"{w} {n}": list(v) for (w, n), v in gated.items()},
            "report": (NOT_AHEAD if a_to_e and not ff else None)}


def decide(per: dict, arms: tuple[str, ...]) -> dict:
    """The rule fixed in the row: at most one name, only an arm that meets (a) to (f)."""
    status = {t: _status(per, t, arms) for t in NAMEABLE}
    q = [t for t in NAMEABLE if status[t]["qualified"]]
    rank = None
    if not q:
        name, rule = None, "nothing qualifies on (a) to (f)"
    elif len(q) == 1:
        name, rule = q[0], "one arm qualified"
    else:
        base_g = {k: per[k[0]][k[1]]["arms"][BASE]["severity_gated"] for k in _cells()}
        same = [k for k in _cells() if all(per[k[0]][k[1]]["arms"][x]["severity_gated"] == base_g[k]
                                           for x in (*q, REPRO))]
        cells = same or _cells()
        means = {t: float(np.mean([per[w][n]["paired"][f"{t} - {REPRO}"]["mean"] for w, n in cells])) for t in q}
        best = max(means.values())
        name = next(t for t in NAMEABLE if t in q and means[t] == best)
        rank = {"cells": [f"{w} {n}" for w, n in cells], "gate_equal_cells": bool(same), "means": means}
        rule = ("several qualified; the largest mean of X - x_e2 over the cells where every qualified arm and x_e2 "
                "have barred_form's gated count" if same else
                "several qualified; no gate-equal cell, the largest mean of X - x_e2 over the four cells")
    names = [name] if name else []
    kind = KIND_OF.get(name) if name else None
    outcome = (f"names {name}: phi_M ** x * the scaled table of {kind} on base B (build_k24_form.py --form k24{kind}-"
               "<B>), e_T and c from K24's committed pregate.json; after C1, C2, the K23 file and CK" if name
               else NO_FILE)
    return {"status": status, "qualified": q, "rank": rank, "names": names, "rule": rule,
            "files": [f"k24{kind}-B (base B of decision log 2026-11-20 (q))"] if name else [], "outcome": outcome,
            "reported": {t: s.get("report") for t, s in status.items() if s.get("report")}}


def expected_files() -> list[str]:
    return k21.expected_files()


INCOMPLETE = ("not a reading: pregate.json, guard.json or a seed's task file is missing (a crash, a refusal such as a "
              "dirty tree, out of memory or a kill); nothing was read. Row K24: if the cause is outside the code, "
              "rerun the missing seeds into this folder at the run's commit on the same machine (the runs skip what "
              "is there); after a committed code fix, rerun everything at the new commit into a fresh folder and "
              "never read this one")


def reading_code() -> dict:
    st = k18.stamp()
    h = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()  # noqa: E731
    return {"script": SCRIPT, "script_sha256": h(Path(__file__)), "outside_forms_sha256": h(HERE / "outside_forms.py"),
            "k22_sha256": h(HERE / "k22_outside_forms.py"), "k23_sha256": h(HERE / "k23_outside_exponents.py"),
            "k21_sha256": h(HERE / "k21_factor_forms.py"), "commit": st["commit"], "tree_clean": st["tree_clean"],
            "machine": st["machine"]}


def _tables_ok(erec: dict, k23_read: dict | None) -> list[str]:
    why = []
    pins = pinned_ext()
    if {"listing_sha256": erec.get("listing_sha256"), "tables": erec.get("tables")} != pins:
        why.append("pregate.json: its tables are not the pinned ones (listing or table sha256)")
    if k23_read is not None:
        ke = k23_read.get("ext") or {}
        if (ke.get("tables") or {}).get("all") != (erec.get("tables") or {}).get("two") or \
                ke.get("combine_sha256") != (erec.get("combine_sha256") or {}).get("two"):
            why.append("pregate.json: `two` is not the table K23's committed reading.json records")
    else:
        why.append("K23's committed reading.json is missing")
    return why


def reading(runs: Path, refs: Path = K18_TASKS, u3d: Path = k18.U3D_RECORDS, k23_dir: Path = K23_RESULTS) -> dict:
    """Apply the reading of row K24 to a run folder."""
    runs = Path(runs)
    if not (runs / "pregate.json").exists() or not (runs / GUARD_FILE).exists():
        miss = [f for f in ("pregate.json", GUARD_FILE) if not (runs / f).exists()]
        return {"incomplete": True, "void": False, "missing": miss, "outcome": INCOMPLETE}
    why: list[str] = []
    sha, prov = {}, {}
    pre_file, sha["pregate.json"] = k21._read_json(runs / "pregate.json")
    prov["pregate.json"] = pre_file["provenance"]
    pre = pre_file["pregate"]
    arms = active_arms(pre)
    if list(arms) != pre_file["arms"]:
        why.append("pregate.json: its arms are not the pre-gate's")
    erec = pre_file.get("ext") or {}
    try:
        _, k23_read = k23_committed(k23_dir)
    except Refused:
        k23_read = None
    why += _tables_ok(erec, k23_read)
    want_rec = {"dropped": dropped_arms(pre), "arms": list(arms), "ext": erec, "constants": constants_of(pre)}
    g, sha[GUARD_FILE] = k21._read_json(runs / GUARD_FILE)
    prov[GUARD_FILE] = g["provenance"]
    gmiss = [f for f in guard_files() if not (runs / f).exists()]
    if gmiss:
        return {"incomplete": True, "void": False, "missing": gmiss, "outcome": INCOMPLETE}
    gwhy: list[str] = []
    gsha, gprov = {}, {}
    grows = _rows_of(runs, guard_files(), f"{GUARD_DIR}/", GUARDS, gwhy, gsha, gprov, want_rec)
    if gsha != g.get("files_sha256"):
        why.append("guard.json: its files' sha256 are not the guard files'")
    why += gwhy
    sha |= gsha
    prov |= gprov
    common = {"files_sha256": sha, "provenance": prov, "arms": list(arms), "pregate": pre, "ext": erec,
              "guard_part": {k: g.get(k) for k in ("guard1", "guard2", "failed", "passes")}}
    if g.get("passes") is not True:
        failed = [x for x in (g.get("failed") or []) if x in GUARDS]
        tries = len(failed) or 1
        assert tries <= MAX_TRIES
        return {"incomplete": False, **common, "void": True, "tries_counted": tries, "worlds": {}, "names": [],
                "why": [f"the guard part failed: {failed or 'unrecorded'}", *why],
                "outcome": "void: the guard part failed before any try; no file; the K23 file stands alone"}
    missing = [f for f in expected_files() if not (runs / f).exists()]
    if missing:
        return {"incomplete": True, "void": False, "missing": missing, "outcome": INCOMPLETE}
    if pre_file.get("constants") is not None and pre_file["constants"] != want_rec["constants"]:
        why.append("pregate.json: its constants are not its pre-gate's")
    rows = _rows_of(runs, expected_files()[1:], "", arms, why, sha, prov, want_rec)
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
        per[w][n] = _summarise(rr, (harness.FLOOR, *arms) if w == "plain" else arms)
        for x in GUARDS:
            d = k18.distance(rr, x, grows.get((w, n)) or [], x)
            match[f"{w} {n} {x}"] = d
            if not k18.guard_passes(d):
                why.append(f"guard match: {w} {n}: {x} in the scoring part is off the guard part's (limits "
                           f"{GUARD_MEAN} and {GUARD_SEED})")
    check = k18.dm_check({n: rows[("dm", n)] for n in TASKS}, u3d)
    tries = try_count(arms)
    base = {"incomplete": False, **common, "guard_match": match, "dm_check_u3d_report_only": check,
            "tries_counted": tries}
    if why:
        return base | {"void": True, "why": why, "worlds": per, "names": [],
                       "outcome": "void: no file; the K23 file stands alone"}
    return (base | {"void": False, "why": [], "commit": next(iter(commits)), "machine": next(iter(machines)),
                    "worlds": per} | decide(per, arms))


def text(res: dict) -> str:
    lines = [f"Try K24: {res['outcome']}."]
    if res.get("incomplete"):
        miss = res["missing"]
        return "\n".join([*lines, f"missing: {len(miss)} file(s), e.g. " + ", ".join(miss[:6])]) + "\n"
    lines.append(f"arms: {', '.join(res['arms'])}; tries counted: {res['tries_counted']}")
    lines.append(f"tables: listing {res['ext'].get('listing_sha256')}, " + ", ".join(
        f"{n} {s}" for n, s in (res["ext"].get("tables") or {}).items()))
    lines += pregate_lines(res["pregate"])
    gp = res.get("guard_part") or {}
    lines.append(f"guard part: passes {gp.get('passes')}, failed {gp.get('failed') or 'none'}")
    for name in ("guard1", "guard2"):
        for k, v in (gp.get(name) or {}).items():
            extra = (f", mean {v['mean']:.2f} against {v['k23_mean']:.2f}, gated {v['gated']} against "
                     f"{v['k23_gated']}" if name == "guard2" and v.get("mean") is not None else "")
            lines.append(f"  {name} {k}: passes {v['passes']}; {k21._dist(v.get('distance'))}{extra}")
    if res["void"]:
        lines.append("VOID: " + "; ".join(res["why"]))
    else:
        lines.append(f"commit {res['commit']}, machine {res['machine']}.")
        for t, s in res["status"].items():
            if not s["ran"]:
                lines.append(f"  {t}: not run (dropped by the pre-gate)")
                continue
            lines.append(f"  {t}: qualified {'yes' if s['qualified'] else 'no'} ((a) {s['a']}, (b) {s['b']}, (c) "
                         f"{s['c']}, (d) {s['d']}, (e) {s['e']}, (f) {s['f']}); ahead in {s['ahead_cells'] or 'none'}"
                         f"; behind in {s['behind_cells'] or 'none'}; X - x_e2 {s['head_to_head']}; gated (arm, "
                         f"base) {s['gated_vs_base']}" + (f"; {s['report']}" if s.get("report") else ""))
        lines.append(f"  qualified: {res['qualified'] or 'none'}; rule: {res['rule']}; rank {res['rank']}; named: "
                     f"{res['files'] or 'none'}")
    for k, d in (res.get("guard_match") or {}).items():
        lines.append(f"  scoring part against the guard part, {k}: {k21._dist(d)}")
    lines.append("U3d's check arm in the dm world, reported only: " + "; ".join(
        f"{n} {k21._dist(c['distance'])}" for n, c in (res.get("dm_check_u3d_report_only") or {}).items()))
    lines.append("task files: " + ", ".join(f"{f} {h[:12]}" for f, h in res["files_sha256"].items()))
    lines.append("g1 and g2 are two replicates of one knockout run; the controls, the twins, x_e2 and barred_form are "
                 "never named.")
    for w in WORLDS:
        lines.append(f"\nWORLD {w} ({k18.WORLD_NAMES[w]}):")
        for n, t in TASKS.items():
            v = (res.get("worlds") or {}).get(w, {}).get(n)
            if v is None:
                lines.append(f"{w} {n} {t.key}: not read")
                continue
            lines.append(f"{w} {n} {t.key}:")
            for a, s in v["arms"].items():
                lines.append(f"    {a:<13} {s['mean']:.2f} (sd {s['sd']:.2f}); severity scored {s['severity_scored']}"
                             f", gated {s['severity_gated']} {s['severity_states']}, median beta "
                             f"{k21._f(s['median_beta'])}, r2 {k21._f(s['median_r2'])}; skills "
                             + ", ".join(f"{m} {x:.3f}" for m, x in s["skills"].items()))
            for k, d in v["paired"].items():
                lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})")
            for a, terms in v["per_term"].items():
                lines.append(f"    per-term skill, {a} - {BASE}: " + ", ".join(f"{m} {x:+.3f}" for m, x in
                                                                              terms.items()))
    if res.get("reading_code"):
        rc = res["reading_code"]
        lines.append(f"\nread by {rc['script']} {rc['script_sha256'][:12]} (outside_forms "
                     f"{rc['outside_forms_sha256'][:12]}, k22 {rc['k22_sha256'][:12]}, k23 {rc['k23_sha256'][:12]}, "
                     f"k21 {rc['k21_sha256'][:12]}) at commit {rc['commit']}, clean {rc['tree_clean']}, machine "
                     f"{rc['machine']}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    out = Path(a.out)
    targets = [out / "reading.json", out / "reading.txt"]
    if any(p.exists() for p in targets):
        print(f"refused: {out} holds a reading already; it is never overwritten (move it, or read into a fresh "
              "folder)", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    try:
        res = reading(Path(a.runs))
        if not res.get("incomplete"):
            res["reading_code"] = reading_code()
        t = text(res)
        js = json.dumps(res, indent=1, default=float)
    except Exception:  # noqa: BLE001  (masked traceback, nothing written)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    if res.get("incomplete"):
        print(t, end="")
        sys.exit(EXIT_INCOMPLETE)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(js)
    (out / "reading.txt").write_text(t)
    print(t, end="")


# ---------------------------------------------------------------------------------------------------------------------
# The listing and the pins (for the main session)


def listing_text(env: dict | None = None) -> str:
    """TABLES_K24.sha256's text: each key's file hashed again from disk and checked against its COMBINE record."""
    lines = []
    for key in KEYS:
        var, name, cfile = SOURCES[key]
        d = folder(var, env)
        comb, _ = combine_record(d, cfile)
        if name not in comb:
            raise Refused(f"{key}: {cfile} names no table {name}: refused")
        for q in (d / name, d / name / "factors"):
            if q.is_symlink() or not q.is_dir():
                raise Refused(f"{key}: a folder on its path is missing or a symlink: refused")
        raw = _plain_file(d / name / "factors" / f"{name}.csv", f"{key}'s pooled table").read_bytes()
        sha = sha256_bytes(raw)
        if sha != comb[name]:
            raise Refused(f"{key}: the file hashes to {sha[:12]}, {cfile} says {comb[name][:12]}: refused")
        lines.append(f"{sha}  {key}\n")
    return "".join(lines)


def pins_text(path: Path | None = None) -> str:
    path = LISTING_FILE if path is None else Path(path)
    raw = path.read_bytes()
    listing = ext.parse_listing(raw)
    if sorted(listing) != list(KEYS):
        raise Refused(f"{path.name} lists {sorted(listing)}, not exactly {list(KEYS)}: refused")
    items = ", ".join(f'"{k}": "{listing[k]}"' for k in ("two", "sc3", "bulk"))
    return (f'K24_LISTING_SHA256: str | None = "{sha256_bytes(raw)}"\n'
            f"K24_TABLE_SHA256: dict[str, str | None] = {{{items}}}\n")


def main(argv=None) -> None:
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv in (["listing"], ["pins"]):
        try:
            sys.stdout.write(listing_text() if argv == ["listing"] else pins_text())
        except (Refused, ext.TableRefused) as e:
            refuse(str(e))
        return
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
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
    for q in (p, r):
        q.add_argument("--data", default=str(harness.DATA_DIR))
        q.add_argument("--run-dir", default=str(k18.RUN_DIR))
        q.add_argument("--decision", default=str(k18.DECISION))
    gd = sub.add_parser("guard")
    gd.add_argument("--runs", required=True)
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"pregate": cmd_pregate, "run": cmd_run, "guard": cmd_guard, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
