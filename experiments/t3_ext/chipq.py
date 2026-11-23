"""CHIPQ-S, the database part: the kept list, the binding table, the placebo pool and the uncounted pre-gate P1.

The specification (the cold challenger's final text for T3 candidate CHIPQ-S, register row X018), steps 2 to 7:
  2. The two readers' kept experiment ids per query antigen are committed as chipq_kept.yaml with its sha256 in
     CHIPQ_KEPT.sha256 (one line, sha256sum form). Every step refuses (exit 2) while either file is absent, untracked,
     changed against HEAD, or the two disagree; committing them changes no line of code.
  3. E_q = the kept experiments' columns of q's +/-5 kb target-gene file, genes mapped to the panel by exact symbol in
     the barred factor's gene order (transfer_mab21l2_all.csv, loaded as row K24's code loads it); a panel gene absent
     from the file scores 0 in every experiment; s_g = the mean over E_q of the per-experiment scores.
  4. B_q = panel genes with a score above 0 in at least half of E_q (at least one), ranked by s_g;
     top K = min(40, |B_q|).
  5. w_g falls linearly from 1 (top gene) to 0.5 (K-th); every other gene 0.
  6. C_q = exp(-w) as a raw table, finished through t3x_tables.finish (exactly 1 at w = 0, clipped to [0.25, 4]).
  7. e = S / SD(ln C_q), population SD over the panel genes in file order, at S* and 4 S*.
The placebo pool: antigens with at least 2 screened-clean experiments (the keyword filter of chipq_fetch.py) and |B|
between 0.75 K and 1.25 K, excluding q, Ctnnb1, Gata4 and Gata6; for seed s one antigen drawn with seed 20261124 + s,
without replacement through the seeds in order (the rule writer's settled point A1), from a pool of at least 10, built
by steps 3 to 7.

Every constant the specification leaves to a reading is one named constant below, in this one place (the rule writer's
settled points, points_settled.md section B, say for each whether the first reading stands or changes; see the block
"READINGS"). No constant comes from a score.

ONE function opens a target-gene file: open_target. It checks the file's bytes against the fetch manifest (sha256),
takes the experiments from the pinned kept list (a query antigen: chipq_kept.yaml; a placebo antigen: the keyword
filter's placebo_kept.json, itself checked against the manifest) and refuses any other experiment, parses the header,
and then reads only the kept experiments' columns: never the Average or the STRING column. No message carries a gene
symbol, a score or the column of a dropped experiment.

Steps (python -I -B experiments/t3_ext/chipq.py <step>):
  p1      prints exactly two integers for Ctnnb1: the number of kept experiments and |B|. Exit 0 when both bars hold
          (P1_MIN_EXPERIMENTS, P1_MIN_BOUND), 3 when not, 2 refused. It writes nothing.
  tables  Gata4's raw table and the per-seed placebo tables (seeds vec.adoption.SEEDS), named chipq_Gata4 and
          chipq_rand_seed0 to chipq_rand_seed9, into the folder the environment variable VEC_T3_EXT_TABLES_K27 names
          (outside every checkout; no path is taken from the command line), each with its manifest
          (t3x_tables.write_table, the outside tables' writer; a placebo's manifest records its target file's sha256
          and pool index, never its name); prints file names and sha256 only. Never Ctnnb1's table: that belongs to
          the build path after the captain's ruling. The committed listing experiments/t3_ext/TABLES_K27.sha256 is
          written by the read's engineer, not here.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import chipq_fetch as cf  # noqa: E402
import t3x_tables as tb  # noqa: E402

Refused = cf.Refused

KEPT_YAML = HERE / "chipq_kept.yaml"
KEPT_SHA = HERE / "CHIPQ_KEPT.sha256"
TEST_QUERY = "Ctnnb1"
PRACTICE_QUERY = "Gata4"
QUERY_ANTIGENS = cf.QUERY_ANTIGENS
BARRED_NAME = "transfer_mab21l2_all"
PHI_M_SHA256 = "64b1e0b8572c7193ed1e1f0b7339a26081e601f1e706e6b08193e1718dc5c7e8"   # build_outside_form.TABLE_SHA256

# ---------------------------------------------------------------------------------------------------------------------
# READINGS. The specification's constants, and the literal reading of each point it leaves open, in one place.
K_MAX = 40                     # step 4: K = min(40, |B_q|)
BOUND_SHARE = 0.5              # step 4: score > 0 in at least half of E_q (2 * count >= n), and in at least one
W_TOP, W_BOTTOM = 1.0, 0.5     # step 5: linear from 1 (rank 1) to 0.5 (rank K); K = 1 gives w = 1
RANK_TIES = "file order"       # step 4: equal s_g keep the barred factor's gene order (a stable sort); open point
FINISH_DEN = None              # step 6: t3x_tables.finish's den. None = a zero profile, so finish does exp, clip and
#                                absent = 1 and no rescale (the specification lists only those two); open point
S_STAR = 0.04949423311987974   # step 7: K23's committed size (k24_wider_tables.S_STAR; a test pins the equality)
SIZES = {"S*": S_STAR, "4S*": 4 * S_STAR}
POOL_MIN_KEPT = cf.PLACEBO_MIN_KEPT   # "at least 2 screened-clean experiments"
POOL_EXCLUDE = ("Ctnnb1", "Gata4", "Gata6")   # and q itself
POOL_BAND = (0.75, 1.25)       # |B_p| in [0.75 K_q, 1.25 K_q], both ends included, K_q = the query's K; open point
PLACEBO_SEED_BASE = 20261124   # seed 20261124 + s
POOL_MIN_SIZE = 10             # the pool holds at least 10 antigens, or the step refuses (settled A1, D's P2 bar)
PLACEBO_DRAW = ("without replacement over the pool sorted by name: for s in vec.adoption.SEEDS in order, "
                "remaining[numpy.random.default_rng(20261124 + s).integers(len(remaining))], remaining = the sorted "
                "pool minus the antigens drawn for earlier seeds")   # settled A1 (reading B2 changed)
P1_MIN_EXPERIMENTS = 3         # P1: at least 3 tier-U Ctnnb1 experiments
P1_MIN_BOUND = 15              # P1: |B| at least 15
TABLE_QUERY = "chipq_{q}"      # the tables step's file names, kept (settled A5)
TABLE_PLACEBO = "chipq_rand_seed{s}"
ENV_TABLES = "VEC_T3_EXT_TABLES_K27"   # the tables step's folder, ~/v/vec-ext/t3/tables_k27 on the DGX (settled A5,
#                                        reading B6 changed): named only by this variable, outside every checkout

_SHA_LINE = re.compile(r"([0-9a-f]{64})  chipq_kept\.yaml\n?")


# ---------------------------------------------------------------------------------------------------------------------
# The kept list


def _git(path: Path, *args: str) -> int:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["GIT_OPTIONAL_LOCKS"] = "0"
    return subprocess.run(["git", "-C", str(path.parent), *args, "--", path.name], capture_output=True, env=env,
                          check=False).returncode


def committed(path: Path) -> bool:
    """Tracked by git and unchanged against HEAD."""
    return _git(path, "ls-files", "--error-unmatch") == 0 and _git(path, "diff", "--quiet", "HEAD") == 0


class Kept:
    """The pinned kept list as load_kept returns it (open_target takes nothing else): ids per query antigen, and the
    sha256 the list was checked against."""

    __slots__ = ("_ids", "sha256")

    def __init__(self, ids: dict[str, list[str]], sha256: str):
        self._ids = {a: tuple(v) for a, v in ids.items()}
        self.sha256 = sha256

    def __contains__(self, a) -> bool:
        return a in self._ids

    def __getitem__(self, a) -> tuple[str, ...]:
        return self._ids[a]


def load_kept(yaml_path: Path | None = None, sha_path: Path | None = None, check_committed: bool = True) -> Kept:
    """{query antigen: kept experiment ids}, from the committed chipq_kept.yaml at the sha256 of CHIPQ_KEPT.sha256."""
    import yaml

    yp = KEPT_YAML if yaml_path is None else Path(yaml_path)
    sp = KEPT_SHA if sha_path is None else Path(sha_path)
    for p in (yp, sp):
        if p.is_symlink() or not p.is_file():
            raise Refused(f"{p.name} is absent (the main session commits it after the screen): refused")
        if check_committed and not committed(p):
            raise Refused(f"{p.name} is not committed (untracked, or changed against HEAD): refused")
    m = _SHA_LINE.fullmatch(sp.read_text(encoding="ascii", errors="replace"))
    if m is None:
        raise Refused(f"{sp.name} is not one line '<sha256>  chipq_kept.yaml': refused")
    raw = yp.read_bytes()
    if cf.sha256_bytes(raw) != m.group(1):
        raise Refused(f"{yp.name} does not hash to the sha256 in {sp.name}: refused")
    d = yaml.safe_load(raw)
    if not isinstance(d, dict) or not d or not set(d) <= set(QUERY_ANTIGENS):
        raise Refused(f"{yp.name}: not a mapping of query antigens {QUERY_ANTIGENS} to experiment ids: refused")
    out = {}
    for a, ids in d.items():
        if not isinstance(ids, list) or not all(isinstance(i, str) and cf.EXP_ID.fullmatch(i) for i in ids):
            raise Refused(f"{yp.name}: {a}'s entry is not a list of experiment ids: refused")
        if ids != sorted(set(ids)):
            raise Refused(f"{yp.name}: {a}'s ids are not sorted and unique: refused")
        out[str(a)] = list(ids)
    return Kept(out, m.group(1))


# ---------------------------------------------------------------------------------------------------------------------
# The one reader of a target-gene file

_EXP_COL = re.compile(r"([SDE]RX\d+)\|.*")


def open_target(antigen: str, genes: list[str], root=None, kept: Kept | None = None, ids: list[str] | None = None,
                man: dict | None = None) -> tuple[np.ndarray, list[str]]:
    """(E, ids): E[g, j] the score of panel gene g (in `genes`' order) in kept experiment ids[j]; 0 where the file has
    no row for the gene. `ids` (default: every pinned kept id) must lie in the pinned kept list.

    A query antigen's pinned ids (the readers' list) must each have a column, or the call is refused. For a placebo
    antigen with the default `ids`, the experiments used are the filter-kept ids that have a column in its target file
    (the specification's step 3: "E_q = the kept SRX columns of q's +/-5 kb file"), in the pinned order; that set may
    be smaller than the kept list, or empty (E then has no column). Ids asked for by name must each have a column."""
    root = cf.root_of(root)
    man = cf.load_manifest(root) if man is None else man
    if antigen in QUERY_ANTIGENS:
        kept = load_kept() if kept is None else kept
        if not isinstance(kept, Kept):
            raise Refused("the kept list must come from load_kept (the pinned chipq_kept.yaml): refused")
        if antigen not in kept:
            raise Refused("the query antigen has no entry in the pinned kept list: refused")
        pinned = list(kept[antigen])
    else:
        if kept is not None:
            raise Refused("a placebo antigen takes its kept ids from placebo_kept.json only: refused")
        pk = cf.placebo_kept(root, man)
        if antigen not in pk:
            raise Refused("the antigen has no entry in placebo_kept.json: refused")
        pinned = list(pk[antigen])
    columned_only = antigen not in QUERY_ANTIGENS and ids is None
    ids = list(pinned) if ids is None else list(ids)
    if any(i not in pinned for i in ids):
        raise Refused("an experiment not in the pinned kept list was asked for: refused")
    if not ids or len(set(ids)) != len(ids):
        raise Refused("no kept experiment, or one asked for twice: refused")
    raw = cf.checked_bytes(root, cf.target_rel(antigen), man)
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        raise Refused("the target file is not UTF-8: refused") from None
    lines = text.split("\n")
    head = lines[0].rstrip("\r").split("\t")
    if (len(head) < 4 or head[0] != "Target_genes" or head[1] != f"{antigen}|Average" or head[-1] != "STRING"):
        raise Refused("the target file's header is not of the fixed form: refused")
    exp_ids = []
    for h in head[2:-1]:
        m = _EXP_COL.fullmatch(h)
        if m is None:
            raise Refused("a target-file experiment column is not of the form <id>|<cell-type class>: refused")
        exp_ids.append(m.group(1))
    if len(set(exp_ids)) != len(exp_ids):
        raise Refused("the target file holds an experiment column twice: refused")
    col = {e: 2 + j for j, e in enumerate(exp_ids)}
    if columned_only:
        ids = [i for i in ids if i in col]
    missing = sum(1 for i in ids if i not in col)
    if missing:
        raise Refused(f"{missing} kept experiment(s) have no column in the target file: refused")
    cols = [col[i] for i in ids]
    pos = {g: k for k, g in enumerate(genes)}
    if len(pos) != len(genes):
        raise Refused("the panel gene list holds a gene twice: refused")
    E = np.zeros((len(genes), len(ids)), dtype=np.float64)
    seen = np.zeros(len(genes), dtype=bool)
    for line in lines[1:]:
        line = line.rstrip("\r")
        if not line:
            continue
        f = line.split("\t")
        if len(f) != len(head):
            raise Refused("a target-file row has another number of fields than the header: refused")
        k = pos.get(f[0])
        if k is None:
            continue
        if seen[k]:
            raise Refused("a panel gene has two rows in the target file: refused")
        seen[k] = True
        try:
            v = [float(f[c]) for c in cols]
        except ValueError:
            raise Refused("a kept experiment's value is not a number: refused") from None
        if not all(np.isfinite(x) and x >= 0 for x in v):
            raise Refused("a kept experiment's value is not finite and non-negative: refused")
        E[k] = v
    return E, ids


# ---------------------------------------------------------------------------------------------------------------------
# Steps 3 to 7


def bound_rank(E: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(s_g, the bound mask, the bound genes' indices ranked by s_g descending, ties in file order)."""
    E = np.asarray(E, dtype=np.float64)
    n = E.shape[1]
    if n < 1:
        raise Refused("no kept experiment: refused")
    s = E.mean(axis=1)
    hits = (E > 0).sum(axis=1)
    bound = (hits >= BOUND_SHARE * n) & (hits >= 1)
    idx = np.flatnonzero(bound)
    order = idx[np.lexsort((idx, -s[idx]))]
    return s, bound, order


def weights(order: np.ndarray, n_genes: int) -> tuple[np.ndarray, int]:
    """(w, K): w falls linearly from W_TOP to W_BOTTOM over the top K = min(K_MAX, |B|), 0 elsewhere."""
    K = int(min(K_MAX, len(order)))
    w = np.zeros(n_genes, dtype=np.float64)
    if K == 1:
        w[order[0]] = W_TOP
    elif K > 1:
        w[order[:K]] = W_TOP - (W_TOP - W_BOTTOM) * np.arange(K) / (K - 1)
    return w, K


def raw_table(w: np.ndarray) -> np.ndarray:
    """C = exp(-w), finished through t3x_tables.finish: genes with w = 0 are absent (exactly 1), clip [0.25, 4]."""
    w = np.asarray(w, dtype=np.float64)
    L = np.where(w > 0, -w, np.nan)
    den = np.zeros(len(w)) if FINISH_DEN is None else np.asarray(FINISH_DEN, dtype=np.float64)
    return tb.finish(L, den)


def exponent(C, S: float) -> float | None:
    """e = S / SD(ln C), population SD over all genes in file order; None when that SD is 0 (as k24's exponent)."""
    lt = np.log(np.asarray(C, dtype=np.float64))
    sd = float(np.std(lt))
    return S / sd if sd > 0 and float(np.ptp(lt)) > 0 else None


def build(antigen: str, genes: list[str], root=None, kept: Kept | None = None, man: dict | None = None) -> dict:
    """Steps 3 to 7 for one antigen."""
    E, ids = open_target(antigen, genes, root, kept, None, man)
    s, bound, order = bound_rank(E)
    w, K = weights(order, len(genes))
    C = raw_table(w)
    return {"n_kept": len(ids), "n_bound": int(bound.sum()), "K": K, "w": w, "C": C,
            "exponents": {k: exponent(C, S) for k, S in SIZES.items()}}


# ---------------------------------------------------------------------------------------------------------------------
# The placebo pool and the per-seed placebo


def placebo_pool(q: str, K_q: int, genes: list[str], root=None, man: dict | None = None) -> list[str]:
    """Antigens of placebo_kept.json with at least POOL_MIN_KEPT kept experiments that have a column in their target
    file (the experiments open_target uses for a placebo) and |B| in POOL_BAND * K_q, excluding q and POOL_EXCLUDE;
    sorted by name."""
    root = cf.root_of(root)
    man = cf.load_manifest(root) if man is None else man
    lo, hi = POOL_BAND[0] * K_q, POOL_BAND[1] * K_q
    pool = []
    for a, ids in sorted(cf.placebo_kept(root, man).items()):
        if a == q or a in POOL_EXCLUDE or a in QUERY_ANTIGENS or len(ids) < POOL_MIN_KEPT:
            continue   # fewer than two kept ids: its target file was never fetched
        E, used = open_target(a, genes, root, None, None, man)
        if len(used) < POOL_MIN_KEPT:
            continue
        n_b = int(bound_rank(E)[1].sum())
        if lo <= n_b <= hi:
            pool.append(a)
    return pool


def placebo_draws(pool: list[str]) -> list[tuple[int, str, int]]:
    """[(seed, antigen, its index in the sorted pool)] for every seed of vec.adoption.SEEDS, in that order (settled
    point A1): for seed s the antigen is remaining[default_rng(PLACEBO_SEED_BASE + s).integers(len(remaining))], where
    remaining is the pool sorted by name minus the antigens drawn for earlier seeds. Ten seeds give ten different
    placebo antigens. Refused unless the pool is sorted, unique and holds at least POOL_MIN_SIZE antigens."""
    from vec.adoption import SEEDS

    pool = list(pool)
    if pool != sorted(set(pool)):
        raise Refused("the placebo pool is not sorted by name and unique: refused")
    if len(pool) < max(POOL_MIN_SIZE, len(SEEDS)):
        raise Refused(f"the placebo pool holds {len(pool)} antigens, fewer than {POOL_MIN_SIZE}: refused")
    remaining = list(pool)
    out = []
    for s in SEEDS:
        a = remaining.pop(int(np.random.default_rng(PLACEBO_SEED_BASE + int(s)).integers(len(remaining))))
        out.append((int(s), a, pool.index(a)))
    return out


def placebo_for_seed(pool: list[str], s: int) -> tuple[str, int]:
    """(antigen, its index in the sorted pool) for seed s, from placebo_draws (the draws of the earlier seeds come
    first, so a seed's antigen does not depend on which other seeds are asked for)."""
    for seed, a, i in placebo_draws(pool):
        if seed == int(s):
            return a, i
    raise Refused("the seed is not one of vec.adoption.SEEDS: refused")


# ---------------------------------------------------------------------------------------------------------------------
# The panel genes in the barred factor's file order, through the loader row K24's code uses


def barred_genes(data: Path, run_dir: Path, decision: Path) -> list[str]:
    """WT E8.75's genes as k18.load_barred_for_plain loads them (load_stage, duplicates dropped), with the barred
    factor (transfer_mab21l2_all) loaded through build_submission.load_factor in that order and at PHI_M_SHA256."""
    t3 = REPO / "experiments" / "t3_decomp"
    if str(t3) not in sys.path:
        sys.path.insert(0, str(t3))
    import k18_built_files as k18

    rec = k18.load_barred_for_plain(Path(data), Path(run_dir), Path(decision))
    if rec.get("name") != BARRED_NAME or rec.get("sha256") != PHI_M_SHA256:
        raise Refused("the barred factor is not transfer_mab21l2_all at its pinned sha256: refused")
    return list(k18._BARRED["genes"])


# ---------------------------------------------------------------------------------------------------------------------
# Steps


def p1(genes: list[str], root=None, kept: Kept | None = None) -> tuple[int, int, bool]:
    """(kept experiments, |B|, both bars hold) for the test query. With no kept experiment the file is not opened."""
    kept = load_kept() if kept is None else kept
    if TEST_QUERY not in kept:
        raise Refused("the pinned kept list has no entry for the test query: refused")
    if not kept[TEST_QUERY]:
        return 0, 0, False
    r = build(TEST_QUERY, genes, root, kept)
    return r["n_kept"], r["n_bound"], r["n_kept"] >= P1_MIN_EXPERIMENTS and r["n_bound"] >= P1_MIN_BOUND


def tables_folder(env: dict | None = None) -> Path:
    """The tables step's folder: the environment variable ENV_TABLES, set, absolute, no symlink, an existing directory,
    outside every checkout (t3x_common.outside_repo). Returns its realpath. No other path is taken."""
    env = os.environ if env is None else env
    v = env.get(ENV_TABLES)
    if not v:
        raise Refused(f"the environment variable {ENV_TABLES} is not set: refused (no table path is taken otherwise)")
    p = Path(v)
    if not p.is_absolute():
        raise Refused(f"{ENV_TABLES} is not an absolute path: refused")
    if p.is_symlink():
        raise Refused(f"{ENV_TABLES} names a symlink: refused")
    if not p.is_dir():
        raise Refused(f"{ENV_TABLES} does not name a directory: refused")
    try:
        return tb.cm.outside_repo(Path(os.path.realpath(p)))
    except tb.cm.Refused as e:
        raise Refused(f"{ENV_TABLES}: {e}") from None


def tables(genes: list[str], out=None, seeds=None, root=None, kept: Kept | None = None) -> list[tuple[str, str]]:
    """Gata4's raw table and the per-seed placebo tables, written with their manifests; [(file name, sha256)]. `out`
    None (the step's only form) is tables_folder(); a folder given here (tests) must lie outside every checkout."""
    from vec.adoption import SEEDS

    seeds = SEEDS if seeds is None else seeds
    root = cf.root_of(root)
    man = cf.load_manifest(root)
    kept = load_kept() if kept is None else kept
    out = tables_folder() if out is None else tb.cm.outside_repo(out)
    q = PRACTICE_QUERY
    rq = build(q, genes, root, kept, man)
    pool = placebo_pool(q, rq["K"], genes, root, man)
    jobs = [(TABLE_QUERY.format(q=q), rq, {"antigen": q})]
    draws = {seed: (a, i) for seed, a, i in placebo_draws(pool)}
    for s in seeds:
        if int(s) not in draws:
            raise Refused("the seed is not one of vec.adoption.SEEDS: refused")
        a, i = draws[int(s)]
        rp = build(a, genes, root, None, man)
        jobs.append((TABLE_PLACEBO.format(s=s), rp,
                     {"seed": int(s), "pool_size": len(pool), "pool_index": i,
                      "target_sha256": man["files"][cf.target_rel(a)]["sha256"]}))
    for name, _, _ in jobs:
        if (out / "factors" / f"{name}.csv").exists() or (out / f"{name}_manifest.json").exists():
            raise Refused(f"{name}: already written in the tables folder: nothing is overwritten")
    written = []
    for name, r, extra in jobs:
        manifest = {"candidate": "CHIPQ-S", "register_row": "X018", "query": q, **extra,
                    "n_kept": r["n_kept"], "n_bound": r["n_bound"], "K": r["K"], "exponents": r["exponents"],
                    "readings": readings(), "fetch_manifest_sha256": cf.sha256_file(root / "MANIFEST.json")}
        files = tb.write_table(out, name, genes, r["C"], None, manifest)
        written.append((f"{name}.csv", files[name]["sha256"]))
    return written


def readings() -> dict:
    return {"K_MAX": K_MAX, "BOUND_SHARE": BOUND_SHARE, "W": [W_TOP, W_BOTTOM], "RANK_TIES": RANK_TIES,
            "FINISH_DEN": "zero profile (no rescale)" if FINISH_DEN is None else "given", "SIZES": SIZES,
            "POOL_MIN_KEPT": POOL_MIN_KEPT, "POOL_EXCLUDE": list(POOL_EXCLUDE), "POOL_BAND": list(POOL_BAND),
            "POOL_MIN_SIZE": POOL_MIN_SIZE, "PLACEBO_DRAW": PLACEBO_DRAW}


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="step", required=True)
    for name in ("p1", "tables"):
        q = sub.add_parser(name)
        q.add_argument("--data", default=None)
        q.add_argument("--run-dir", default=None)
        q.add_argument("--decision", default=None)
    a = ap.parse_args(argv)
    import contextlib

    try:
        # anything a loader prints goes to stderr: stdout carries only this step's own lines
        with contextlib.redirect_stdout(sys.stderr):
            kept = load_kept()
            t3 = REPO / "experiments" / "t3_decomp"
            if str(t3) not in sys.path:
                sys.path.insert(0, str(t3))
            import k18_built_files as k18

            from vec import harness

            genes = barred_genes(Path(a.data or harness.DATA_DIR), Path(a.run_dir or k18.RUN_DIR),
                                 Path(a.decision or k18.DECISION))
            if a.step == "p1":
                n, b, ok = p1(genes, None, kept)
            else:
                written = tables(genes, tables_folder(), None, None, kept)
    except (Refused, tb.Refused) as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2
    except Exception as e:  # noqa: BLE001  (another module's message could name a gene: only its type is shown)
        print(f"REFUSED: {type(e).__name__} in a loader", file=sys.stderr)
        return 2
    if a.step == "p1":
        print(f"{n} {b}")
        return 0 if ok else 3
    for fname, sha in written:
        print(f"{sha}  {fname}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
