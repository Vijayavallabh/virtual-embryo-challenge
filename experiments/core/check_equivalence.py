#!/usr/bin/env python
"""feat-020 equivalence check on released data: T1 and T2 builder outputs, hashed, from whichever vec package is first
on the path. Run it once with the code from before the move into vec.core and once with the code after it, then
compare the two JSONs (experiments/core/REPORT.md has the commands and the results).

  PYTHONPATH=<copy>/src <python> experiments/core/check_equivalence.py run --data <repo>/data/raw --out OUT.json
  <python> experiments/core/check_equivalence.py compare OLD.json NEW.json
  <python> experiments/core/check_equivalence.py source --base eb946ac      (in a git checkout, no data)

`run` prints the vec package it imported and records the sha256 of its builder files. It reads the released files
read-only and writes only OUT.json: the T1 builder writes no file, and the T2 builder checks its file in a temporary
directory without a build stamp. Board specs come through vec.boards (gene order and cell limits; no anchor is read).
The T1 shares are an arbitrary input for this check (the reference's own category shares, re-weighted so that some
quotas exceed their pools), not a prediction of any stage. Every hash is bit for bit: X, obs, var, obsm, layers, the
sampled rows and coordinates, and the reports with their floats as written by json (repr, exact).

`source` compares the source of every function in src/vec/core.py with where it came from: the moved functions
with their full source in vec.t1_scope / vec.t2_build at the base commit (identical except the name), and the T3
primitives with the frozen copies in experiments/t3_decomp/t3_decompose.py (identical def line except the name, and
identical body; only the docstrings differ).
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import platform
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

T1_BOARD = "T1:val"
T1_REF, T1_START = "E9.5_RNA", "E8.5_RNA"
T1_CELLS = 5000
T1_DAMPING = 0.5
T2_POOLS = {"T2:heart:val_extrap": "E9.5", "T2:heart:val_interp": "E8.75"}  # board -> released heart pool
SEEDS = (0, 1)
BUILDER_FILES = ("t1_scope.py", "t2_build.py", "core.py", "prep.py", "io.py", "boards.py", "spatial.py")
REPO = Path(__file__).resolve().parents[2]
MOVED = (("src/vec/t1_scope.py", "quotas", "quotas"), ("src/vec/t1_scope.py", "_draw", "draw_least_used"),
         ("src/vec/t1_scope.py", "counts_of", "counts_of"), ("src/vec/t1_scope.py", "mean_cp10k", "mean_cp10k"),
         ("src/vec/t1_scope.py", "apply_factors", "apply_factors"),
         ("src/vec/t2_build.py", "_balanced_rows", "draw_balanced"))
FROZEN = (("alloc", "alloc"), ("draw_rows", "draw_then_repeat"), ("thin", "thin"), ("rescale_depth", "rescale_depth"),
          ("apply_gene_factor", "apply_gene_factor"))  # name in experiments/t3_decomp/t3_decompose.py -> in the core


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def plain(o):
    """JSON-able values, numpy scalars and arrays as Python values; floats stay exact (json writes their repr)."""
    import numpy as np

    if isinstance(o, dict):
        return {str(k): plain(v) for k, v in o.items()}
    if isinstance(o, list | tuple):
        return [plain(v) for v in o]
    if isinstance(o, np.ndarray):
        return plain(o.tolist())
    if isinstance(o, np.generic):
        return o.item()
    return o


def h_json(o) -> str:
    return sha(json.dumps(plain(o), sort_keys=True).encode())


def h_array(a) -> str:
    import numpy as np
    from scipy import sparse

    if sparse.issparse(a):
        assert a.format in ("csr", "csc"), a.format
        return h_json([a.format, a.dtype.str, a.shape, h_array(a.data), h_array(a.indices), h_array(a.indptr)])
    a = np.asarray(a)
    assert a.dtype != object
    return h_json([a.dtype.str, a.shape, sha(np.ascontiguousarray(a).tobytes())])


def h_frame(df) -> str:
    cols = [[str(c), repr(df[c].dtype), [repr(v) for v in df[c].tolist()]] for c in df.columns]
    return h_json([[repr(i) for i in df.index], repr(df.index.dtype), cols])


def h_adata(a) -> dict[str, str]:
    layers = {k: h_array(a.layers[k]) for k in sorted(k for k in a.layers.keys() if k is not None)}
    return {"X": h_array(a.X), "obs": h_frame(a.obs), "var": h_frame(a.var),
            "obsm": h_json({k: h_array(a.obsm[k]) for k in sorted(a.obsm.keys())}), "layers": h_json(layers),
            "n_cells": int(a.n_obs)}


def t1_shares(ref, spec) -> dict[str, float]:
    """The reference's category shares with the chamber categories weighted x3 and epicardium x8, renormalised."""
    from vec import t1_scope as ts

    comp = ts.composition(ref)
    chamber = set(spec.in_scope("chamber"))
    w = {c: s * (8.0 if c == "epicardium" else 3.0 if c in chamber else 1.0) for c, s in comp.items()}
    total = sum(w.values())
    return {c: v / total for c, v in w.items()}


def run_t1(data: Path) -> dict:
    from vec import t1_scope as ts
    from vec.boards import get_board
    from vec.prep import load_stage

    board = get_board(T1_BOARD)
    ref0 = load_stage(data / "T1" / f"{T1_REF}.h5ad")
    start0 = load_stage(data / "T1" / f"{T1_START}.h5ad")
    probe = ref0.copy()
    spec = ts.load_spec()
    ts.classify(probe, source=T1_REF, spec=spec)
    shares = t1_shares(probe, spec)
    out = {"board": T1_BOARD, "ref": T1_REF, "start": T1_START, "ref_cells": int(ref0.n_obs),
           "start_cells": int(start0.n_obs), "n_cells": T1_CELLS, "damping": T1_DAMPING, "shares": shares}
    for seed in SEEDS:
        for name, kw in (("shift", {"start": start0.copy(), "damping": T1_DAMPING, "start_source": T1_START}),
                         ("no_shift", {})):
            t0 = time.perf_counter()
            sub, rep = ts.build_submission(ref0.copy(), board, shares, T1_CELLS, seed, ref_source=T1_REF, **kw)
            rep = {**rep, "config": Path(rep["config"]).name}  # the copy's absolute path otherwise
            reused = [r["category"] for r in rep["reweight"]["rows"] if r["reused"] > 0]
            shifted = [r["type"] for r in (rep["shift"] or {}).get("rows", []) if r["status"] == "shifted"]
            out[f"{name}_seed{seed}"] = {**h_adata(sub), "report": h_json(rep), "reused_categories": reused,
                                         "shifted_types": shifted, "seconds": round(time.perf_counter() - t0, 1)}
            print(f"T1 {name} seed {seed}: {sub.n_obs} cells x {sub.n_vars} genes, "
                  f"X {out[f'{name}_seed{seed}']['X'][:16]}, reused {reused}, shifted {shifted}", flush=True)
    return out


def run_t2(data: Path) -> dict:
    from vec import t2_build as tb
    from vec.boards import get_board
    from vec.prep import load_stage

    out = {}
    for key, stage in T2_POOLS.items():
        board = get_board(key)
        pool = load_stage(data / "T2" / f"{stage}.h5ad", genes=board.genes, drop_duplicates=True)
        n, count_record = tb.n_cells(board, pool.n_obs)
        res = {"pool": stage, "pool_cells": int(pool.n_obs), "dropped_records": int(pool.uns["qc"]["n_dropped"]),
               "n": int(n), "count_record": h_json(count_record)}
        for seed in SEEDS:
            s = tb.sample_cells(pool, n, seed)
            sub, rep = tb.build_submission(board, pool.X[s.rows], s.coords, genes=list(pool.var_names))
            guard = tb.cloud_guard(s.coords, pool.obsm["spatial_3D"], seed=seed)
            res[f"seed{seed}"] = {"rows": h_array(s.rows), "coords": h_array(s.coords),
                                  "sample_record": h_json(s.record), **h_adata(sub), "check": h_json(asdict(rep)),
                                  "check_ok": rep.ok, "guard": h_json(guard.as_dict()), "guard_ok": guard.ok}
            print(f"T2 {key} ({stage}, {pool.n_obs} cells) seed {seed}: n {n}, X {res[f'seed{seed}']['X'][:16]}, "
                  f"obsm {res[f'seed{seed}']['obsm'][:16]}, check ok {rep.ok}, guard ok {guard.ok}", flush=True)
        out[key] = res
    return out


def run(data: Path, out: Path) -> None:
    import numpy as np

    import vec

    pkg = Path(vec.__file__).resolve().parent
    print(f"vec imported from {vec.__file__}", flush=True)
    t0 = time.perf_counter()
    res = {"vec_file": vec.__file__, "host": platform.node(), "python": sys.version.split()[0], "numpy": np.__version__,
           "has_vec_operator": (pkg / "core.py").exists(),
           "builder_files_sha256": {f: sha((pkg / f).read_bytes()) for f in BUILDER_FILES if (pkg / f).exists()},
           "data": str(data), "t1": run_t1(data), "t2": run_t2(data)}
    res["seconds"] = round(time.perf_counter() - t0, 1)
    out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out} ({res['seconds']} s)")


def leaves(d: dict, prefix: str = ""):
    for k, v in d.items():
        if isinstance(v, dict):
            yield from leaves(v, f"{prefix}{k}.")
        else:
            yield f"{prefix}{k}", v


def compare(old: Path, new: Path) -> int:
    a, b = json.loads(old.read_text()), json.loads(new.read_text())
    print(f"old: vec from {a['vec_file']} (vec.core present: {a['has_vec_operator']})")
    print(f"new: vec from {b['vec_file']} (vec.core present: {b['has_vec_operator']})")
    bad = 0
    for part in ("t1", "t2"):
        la, lb = dict(leaves(a[part])), dict(leaves(b[part]))
        for k in sorted(set(la) | set(lb)):
            if k.endswith(".seconds"):
                continue
            same = la.get(k) == lb.get(k)
            bad += not same
            if k.rsplit(".", 1)[-1] in ("X", "obs", "obsm", "var", "layers", "rows", "coords", "report", "check",
                                         "guard", "sample_record") or not same:
                print(f"{'same' if same else 'DIFFERENT'}  {part}.{k}  {str(la.get(k))[:16]}  {str(lb.get(k))[:16]}")
    print("ALL IDENTICAL" if not bad else f"{bad} DIFFERENCES")
    return 1 if bad else 0


def functions(src: str) -> dict[str, tuple[str, str, str]]:
    """name -> (def line, body without the docstring, full source) of each top-level function."""
    lines = src.splitlines()
    out = {}
    for node in ast.parse(src).body:
        if isinstance(node, ast.FunctionDef):
            body = node.body
            if isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                body = body[1:]
            out[node.name] = (lines[node.lineno - 1], "\n".join(lines[body[0].lineno - 1:body[-1].end_lineno]),
                              "\n".join(lines[node.lineno - 1:node.end_lineno]))
    return out


def source(base: str) -> int:
    def at_base(path):
        return subprocess.run(["git", "-C", str(REPO), "show", f"{base}:{path}"], capture_output=True, text=True,
                              check=True).stdout

    core = functions((REPO / "src" / "vec" / "core.py").read_text())
    frozen = functions((REPO / "experiments" / "t3_decomp" / "t3_decompose.py").read_text())
    bad = 0
    for path, old, new in MOVED:
        same = functions(at_base(path))[old][2].replace(f"def {old}(", f"def {new}(", 1) == core[new][2]
        bad += not same
        print(f"{'same' if same else 'DIFFERENT'}  {path}@{base} {old} -> vec.core.{new}: full source, name aside")
    for old, new in FROZEN:
        same = (frozen[old][0].replace(f"def {old}(", f"def {new}(", 1) == core[new][0]
                and frozen[old][1] == core[new][1])
        bad += not same
        print(f"{'same' if same else 'DIFFERENT'}  t3_decompose.{old} -> vec.core.{new}: def line (name aside) and "
              "body; docstrings differ")
    print("ALL IDENTICAL" if not bad else f"{bad} DIFFERENCES")
    return 1 if bad else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="build and hash with the vec package on the path")
    r.add_argument("--data", type=Path, required=True, help="the released data's raw/ directory (read only)")
    r.add_argument("--out", type=Path, required=True)
    c = sub.add_parser("compare", help="compare two run JSONs")
    c.add_argument("old", type=Path)
    c.add_argument("new", type=Path)
    o = sub.add_parser("source", help="compare the core's functions with their sources (git checkout)")
    o.add_argument("--base", default="eb946ac", help="the commit before the move")
    args = ap.parse_args(argv)
    if args.cmd == "run":
        run(args.data, args.out)
        return 0
    if args.cmd == "source":
        return source(args.base)
    return compare(args.old, args.new)


if __name__ == "__main__":
    sys.exit(main())
