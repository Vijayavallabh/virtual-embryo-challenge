"""Row R72's sign cube builder (experiments/t3_decomp/build_sign_cube.py) on synthetic arrays and hand-made ledgers:
its refusals (form, vertex, seed, count, ruling, row R72's whole option strings, row K31's manifest sha256), the
vertex factor from the tables (ppp C2's own array, the others held), the loader hook, the vertex list read from row
K31's reading against G0's, and the command line's refusals before anything loads. Nothing is built and no
competition file is read."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "t3_decomp"))
import build_sign_cube as bsc  # noqa: E402
import sign_cube as sc  # noqa: E402

VERTS = sc.NAMED
SHA = "a" * 64


def _r72(verts=VERTS, ruling="`--captain-ruling R72`") -> str:
    return ("| R72 (the captain's order of 2026-11-25, ruling, no run: \"Yes, ... each only with " + ruling + " ...: "
            + ", ".join(f"T3_cube_{v}.h5ad `--form c2-s2 --vertex {v}`" for v in verts) + "\") | T3 | x |\n")


def _k31(sha=SHA) -> str:
    return f"| K31 (one row `experiments/t3_decomp/k31_cube_read.py`; cube tables manifest sha256 `{sha}`) | T3 | x |\n"


def _ledger(tmp_path, *rows) -> Path:
    p = tmp_path / "TRIES.md"
    p.write_text("| K30 (x) | T3 | y |\n" + "".join(rows))
    return p


@pytest.mark.parametrize("args,why", [
    (("c2-s1", "ppn", "R72"), "--form"), (("c2-s2 ", "ppn", "R72"), "--form"), (("c2-s2", "pp", "R72"), "--vertex"),
    (("c2-s2", "PPN", "R72"), "--vertex"), (("c2-s2", "ppn", None), "captain-ruling"),
    (("c2-s2", "ppn", "R22"), "captain-ruling"),
])
def test_the_command_is_refused_unless_c2s2_a_vertex_and_r72(tmp_path, args, why):
    with pytest.raises(bsc.Refused, match=why):
        bsc.check_command(*args, tries=_ledger(tmp_path, _r72()))


def test_seed_and_count_are_c2s():
    with pytest.raises(bsc.Refused, match="1000001"):
        bsc.check_command("c2-s2", "ppn", "R72", seed=7)
    with pytest.raises(bsc.Refused, match="max_cells"):
        bsc.check_command("c2-s2", "ppn", "R72", n_cells="5000")


def test_row_r72_must_name_the_vertex_as_one_whole_string(tmp_path):
    rec = bsc.check_command("c2-s2", "pnn", "R72", tries=_ledger(tmp_path, _r72()))
    assert rec["option"] == "`--form c2-s2 --vertex pnn`" and rec["file"] == "T3_cube_pnn.h5ad"
    with pytest.raises(bsc.Refused, match="one whole string"):
        bsc.check_command("c2-s2", "pnn", "R72", tries=_ledger(tmp_path, _r72(verts=("ppn", "npp"))))
    bad = _r72().replace("`--form c2-s2 --vertex pnn`", "`--form c2-s2 --vertex pnn --size-scale 0.67`")
    with pytest.raises(bsc.Refused, match="one whole string"):
        bsc.check_command("c2-s2", "pnn", "R72", tries=_ledger(tmp_path, bad))
    for led in (_r72(ruling="--captain-ruling R72"), _r72().replace("the captain's order", "a note"),
                _r72().replace("ruling, no run", "ruled"), _r72() * 2, ""):
        with pytest.raises(bsc.Refused, match="captain's order"):
            bsc.check_command("c2-s2", "pnn", "R72", tries=_ledger(tmp_path, led))


def test_ppp_the_identity_gate_needs_the_order_only(tmp_path):
    rec = bsc.check_command("c2-s2", "ppp", "R72", tries=_ledger(tmp_path, _r72()))
    assert rec["vertex"] == "ppp" and "never uploaded" in rec["file"]
    with pytest.raises(bsc.Refused):
        bsc.check_command("c2-s2", "ppp", None, tries=_ledger(tmp_path, _r72()))


def test_the_manifest_must_be_the_one_row_k31_records(tmp_path):
    man = tmp_path / "MANIFEST.json"
    man.write_text('{"x": 1}\n')
    sha = hashlib.sha256(man.read_bytes()).hexdigest()
    assert bsc.check_manifest(man, _ledger(tmp_path, _r72(), _k31(sha))) == sha
    with pytest.raises(bsc.Refused, match="K31"):
        bsc.check_manifest(man, _ledger(tmp_path, _r72(), _k31("b" * 64)))
    with pytest.raises(Exception):  # noqa: B017  (no row K31: k31_cube_read refuses)
        bsc.check_manifest(man, _ledger(tmp_path, _r72()))


def _synthetic(seed=0, G=300, L=8, N=2000):
    rng = np.random.default_rng(seed)
    lab = rng.integers(0, L, N)
    prof = rng.dirichlet(np.ones(G) * 0.3, size=L)
    C = np.vstack([rng.poisson(prof[t] * 150) for t in lab])
    l_m, l_g = rng.normal(0, 0.5, G), rng.normal(0, 0.4, G)
    tabs = {"S_M": np.exp(0.5 * l_m + rng.normal(0, 0.05, G)), "Y_M": np.exp(-0.5 * l_m + rng.normal(0, 0.05, G)),
            "S_G": np.exp(0.3 * l_g + rng.normal(0, 0.05, G)), "Y_G": np.exp(-0.3 * l_g + rng.normal(0, 0.05, G))}
    clip = np.log(4)
    return (np.exp(np.clip(l_m, -clip, clip)), np.exp(np.clip(l_g, -clip, clip)), tabs, C, lab, C.sum(0) / C.sum())


def test_ppp_is_c2s_combined_factor_exactly():
    import build_two_knockout as btk

    phi_m, phi_g, tabs, C, lab, den = _synthetic()
    phi, rec, phi_c2 = bsc.vertex_of(phi_m, phi_g, tabs, C, lab, den, "ppp")
    want, _ = btk.combine("c2", 1.0, phi_m, phi_g)
    assert phi.tobytes() == want.tobytes() == phi_c2.tobytes() and rec["identity"] and rec["k"] == 1.0


@pytest.mark.parametrize("v", VERTS)
def test_a_named_vertex_is_held_within_the_band(v):
    phi_m, phi_g, tabs, C, lab, den = _synthetic(1)
    phi, rec, phi_c2 = bsc.vertex_of(phi_m, phi_g, tabs, C, lab, den, v)
    assert not rec["identity"] and sc.within_band(rec) and rec["hold"]["converged"]
    assert abs(den @ phi / (den @ phi_c2) - 1) <= sc.DEN_TOL


def test_the_loader_hook_hands_the_checked_factor_once(monkeypatch):
    import build_two_knockout as btk

    calls = []

    def frozen(run_dir, name, genes, dec):
        calls.append(name)
        return np.ones(len(genes)), {"file": f"{name}.csv", "sha256": btk.TABLE_SHA256.get(name, "")}

    monkeypatch.setattr(bsc.bs, "load_factor", frozen)
    phi, out = np.array([0.5, 2.0, 1.0]), {}
    with bsc.vertex_loader(phi, ["a", "b", "c"], out):
        got, rec = bsc.bs.load_factor(None, btk.TABLE_M, ["a", "b", "c"], {})
        assert np.array_equal(got, phi) and out["n_loads"] == 1 and out["factor_sha256"] == sc.factor_sha256(phi)
        with pytest.raises(bsc.FactorMismatch):
            bsc.bs.load_factor(None, btk.TABLE_M, ["a", "c", "b"], {})
        with pytest.raises(bsc.Refused):
            bsc.bs.load_factor(None, btk.TABLE_G, ["a", "b", "c"], {})
    assert bsc.bs.load_factor is frozen and calls == [btk.TABLE_M]


def _readings(tmp_path, outcome="CONTINUE", build=("ppn", "npp"), left=VERTS, g_outcome="CONTINUE", sha=None):
    g0 = tmp_path / "g0.json"
    g0.write_text(json.dumps({"outcome": g_outcome, "vertices_left": list(left)}))
    r = tmp_path / "reading.json"
    r.write_text(json.dumps({"outcome": outcome, "build": list(build),
                             "g0": {"sha256": sha or hashlib.sha256(g0.read_bytes()).hexdigest()}}))
    return r, g0


def test_the_vertices_built_are_the_readings_checked_against_g0s(tmp_path):
    assert bsc.vertices_to_build(*_readings(tmp_path)) == ["ppn", "npp"]
    assert bsc.vertices_to_build(*_readings(tmp_path, outcome="STOP", build=())) == []
    for kw in ({"sha": "b" * 64}, {"left": ("ppn",)}, {"build": ("ppp",)}, {"build": ("npp", "npp")},
               {"g_outcome": "STOP"}):
        with pytest.raises(bsc.Refused):
            bsc.vertices_to_build(*_readings(tmp_path, **kw))


def _cli(*args, cwd=REPO) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-I", "-B", str(REPO / "experiments/t3_decomp/build_sign_cube.py"), *args],
                          cwd=cwd, capture_output=True, text=True, timeout=300)


@pytest.mark.parametrize("change", [("--form", "c2-s1"), ("--vertex", "ppx"), ("--captain-ruling", "R22"),
                                    ("--seed", "7"), ("--n-cells", "100")])
def test_the_command_line_refuses_before_anything_loads(tmp_path, change):
    base = {"--form": "c2-s2", "--vertex": "ppn", "--wt": str(tmp_path / "none.h5ad"), "--board": "T3:bcatenin",
            "--out": str(tmp_path / "out.h5ad"), "--cube-tables": str(tmp_path / "none.json"),
            "--captain-ruling": "R72"}
    base[change[0]] = change[1]
    args = ["build", str(tmp_path / "dec.json")] + [x for kv in base.items() for x in kv]
    r = _cli(*args)
    assert r.returncode == 2 and "refused" in (r.stdout + r.stderr)
    assert not (tmp_path / "out.h5ad").exists()


def test_the_command_line_lists_vertices_or_says_none(tmp_path):
    r, g0 = _readings(tmp_path)
    got = _cli("vertices", "--reading", str(r), "--g0", str(g0))
    assert got.returncode == 0 and got.stdout.split() == ["ppn", "npp"]
    r, g0 = _readings(tmp_path, outcome="VOID", build=())
    assert _cli("vertices", "--reading", str(r), "--g0", str(g0)).returncode == 3
    r, g0 = _readings(tmp_path, sha="c" * 64)
    assert _cli("vertices", "--reading", str(r), "--g0", str(g0)).returncode == 2


def test_the_builder_uses_the_cube_code_that_made_the_tables():
    import cube_tables as ct

    man, _ = ct.load_manifest(ct.DEFAULT_OUT / ct.MANIFEST)
    got = bsc.check_code(man)  # the committed code is the manifest's
    assert got["sign_cube.py"] == man["code_sha256"]["sign_cube.py"] and got["cube_tables.py"] == man["script_sha256"]
    with pytest.raises(bsc.Refused, match="not the code"):
        bsc.check_code({**man, "code_sha256": {**man["code_sha256"], "sign_cube.py": "0" * 64}})
    with pytest.raises(bsc.Refused, match="not the code"):
        bsc.check_code({**man, "script_sha256": "0" * 64})
