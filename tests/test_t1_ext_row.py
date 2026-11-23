"""Row T72's data steps (candidate MX2s), on a synthetic checkout, synthetic records and synthetic files only (no
network, no released or outside data, no git repository touched):

- the guards of cell 8 (row_t72.guards), each refusing on its own: the row in the ledger, the register strings, the
  frozen yaml's pinned sha256 and its version at HEAD, the derived classes' sha256, HEAD at the row's commit (the commit
  that adds the row), a clean tree, the output folders absent and inside the checkout's own data/;
- the three steps end to end (fetch, load, tables): the folders created exclusively (a rerun refuses), the fetched
  barcodes and features files tied to the screened copies, a failed load recorded as a failed (p1) and carried to the
  manifest, the member tables in section B's format with tasks ["T1"].
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t1_ext"
if str(EXT) not in sys.path:
    sys.path.insert(0, str(EXT))
import t1x_common as cm  # noqa: E402


def _load_path(key, path):
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[key] = m
    spec.loader.exec_module(m)
    return m


rt = cm.sibling("row_t72")
fetch = cm.sibling("fetch")
ld = cm.sibling("load")
tb = cm.sibling("tables")
gr = cm.sibling("groups")
H = _load_path("t1x_test_mx2s_helpers", ROOT / "tests" / "test_t1_ext_mx2s.py")  # the synthetic series and files
C40 = "c" * 40
TRIES_ROW = "| T72 | T1 | synthetic |\n"
REGISTER = ("| X002 | Feng; consensus member of T1 candidate MX2s, practice read under row T72 |\n"
            "| X010 | marker statements; skeletal-muscle competing class of row T71 |\n"
            "| X019 | GSE244783; Candidate MX2s, row T72 |\n| X020 | GSE230531; Candidate MX2s, row T72 |\n")


# ------------------------------------------------------------------------------------------------ the guards


def checkout(tmp_path):
    """A synthetic checkout: the ledger with row T72, the register, a frozen copy of atlases.yaml with its derived
    file, and data/. Returns (repo, yaml sha256, git answers)."""
    repo = tmp_path / "checkout"
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "experiments" / "harness" / "TRIES.md").write_text("| Row | Task |\n" + TRIES_ROW)
    (repo / "EXTERNAL_DATA.md").write_text(REGISTER)
    der = repo / "experiments" / "t1_ext" / "results" / "extra_competing.json"
    der.parent.mkdir(parents=True)
    der.write_text('{"stage": "E8.5", "classes": {}}\n')
    raw = yaml.safe_load((EXT / "atlases.yaml").read_text())
    raw["status"] = "frozen"
    raw["extra_competing"]["derived_sha256"] = cm.sha256_file(der)
    yp = repo / rt.ATLASES_REL
    yp.write_text(yaml.safe_dump(raw))
    (repo / "data").mkdir()
    answers = {("show", f"HEAD:{rt.ATLASES_REL}"): yp.read_text(),
               ("show", f"{C40}:{rt.TRIES_REL}"): "| Row | Task |\n" + TRIES_ROW,
               ("show", f"{C40}^:{rt.TRIES_REL}"): "| Row | Task |\n"}
    return repo, cm.sha256_file(yp), answers


def run_guards(repo, sha, answers, step="fetch", commit=C40, state=None):
    return rt.guards(step, commit, repo, git=lambda *a: answers.get(a), state=state or (lambda r: (C40, [])),
                     yaml_sha256=sha)


def test_the_guards_pass_on_a_placed_row(tmp_path):
    repo, sha, ans = checkout(tmp_path)
    g = run_guards(repo, sha, ans)
    assert g == {"row": "T72", "step": "fetch", "commit": C40, "clean": True, "atlases_yaml_sha256": sha,
                 "derived_sha256": cm.sha256_file(repo / "experiments/t1_ext/results/extra_competing.json"),
                 "tasks": ["T1"]}
    assert rt.FROZEN_YAML_SHA256 == "PENDING" or len(rt.FROZEN_YAML_SHA256) == 64  # filled at placement
    assert rt.REGISTER_PHRASES == {"X002": ("consensus member of T1 candidate MX2s", "row T72"),
                                   "X019": ("Candidate MX2s", "row T72"), "X020": ("Candidate MX2s", "row T72"),
                                   "X010": ("skeletal-muscle competing class of row T71",)}


def _edit(path, old, new):
    path.write_text(path.read_text().replace(old, new))


@pytest.mark.parametrize("breakage,why", [
    ("no_row", "is not in"), ("register_x002", "X002"), ("register_x010", "X010"), ("register_x019", "X019"),
    ("yaml_sha", "sha256 the code pins"), ("yaml_pending", "not filled"), ("not_frozen", "sha256 the code pins"),
    ("head_differs", "differs from its version at HEAD"), ("derived", "derived"), ("short_commit", "40-digit"),
    ("head_not_commit", "not the row's commit"), ("unclean", "not clean"), ("unknown_state", "not clean"),
    ("parent_has_row", "not the commit that adds row"), ("folder_exists", "rerun"), ("data_link", "link"),
])
def test_each_guard_refuses(tmp_path, breakage, why):
    repo, sha, ans = checkout(tmp_path)
    kw = {}
    if breakage == "no_row":
        _edit(repo / rt.TRIES_REL, TRIES_ROW, "| T71 | T1 | x |\n")
    elif breakage.startswith("register_"):
        row = breakage.split("_")[1].upper()
        lines = [ln for ln in REGISTER.splitlines(True) if not ln.startswith(f"| {row} |")]
        (repo / "EXTERNAL_DATA.md").write_text("".join(lines) + f"| {row} | something else |\n")
    elif breakage == "yaml_sha":
        sha = "0" * 64
    elif breakage == "yaml_pending":
        sha = "PENDING"
    elif breakage == "not_frozen":
        _edit(repo / rt.ATLASES_REL, "status: frozen", "status: draft")
    elif breakage == "head_differs":
        ans[("show", f"HEAD:{rt.ATLASES_REL}")] += "# another line\n"
    elif breakage == "derived":
        (repo / "experiments/t1_ext/results/extra_competing.json").write_text("{}\n")
    elif breakage == "short_commit":
        kw["commit"] = "c" * 7
    elif breakage == "head_not_commit":
        kw["state"] = lambda r: ("d" * 40, [])
    elif breakage == "unclean":
        kw["state"] = lambda r: (C40, ["experiments/t1_ext/fetch.py"])
    elif breakage == "unknown_state":
        kw["state"] = lambda r: (C40, None)
    elif breakage == "parent_has_row":
        ans[("show", f"{C40}^:{rt.TRIES_REL}")] = "| Row | Task |\n" + TRIES_ROW
    elif breakage == "folder_exists":
        (repo / rt.ROW_DERIVED_REL).mkdir(parents=True)
    elif breakage == "data_link":
        shutil.rmtree(repo / "data")
        (tmp_path / "elsewhere").mkdir()
        (repo / "data").symlink_to(tmp_path / "elsewhere", target_is_directory=True)
    with pytest.raises(cm.Refused, match=why):
        run_guards(repo, sha, ans, **kw)


def test_each_step_refuses_its_own_existing_output(tmp_path):
    repo, sha, ans = checkout(tmp_path)
    ext = repo / rt.ROW_EXTERNAL_REL
    (ext / "reduced").mkdir(parents=True)
    with pytest.raises(cm.Refused, match="rerun"):
        run_guards(repo, sha, ans, "fetch")
    with pytest.raises(cm.Refused, match="rerun"):
        run_guards(repo, sha, ans, "load")
    run_guards(repo, sha, ans, "tables")  # tables needs only data/derived/t1_ext_t72/ absent
    with pytest.raises(cm.Refused, match="rerun"):
        rt.make_exclusive(ext / "reduced")
    with pytest.raises(ValueError):
        run_guards(repo, sha, ans, "pregate")


# ------------------------------------------------------------------------------------------------ the steps


def mtx(n_features, n_barcodes, seed):
    rng = np.random.default_rng(seed)
    ent = [(i + 1, j + 1, int(rng.integers(1, 9))) for i in range(n_features) for j in range(n_barcodes)
           if i < 20]  # counts on the Gene Expression rows only
    head = f"%%MatrixMarket matrix coordinate integer general\n{n_features} {n_barcodes} {len(ent)}\n"
    return head + "".join(f"{i} {j} {v}\n" for i, j, v in ent)


FILES = {**H.FILES, f"{H.EARLY}_matrix.mtx.gz": mtx(25, 12, 1), f"{H.LATE}_matrix.mtx.gz": mtx(29, 12, 2)}


def placed(tmp_path, files=FILES):
    """The screens' fetch into a screens folder, the readers' screen block and pins in the config (as the placement
    commit writes them), and a synthetic checkout for the row's run."""
    fm, cfg, _ = H.meta_ready(tmp_path / "screen", files)
    scr = {g: {r: {"name": f["name"], "sha256": f["sha256"]} for r, f in s["files"].items()}
           for g, s in fm["samples"].items()}
    cfg = H.config(screen={"verdict": "kept", "files": scr, "read_by": ["Reader A", "Reader B"],
                           "date": "2026-11-24"})
    _, cfg2 = H.pinned(tmp_path / "pins", cfg)
    repo = tmp_path / "checkout"
    (repo / "data").mkdir(parents=True)
    return repo, cfg2


def guard(step, commit, repo):
    return {"row": "T72", "step": step, "commit": commit, "clean": True, "atlases_yaml_sha256": "y" * 64,
            "derived_sha256": "d" * 64, "tasks": ["T1"]}


def stub_labeller(monkeypatch):
    monkeypatch.setattr(gr, "label_configs", lambda a, cfg, extra, early, late: [])
    monkeypatch.setattr(gr, "label_atlas", lambda red, runs: {
        "labels": np.array(["r_cm"] * len(red["gsm"])), "skm": np.zeros(len(red["gsm"])), "runs": []})


def test_the_three_steps_end_to_end(tmp_path, monkeypatch):
    repo, cfg = placed(tmp_path)
    log = []
    out = fetch.row_fetch(C40, H.getter_for(H.soft(H.SAMPLES), log, FILES), repo, guard, cfg, members=(H.ACC,))
    ext = repo / rt.ROW_EXTERNAL_REL
    assert (ext / "FETCH.json").is_file() and (ext / "STARTED.json").is_file()
    assert set(out["atlases"][H.ACC]["files"]) == {H.EARLY, H.LATE} and out["tasks"] == ["T1"]
    data = [u for u in log if not u.endswith("_family.soft.gz")]
    assert len(data) == 6 and all(H.EARLY in u or H.LATE in u for u in data)  # the pinned samples' three files
    assert not any(n in u for u in log for n in cm.NEVER_FETCHED_ALL)
    with pytest.raises(cm.Refused, match="rerun"):  # the fetch folder blocks a rerun
        fetch.row_fetch(C40, H.getter_for(H.soft(H.SAMPLES), [], FILES), repo, guard, cfg, members=(H.ACC,))
    lj = ld.row_load(C40, repo, guard, cfg)
    assert lj["p1"] == {H.ACC: {"passed": True, "rechecked": {"suffix": True, "own_peaks": True}}}
    red = ld.read_reduced(H.ACC, ext, cfg)
    assert np.all(red["total_ge"] > 0) and red["C"].shape[0] == 24
    with pytest.raises(cm.Refused, match="rerun"):
        ld.row_load(C40, repo, guard, cfg)
    stub_labeller(monkeypatch)
    t1 = [f"Gene{i}" for i in range(20)] + ["mt-Co1"]
    man = tb.row_tables(C40, tmp_path, repo, guard, cfg, t1_genes=t1, extra={})
    der = repo / rt.ROW_DERIVED_REL
    e = man["atlases"][H.ACC]["groups"]["CM"]
    assert e["table"] == f"tables/{H.ACC}/CM.npz" and cm.sha256_file(der / e["table"]) == e["sha256"]
    z = np.load(der / e["table"])
    assert sorted(z.files) == sorted(tb.TABLE_KEYS) and list(z["genes"]) == t1
    assert np.isfinite(z["D"][:20]).all() and np.isnan(z["D"][20])
    assert e["cells_per_end"] == {"early": 12, "late": 12}
    j = json.loads((der / "MANIFEST.json").read_text())
    assert j["tasks"] == ["T1"] and j["commit"] == C40 and j["p1"][H.ACC]["passed"]
    assert j["atlases"][H.ACC]["skm_purity_late"]["late_r_cm_cells"] == 12
    assert j["constants"]["pseudocount"] == 0.5 and j["atlases"][H.ACC]["factor"] == 5.0 / 8.0
    text = json.dumps(j)
    assert "composition" not in text and "consensus" not in text
    with pytest.raises(cm.Refused, match="rerun"):
        tb.row_tables(C40, tmp_path, repo, guard, cfg, t1_genes=t1, extra={})


def test_a_fetched_file_that_is_not_the_screened_copy_is_a_failed_p1(tmp_path, monkeypatch):
    repo, cfg = placed(tmp_path)
    other = {**FILES, f"{H.LATE}_barcodes.tsv.gz": H.BARCODES + "TTTTGGGG-1\n"}  # the repository changed the file
    fetch.row_fetch(C40, H.getter_for(H.soft(H.SAMPLES), [], other), repo, guard, cfg, members=(H.ACC,))
    lj = ld.row_load(C40, repo, guard, cfg)
    assert lj["p1"][H.ACC]["passed"] is False
    assert not (repo / rt.ROW_EXTERNAL_REL / "reduced" / H.ACC).exists()
    stub_labeller(monkeypatch)
    man = tb.row_tables(C40, tmp_path, repo, guard, cfg, t1_genes=["Gene0"], extra={})
    assert man["atlases"][H.ACC] == {"p1": lj["p1"][H.ACC], "groups": {}}


def test_a_check_that_fails_again_on_the_fetched_files_is_a_failed_p1(tmp_path):
    repo, cfg = placed(tmp_path)
    fetch.row_fetch(C40, H.getter_for(H.soft(H.SAMPLES), [], FILES), repo, guard, cfg, members=(H.ACC,))
    # the late features file on disk replaced by the early one's bytes (as if both libraries shared one peak set):
    # its sha256 no longer matches the manifest, so the load refuses before any check reads it
    raw = repo / rt.ROW_EXTERNAL_REL / "raw" / H.ACC
    (raw / H.LATE / f"{H.LATE}_features.tsv.gz").write_bytes(
        (raw / H.EARLY / f"{H.EARLY}_features.tsv.gz").read_bytes())
    lj = ld.row_load(C40, repo, guard, cfg)
    assert lj["p1"][H.ACC]["passed"] is False and "not the file the fetch manifest records" in \
        lj["p1"][H.ACC]["reason"]


def test_recheck_screens_says_no_on_one_peak_set(tmp_path):
    a = cm.atlas(H.ACC, H.config())
    d = tmp_path / "raw"
    man = {"samples": {}}
    for g in (H.EARLY, H.LATE):
        (d / g).mkdir(parents=True)
        for r, body in (("barcodes", H.BARCODES), ("features", H.features(range(0, 5000, 1000)))):
            p = d / g / f"{g}_{r}.tsv.gz"
            p.write_bytes(gzip.compress(body.encode(), mtime=0))
            man["samples"].setdefault(g, {"files": {}})["files"][r] = {"name": p.name, "sha256": cm.sha256_file(p)}
    with pytest.raises(cm.Refused, match="own_peaks"):
        ld.recheck_screens(H.ACC, a, cm.load_atlases(), man, d)


def test_a_dropped_member_is_never_requested(tmp_path):
    repo, cfg = placed(tmp_path)
    cfg["atlases"][H.ACC]["status"] = "dropped"
    cfg["atlases"][H.ACC]["screen"] = {"verdict": "dropped", "failed": ["b"]}
    log = []
    with pytest.raises(cm.Refused, match="no member is kept"):
        fetch.row_fetch(C40, H.getter_for(H.soft(H.SAMPLES), log, FILES), repo, guard, cfg, members=(H.ACC,))
    assert log == []
