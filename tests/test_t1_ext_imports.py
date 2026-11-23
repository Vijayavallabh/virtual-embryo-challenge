"""Row T72 (candidate MX2s), engineer_notes A.11: the outside-atlas tables are Task 1 only. E8.5 lies inside the T2
heart's interpolation window, so no module under src/vec/ (the shared core and the T1 and T2 heads), no T2 or T3 path
(experiments/t2*, experiments/t3*) and no builder (experiments/builds/, build_test.py among them) imports
experiments/t1_ext or reads data/derived/t1_ext_t72 or data/external/t1_ext_t72. The T1 practice read
(experiments/t1_m1) may: it is the row's own reader. Source text only; no module is imported and no data is read.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# A path into the folder or its module names, or the row's two data folders.
NAMES = re.compile(r"experiments[/.]t1_ext\b|(?<![\w])t1_ext[/.]|[\"']t1_ext[\"']|t1x_common|row_t72|"
                   r"t1_ext_t72|t1_ext_mod_")


def scanned() -> list[Path]:
    roots = [ROOT / "src" / "vec", ROOT / "experiments" / "builds"]
    roots += sorted(p for p in (ROOT / "experiments").glob("t2*") if p.is_dir())
    roots += sorted(p for p in (ROOT / "experiments").glob("t3*") if p.is_dir())
    out = []
    for r in roots:
        for p in [*r.rglob("*.py"), *r.rglob("*.sh"), *r.rglob("*.yaml")]:
            if "results" in p.relative_to(ROOT).parts:
                continue  # results folders hold outputs, some never opened
            out.append(p)
    return out


def test_the_pattern_names_every_route_in():
    for ref in ('EXT = ROOT / "experiments" / "t1_ext"', "experiments/t1_ext/tables.py", "import t1x_common",
                "from experiments.t1_ext import tables", "sys.path.insert(0, str(HERE / 't1_ext'))",
                "data/derived/t1_ext_t72/tables/GSE244783/CM.npz", "import row_t72",
                "spec_from_file_location('t1_ext_mod_tables', p)"):
        assert NAMES.search(ref), ref
    assert not NAMES.search("tests/test_t1_ext_mx2s.py tests/test_t1_ext_load.py")  # a test file's name is not a use


def test_no_shared_core_t2_t3_or_builder_path_imports_the_tables():
    files = scanned()
    assert any(p.name == "build_test.py" for p in files)  # the builder is scanned
    assert any("src/vec" in p.as_posix() for p in files)
    hits = [str(p.relative_to(ROOT)) for p in files if NAMES.search(p.read_text(errors="replace"))]
    assert hits == []
    assert "t1_ext" not in (ROOT / "experiments" / "builds" / "build_test.py").read_text()
