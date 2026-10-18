import hashlib
import json

import pytest

from vec.boards import load_boards, sanitise_index, update_boards


def _spec(key, split, anchors=True):
    s = {"key": key, "task": key.split(":")[0], "split": split, "genes_file": key.replace(":", "__") + ".genes.txt",
         "needs_coords": False, "min_cells": 1000, "max_cells": 5000, "truth_cells": 1500, "ref_cells": 1700,
         "floor_model": "copy_last", "secret_extra": 1}
    if anchors:
        s["anchors"] = {"mmd_u": {"floor": 0.1, "ceiling": 0.01}}
    return s


def test_sanitise_strips_everything_but_the_spec_from_non_validation_boards():
    index = {"T1:val": _spec("T1:val", "val"), "T1:test": _spec("T1:test", "test"),
             "T3:gata4": _spec("T3:gata4", "test"), "T2:heart:new": _spec("T2:heart:new", "val_extrap")}
    clean, reduced = sanitise_index(index)
    assert clean["T1:val"] == index["T1:val"]
    assert sorted(reduced) == ["T1:test", "T2:heart:new", "T3:gata4"]
    for k in reduced:
        assert "anchors" not in clean[k] and "secret_extra" not in clean[k]
        assert clean[k]["truth_cells"] == 1500 and clean[k]["max_cells"] == 5000


def test_update_boards_writes_no_test_anchors_and_quarantines_unreleased_composition(tmp_path):
    genes = b"Gata4\nNkx2-5\n"
    index = {"T1:val": _spec("T1:val", "val"), "T1:test": _spec("T1:test", "test")}
    for s in index.values():
        s["genes_sha256"] = hashlib.sha256(genes).hexdigest()
    comp = {"stages": {"E8.5": {}, "E9.5": {}, "E12.5": {"x": 1}}}
    files = {"index.json": json.dumps(index).encode(), "t1_composition.json": json.dumps(comp).encode(),
             "T1__val.genes.txt": genes, "T1__test.genes.txt": genes}
    summary = update_boards(tmp_path, fetch=lambda url: files[url.rsplit("/", 1)[1]])
    written = json.loads((tmp_path / "index.json").read_text())
    assert "anchors" not in written["T1:test"] and written["T1:val"]["anchors"]
    assert summary["reduced"] == ["T1:test"]
    assert not (tmp_path / "t1_composition.json").exists() and summary["t1_composition_unreleased_stages"] == 1
    boards = load_boards(tmp_path)
    assert boards["T1:test"].anchors == {} and boards["T1:test"].truth_cells == 1500


def test_update_boards_rejects_a_bad_gene_hash(tmp_path):
    index = {"T1:val": _spec("T1:val", "val")}
    index["T1:val"]["genes_sha256"] = "0" * 64
    files = {"index.json": json.dumps(index).encode(), "T1__val.genes.txt": b"A\n"}
    with pytest.raises(ValueError, match="sha256"):
        update_boards(tmp_path, fetch=lambda url: files[url.rsplit("/", 1)[1]])
